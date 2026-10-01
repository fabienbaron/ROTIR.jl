# spectral_cube.jl
# ---------------------------------------------------------------------------
# Velocity-resolved observables: one surface, one intensity provider, and both the
# interferometric visibilities and the disk-integrated line profile out of the SAME sum.
#
# WHAT WAS MISSING. ROTIR's forward model collapses the wavelength axis at load — every demo
# does `data = data_all[1, :]` — and `grep polychromatic src/ demos/ test/` returned nothing.
# There is one geometry per EPOCH and no wavelength axis anywhere in `setup_oi!` or any χ².
# So there were no spectrally-resolved observables at all: no differential visibility across
# a line, no line profile computed from the same surface the visibilities come from.
#
# WHY NO KERNEL HAD TO CHANGE. The natural worry is that a λ-dependent intensity makes the
# per-tessel weight vector `xw` an `(ntessel, nλ)` object, which nothing accepts —
# `setup_polyft_single`, `interferometric_chi2` and `compute_polyflux_and_cvis!` all take
# `xw::Vector`. But the GEOMETRY does not depend on λ. Only `xw` and the uv scaling do. So a
# loop over channels, each iteration calling the existing scalar-`xw` kernel, is exactly
# correct and reuses every backend (`:t3`, `:nufft`, `:turbo`, `:scalar`) and every rrule
# verbatim. `demos/spica_proximity_chara_2027.jl:264-268` already does this by hand, and its
# comment states the constraint correctly. Widening the kernels is an optimisation to be
# justified by profiling, not a prerequisite.
#
# The one trick is that `fused_cvis_parts(x, star, data; intensity_model = :linear)` means
# "use `x` directly as surface brightness" — so passing the per-tessel intensity a provider
# already computed goes straight through it. With `ldtype = 0` the `ldmap` factor it applies
# is unity, so the provider's μ dependence is not counted twice.
#
# THE SIGN. `los_velocity` is positive RECEDING. A receding tessel's spectrum is REDSHIFTED,
# so at observed channel λ_c that tessel contributes its rest-frame wavelength
# λ_rest = λ_c / (1 + v/c), which is BLUEWARD of λ_c. Getting this backwards mirrors the line
# profile about its centre — which for a symmetric rotator looks perfectly plausible, and only
# shows up once the surface is not symmetric (a spot, or a differentially rotating star).
# ---------------------------------------------------------------------------

"""
    rest_lambda(λ_obs, v_rad) -> λ_rest

Rest-frame wavelength that a tessel with receding-positive velocity `v_rad` (km/s)
contributes at the observed wavelength `λ_obs`:

    λ_rest = λ_obs / (1 + v_rad/c)

The INVERSE of [`doppler_lambda`](@ref), and the one a spectral cube needs: the channel is
fixed by the instrument and the question is which part of each tessel's rest-frame spectrum
lands in it. Broadcasts, so `λ_obs` may be a scalar and `v_rad` a per-tessel vector.
"""
@inline function rest_lambda(λ_obs, v_rad)
    # Typed constant, not the bare Float64 `_C_KMS` — see the note in `doppler_lambda`. This
    # one matters more: the spectral cube calls it once per channel, so a Float32 model would
    # allocate a Float64 λ vector per channel and query the grid in double precision.
    T = float(promote_type(eltype(λ_obs), eltype(v_rad)))
    return λ_obs ./ (one(T) .+ v_rad ./ T(_C_KMS))
end

"""
    SurfaceState{T}

Everything per-tessel that a spectral cube needs and that does NOT depend on wavelength:
effective temperature, `logg`, the limb cosine μ, and the line-of-sight velocity.

Built once per epoch by [`surface_state`](@ref) and reused across every channel — which is
the whole reason it exists, since `logg_map` and `los_velocity` would otherwise be recomputed
tens of times per epoch for no change in their answer.
"""
struct SurfaceState{T}
    Teff::Vector{T}
    logg::Vector{T}
    μ::Vector{T}
    v_los::Vector{T}
    polyflux::Vector{T}          # projected areas of the VISIBLE quads (shoelace)
    visible::Vector{Int}         # star.index_quads_visible
    weights::Vector{T}           # vis_weights .* ldmap, on the visible quads
end

"""
    surface_state(star, star_params; temperature_map = nothing) -> SurfaceState

Assemble the wavelength-independent per-tessel state of one epoch.

`temperature_map` defaults to `parametric_temperature_map(star_params, star)`, i.e. the
gravity-darkened map under whichever law `star_params` names. Pass one explicitly to use a
reconstructed or spotted map instead.

Requires a distance: both `logg` and the velocity are meaningless without one (see
`src/stellar_physics.jl`). Warns if `ldtype ≠ 0`, because an atmosphere provider and a
limb-darkening law together count the limb twice.
"""
function surface_state(star, star_params; temperature_map = nothing)
    has_physical_scale(star_params) || throw(ArgumentError(
        "surface_state: `star_params` needs a distance `d` — both `logg` and the velocity " *
        "field are undefined without one. See src/stellar_physics.jl."))
    T = eltype(star.proj_west)
    Tmap = temperature_map === nothing ?
           parametric_temperature_map(star_params, star) : temperature_map
    # Colatitude comes from the BODY frame, which is parameter-independent.
    θ = @view star.vertices_spherical[:, 5, 2]
    lg = logg_map(T(star_params.rpole), T(star_params.d), T(star_params.frac_escapevel),
                  T(star_params.rotation_period), T.(sin.(θ)), T.(cos.(θ)))
    μ = limb_mu.(T.(star.normals[:, 3]))
    v = T.(los_velocity(star, star_params))
    indx = star.index_quads_visible
    pf = setup_polyflux_single(@view(star.proj_west[indx, :]),
                              @view(star.proj_north[indx, :]))
    w = T.(star.vis_weights[indx] .* star.ldmap[indx])
    return SurfaceState{T}(T.(Tmap), lg, μ, v, pf, indx, w)
end

"""
    channel_intensity(provider, state, λ_obs) -> Vector

Per-tessel intensity at one observed wavelength, with each tessel evaluated at its own
Doppler-shifted rest wavelength.

This is where the velocity field enters the observables. With a provider that ignores λ (or
a zero velocity field) it reduces to a single grid lookup per tessel.
"""
channel_intensity(p::IntensityProvider, st::SurfaceState, λ_obs) =
    provider_intensity(p, st.Teff, st.logg, st.μ, rest_lambda(λ_obs, st.v_los))

"""
    line_profile(provider, state, λ_channels) -> Vector

Disk-integrated flux per channel — the spectrum the star would show, from the same surface
and the same intensities the visibilities are built from.

**No Fourier transform is involved**, so this is two orders of magnitude cheaper than
[`spectral_cvis`](@ref) and is the path a purely spectroscopic fit (Doppler imaging) wants.
The flux is `Σ I_i · w_i · dA_i` over visible tessels, with `w` the soft-visibility (times
`ldmap`, unity under `ldtype = 0`) weights and `dA` the projected areas — exactly the
normalisation `fused_cvis_parts` divides its transform by, so the profile and the
visibilities cannot disagree about the star.

Returns absolute flux in the provider's units. Divide by a continuum channel for a
normalised profile; see [`normalize_profile`](@ref).
"""
function line_profile(p::IntensityProvider, st::SurfaceState{T},
                      λ_channels::AbstractVector) where {T}
    F = Vector{T}(undef, length(λ_channels))
    @inbounds for (c, λc) in enumerate(λ_channels)
        I = channel_intensity(p, st, λc)
        F[c] = dot(st.polyflux, T.(I[st.visible]) .* st.weights)
    end
    return F
end

"""
    spectral_cvis(provider, state, star, data_channels; λ = nothing)
        -> (cvis, flux, λ)

Complex visibilities and total flux per wavelength channel, for one epoch.

`data_channels` is a vector of `OIdata` for the SAME epoch at different wavelengths — i.e. a
column of what `readoifits(..., polychromatic = true)` returns, which every demo currently
throws away with `data = data_all[1, :]`. Each channel carries its own `uv`, because `uv` is
`B/λ` and so is channel-specific; that is why the transform has to be redone per channel
while the geometry is not.

`cvis[c]` is normalised (`F/flux`) exactly as `fused_cvis` is, so it is directly
comparable to `poly_to_cvis` on a single channel. `flux[c]` is the same quantity
[`line_profile`](@ref) returns — computed here as a by-product of the normalisation rather
than a second time.
"""
function spectral_cvis(p::IntensityProvider, st::SurfaceState{T}, star,
                       data_channels::AbstractVector; λ = nothing) where {T}
    nc = length(data_channels)
    λs = λ === nothing ? T[T(band_of(d)) for d in data_channels] : T.(λ)
    length(λs) == nc || throw(DimensionMismatch(
        "spectral_cvis: $(length(λs)) wavelengths for $nc data channels"))
    cvis = Vector{Vector{Complex{T}}}(undef, nc)
    flux = Vector{T}(undef, nc)
    @inbounds for c in 1:nc
        I = channel_intensity(p, st, λs[c])
        # `intensity_model = :linear` means "use `I` as the surface brightness", which is
        # precisely what a provider has already produced. `fused_cvis_parts` then applies
        # `vis_weights .* ldmap` itself — and `ldmap` is unity under ldtype = 0, so the
        # provider's μ dependence is not applied twice.
        F, fl = fused_cvis_parts(I, star, data_channels[c]; intensity_model = :linear)
        cvis[c] = F ./ fl
        flux[c] = fl
    end
    return cvis, flux, λs
end

# ===========================================================================
# Differential observables
# ===========================================================================
# OIFITS v2 already carries these: `visamp`/`visphi` tagged `amptyp`/`phityp` =
# "differential", which OITOOLS reads and stores. They sit at weight positions 4 and 5 of the
# 7-element weight vector and are currently switched OFF —
# `OI_DEFAULT_WEIGHTS = [1,1,1,0,0,0,0]` (src/oichi2_spheroid.jl:211). Turning them on is a
# weight change, not new plumbing.

"""
    continuum_mask(λ, window) -> BitVector

`true` for channels OUTSIDE `window = (λ_lo, λ_hi)`, i.e. the continuum reference channels.

Errors if fewer than two channels fall outside: a differential quantity measured against a
single channel carries that channel's noise into every point of the line, and against none
is undefined.
"""
function continuum_mask(λ::AbstractVector, window::Tuple{<:Real,<:Real})
    m = .!(window[1] .<= λ .<= window[2])
    n = count(m)
    n >= 2 || throw(ArgumentError(
        "continuum_mask: only $n channel(s) fall outside the window " *
        "$(window) — a differential observable needs at least 2 continuum channels. " *
        "Widen the wavelength coverage or narrow the window."))
    return m
end

"""
    differential_observables(cvis, λ, window) -> (dvisamp, dvisphi)

Differential visibility amplitude and phase across a line, relative to the continuum.

    dvisamp[c] = |V_c| / |V_cont|
    dvisphi[c] = arg(V_c) − arg(V_cont)      [degrees, wrapped to (−180, 180]]

with `V_cont` the mean of the COMPLEX visibility over the continuum channels
([`continuum_mask`](@ref)). Averaging the complex visibility rather than the phases is what
avoids a wrap-around bias when the continuum phase sits near ±180°; it is also what an
instrument pipeline does.

Each returned array is `(nuv, nchannel)`. Requires every channel to share a `uv` layout —
true for a polychromatic split of one exposure, and checked.
"""
function differential_observables(cvis::AbstractVector{<:AbstractVector{<:Complex}},
                                  λ::AbstractVector, window::Tuple{<:Real,<:Real})
    nc = length(cvis)
    nuv = length(first(cvis))
    all(length(v) == nuv for v in cvis) || throw(DimensionMismatch(
        "differential_observables: channels have different uv counts " *
        "$(unique(length.(cvis))); a differential quantity compares the SAME baseline " *
        "across wavelength, so the channels must share a uv layout"))
    m = continuum_mask(λ, window)
    T = real(eltype(first(cvis)))
    ref = Vector{Complex{T}}(undef, nuv)
    @inbounds for k in 1:nuv
        s = zero(Complex{T})
        for c in 1:nc
            m[c] && (s += cvis[c][k])
        end
        ref[k] = s / count(m)
    end
    dvisamp = Matrix{T}(undef, nuv, nc)
    dvisphi = Matrix{T}(undef, nuv, nc)
    @inbounds for c in 1:nc, k in 1:nuv
        dvisamp[k, c] = abs(cvis[c][k]) / abs(ref[k])
        dvisphi[k, c] = mod360(rad2deg(angle(cvis[c][k]) - angle(ref[k])))
    end
    return dvisamp, dvisphi
end

"""
    normalize_profile(F, λ, window) -> Vector

`F` divided by its mean over the continuum channels — the normalised line profile a
spectroscopic dataset is compared against.

Uses the same [`continuum_mask`](@ref) as the differential visibilities, so the two are
normalised against the same reference and a residual in one is comparable to a residual in
the other.
"""
function normalize_profile(F::AbstractVector, λ::AbstractVector,
                           window::Tuple{<:Real,<:Real})
    m = continuum_mask(λ, window)
    return F ./ (sum(@view F[m]) / count(m))
end

"""
    line_equivalent_width(F, λ, window) -> W

Equivalent width of the line in `F`, in the units of `λ`.

    W = Σ (1 − F_c/F_cont) Δλ_c

over the channels inside `window`, with `Δλ` from the channel spacing. Positive for
absorption. A cheap scalar summary for checking that a profile behaves as the geometry says
it should — it must, for instance, be independent of inclination for a rigidly rotating star
with no spots, because rotation redistributes flux in wavelength without removing any.
"""
function line_equivalent_width(F::AbstractVector, λ::AbstractVector,
                               window::Tuple{<:Real,<:Real})
    R = normalize_profile(F, λ, window)
    inw = findall(window[1] .<= λ .<= window[2])
    isempty(inw) && return zero(eltype(R))
    W = zero(eltype(R))
    for i in inw
        # Central difference on the channel grid; half-width at the ends.
        dλ = if i == firstindex(λ)
            λ[i+1] - λ[i]
        elseif i == lastindex(λ)
            λ[i] - λ[i-1]
        else
            (λ[i+1] - λ[i-1]) / 2
        end
        W += (one(eltype(R)) - R[i]) * dλ
    end
    return W
end

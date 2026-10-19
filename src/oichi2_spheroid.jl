using LinearAlgebra, OptimPackNextGen
# The Fourier transform for polygonal surfaces
function poly_to_cvis(x, polyflux, polyft)
  flux = sum(polyflux.*x); # get the total flux
  cvis_model = polyft * x / flux;
end

"""
    poly_to_cvis(x, star; intensity_model=:linear, band=nothing)

Flux-normalised complex visibilities of one component, from the `polyft` that `setup_oi!`
stored in `star`.

With `intensity_model = :linear` (default) `x` is the surface brightness (Rayleigh–Jeans
proxy, adequate in H and K). With `:planck`, `x` is a temperature map converted to intensity
by [`intensity`](@ref) at wavelength `band` (metres; see [`band_of`](@ref)). In R and V the
linear proxy misstates flux ratios by more than 10 %.
"""
function poly_to_cvis(x, star; intensity_model::Symbol = :linear, band = nothing)
   indx = star.index_quads_visible
   I = intensity_model === :linear ? x : intensity(x, intensity_model, band)
   xw = I[indx] .* star.vis_weights[indx] .* star.ldmap[indx]  # fold in limb darkening
   flux = dot(star.polyflux, xw); # get the total flux
   return star.polyft * xw / flux;
end

"""
    poly_to_flux(x, star; intensity_model=:linear, band=nothing)

Total (un-normalized) flux of a single component; see [`poly_to_cvis`](@ref) for the
intensity-model keywords.
"""
function poly_to_flux(x, star; intensity_model::Symbol = :linear, band = nothing)
   indx = star.index_quads_visible
   I = intensity_model === :linear ? x : intensity(x, intensity_model, band)
   xw = I[indx] .* star.vis_weights[indx] .* star.ldmap[indx]  # fold in limb darkening
   flux = dot(star.polyflux, xw); # get the total flux
   return flux;
end

# One epoch's flux/FT operators. Deliberately a separate function rather than the body of
# the threaded loop below: if `stars` is not concretely typed (e.g. a Vector{Any}), the
# loop body's locals become type-unstable, Julia boxes them into the closure, and the box
# is then SHARED by every thread — so one thread's `indx` can be read by another between
# the proj_west and proj_north lookups, producing a DimensionMismatch. A function call
# gives each invocation its own frame and cannot be boxed that way.
@views function _setup_oi_epoch!(star, dat, ::Type{T}; alt::Bool=false) where {T}
  indx = star.index_quads_visible
  pjx = star.proj_west[indx, :]
  pjy = star.proj_north[indx, :]
  star.polyflux = setup_polyflux_single(pjx, pjy)
  star.polyft = alt ? setup_polyft_single_alt(dat.uv, pjx, pjy; T=T) :
                      setup_polyft_single(dat.uv, pjx, pjy; T=T)
  return nothing
end

"""
    setup_oi!(data, stars)

Compute and store in each `stars[i]` the projected tessel areas (`polyflux`) and the polygon
Fourier transform matrix (`polyft`, Nuv × visible tessels) for `data[i]`. Threaded over
epochs, or within the single epoch when there is one.

The precision is the element type of `stars` (`stellar_geometry{T}`), or of
`stars[1].vertices_xyz` for an untyped container.
"""
function setup_oi!(data, stars::AbstractVector{stellar_geometry{T}}) where {T}
  _setup_oi_dispatch!(data, stars, T)
end

function setup_oi!(data, stars)
  _setup_oi_dispatch!(data, stars, eltype(stars[1].vertices_xyz))
end

function _setup_oi_dispatch!(data, stars, ::Type{T}) where {T}
  nepochs = size(data,1);
  if nepochs>1
    Threads.@threads for i=1:nepochs
      _setup_oi_epoch!(stars[i], data[i], T)
    end
  else # single epoch, thread over calculation
    _setup_oi_epoch!(stars[1], data[1], T; alt=true)
  end
end

# Same boxing hazard as _setup_oi_epoch! — keep the per-epoch work in its own function.
@views function _polygon_ft_epoch(star, dat, ::Type{T}) where {T}
  indx = star.index_quads_visible
  pjx = star.proj_west[indx, :]
  pjy = star.proj_north[indx, :]
  return setup_polyflux_single(pjx, pjy), setup_polyft_single(dat.uv, pjx, pjy; T=T)
end

function setup_polygon_ft(data, star_epoch_geom)
  nepochs = size(data,1);
  T = eltype(star_epoch_geom[1].vertices_xyz);
  polyflux = Array{Array{T,1}}(undef,nepochs);
  polyft = Array{Array{Complex{T},2}}(undef,nepochs);
  Threads.@threads for i=1:nepochs
    polyflux[i], polyft[i] = _polygon_ft_epoch(star_epoch_geom[i], data[i], T)
  end
  return polyflux, polyft
end

@views function setup_polyflux_single(pjx, pjy)
  # Polyflux is the projected area of each pixel (shoelace formula)
  T = eltype(pjx);
  polyflux = T(0.5)*(
    pjx[:,1].*pjy[:,2]
  - pjx[:,2].*pjy[:,1]
  + pjx[:,2].*pjy[:,3]
  - pjx[:,3].*pjy[:,2]
  + pjx[:,3].*pjy[:,4]
  - pjx[:,4].*pjy[:,3]
  + pjx[:,4].*pjy[:,1]
  - pjx[:,1].*pjy[:,4]);
  return polyflux;
end

function stcis(x1,x2,y1,y2,kx,ky)
  return sinc.(kx*(x2-x1) + ky*(y2-y1)).*cis.(-π*(kx*(x2+x1)+ ky*(y2+y1))).*(ky*(x2-x1)-kx*(y2-y1))
end

@views function setup_polyft_single_alt(uv, pjx, pjy; T = float(real(eltype(pjx))))
  kx =  uv[1,:] * T(-pi / (180*3600000))
  ky =  uv[2,:] * T( pi / (180*3600000))
  x = [Array(pjx[:,i])' for i in 1:4]
  y = [Array(pjy[:,i])' for i in 1:4]

  term1 = Threads.@spawn stcis(x[1], x[2], y[1], y[2], kx, ky)
  term2 = Threads.@spawn stcis(x[2], x[3], y[2], y[3], kx, ky)
  term3 = Threads.@spawn stcis(x[3], x[4], y[3], y[4], kx, ky)
  term4 = Threads.@spawn stcis(x[4], x[1], y[4], y[1], kx, ky)

  factor = -im * T(1/(2pi)) ./ (kx .* kx + ky .* ky)
  
  # Wait for the tasks to complete and get their results
  term1 = fetch(term1)
  term2 = fetch(term2)
  term3 = fetch(term3)
  term4 = fetch(term4)

  polyft = factor .* (term1 + term2 + term3 + term4)
  return polyft
end

@views function setup_polyft_single(uv, pjx, pjy; T = float(real(eltype(pjx))))
  kx =  uv[1,:] * T(-pi / (180*3600000));
  ky =  uv[2,:] * T( pi / (180*3600000));
  # note: check definition sinc(x) = sin(pi*x)/(pi*x)
  polyft = -im*T(1/(2pi))*(((sinc.( (kx*transpose(pjx[:,2]-pjx[:,1])) + (ky*transpose(pjy[:,2]-pjy[:,1])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,2]+pjx[:,1])) +  (ky*transpose(pjy[:,2]+pjy[:,1])) )).* ( (ky*transpose(pjx[:,2]-pjx[:,1]))  - (kx*transpose(pjy[:,2]-pjy[:,1])) ) + sinc.( (kx*transpose(pjx[:,3]-pjx[:,2])) + (ky*transpose(pjy[:,3]-pjy[:,2])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,3]+pjx[:,2])) +  (ky*transpose(pjy[:,3]+pjy[:,2])) )).* ( (ky*transpose(pjx[:,3]-pjx[:,2]))  - (kx*transpose(pjy[:,3]-pjy[:,2])) )+ sinc.( (kx*transpose(pjx[:,4]-pjx[:,3])) + (ky*transpose(pjy[:,4]-pjy[:,3])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,4]+pjx[:,3])) +  (ky*transpose(pjy[:,4]+pjy[:,3])) )).* ( (ky*transpose(pjx[:,4]-pjx[:,3]))  - (kx*transpose(pjy[:,4]-pjy[:,3])) ) + sinc.( (kx*transpose(pjx[:,1]-pjx[:,4])) + (ky*transpose(pjy[:,1]-pjy[:,4])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,1]+pjx[:,4])) +  (ky*transpose(pjy[:,1]+pjy[:,4])) )).* ( (ky*transpose(pjx[:,1]-pjx[:,4]))  - (kx*transpose(pjy[:,1]-pjy[:,4])) )))./((kx.*kx+ky.*ky)));
  return polyft;
end




@views function setup_polyft_single(kx, ky, pjx, pjy; T = float(real(eltype(pjx))))
  polyft = -im*T(1/(2pi))*(((sinc.( (kx*transpose(pjx[:,2]-pjx[:,1])) + (ky*transpose(pjy[:,2]-pjy[:,1])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,2]+pjx[:,1])) +  (ky*transpose(pjy[:,2]+pjy[:,1])) )).* ( (ky*transpose(pjx[:,2]-pjx[:,1]))  - (kx*transpose(pjy[:,2]-pjy[:,1])) ) + sinc.( (kx*transpose(pjx[:,3]-pjx[:,2])) + (ky*transpose(pjy[:,3]-pjy[:,2])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,3]+pjx[:,2])) +  (ky*transpose(pjy[:,3]+pjy[:,2])) )).* ( (ky*transpose(pjx[:,3]-pjx[:,2]))  - (kx*transpose(pjy[:,3]-pjy[:,2])) )+ sinc.( (kx*transpose(pjx[:,4]-pjx[:,3])) + (ky*transpose(pjy[:,4]-pjy[:,3])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,4]+pjx[:,3])) +  (ky*transpose(pjy[:,4]+pjy[:,3])) )).* ( (ky*transpose(pjx[:,4]-pjx[:,3]))  - (kx*transpose(pjy[:,4]-pjy[:,3])) ) + sinc.( (kx*transpose(pjx[:,1]-pjx[:,4])) + (ky*transpose(pjy[:,1]-pjy[:,4])) ).*cis.(-T(pi)*( (kx*transpose(pjx[:,1]+pjx[:,4])) +  (ky*transpose(pjy[:,1]+pjy[:,4])) )).* ( (ky*transpose(pjx[:,1]-pjx[:,4]))  - (kx*transpose(pjy[:,1]-pjy[:,4])) )))./((kx.*kx+ky.*ky)));
  return polyft;
end





function mod360(x)
  return mod.(mod.(x.+180,360).+360, 360) .- 180
end

function cvis_to_v2(cvis, indx)
  v2_model = abs2.(cvis[indx]);
end

# T follows the visibilities rather than defaulting to Float32: it sets the precision of the
# 180/π conversion, and Float32(180/π) = 57.2957802 is only 7 digits, which shifted χ²_t3phi
# by ~1e-8 relative against an otherwise-Float64 chain.
function cvis_to_t3(cvis, indx1, indx2, indx3; T = real(float(eltype(cvis))))
  t3 = cvis[indx1].*cvis[indx2].*cvis[indx3];
  t3amp = abs.(t3);
  t3phi = angle.(t3)*T(180/pi);
  return t3, t3amp, t3phi
end

# function cvis_to_t4(cvis, indx1, indx2, indx3, indx4)
#   t4 = cvis[indx1].*cvis[indx2]./(cvis[indx3]*conj(cvis[indx4]));
#   t4amp = abs.(t4);
#   t4phi = angle.(t4)*180.0/pi;
#   return t4, t4amp, t4phi
# end

"""
    cvis_to_obs(cvis, data) -> (v2, t3amp, t3phi)

Convert complex visibilities to interferometric observables (V2, T3amp, T3phi).
Common step for both single-star and binary forward models.
"""
function cvis_to_obs(cvis, data)
  v2_model = cvis_to_v2(cvis, data.indx_v2);
  _, t3amp_model, t3phi_model = cvis_to_t3(cvis, data.indx_t3_1, data.indx_t3_2, data.indx_t3_3);
  return v2_model, t3amp_model, t3phi_model
end

"V², T3amp, T3phi on; visamp, visphi, flux, differential phase off."
const OI_DEFAULT_WEIGHTS = [1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0]

# Pass the visibilities through at whatever complex float type they already are — do NOT
# force ComplexF64. OITOOLS accepts `AbstractVector{<:Complex}` and (in recent versions)
# derives its working eltype from the inputs, so widening here would throw away the point
# of carrying Float32 through the forward model.
"""
    cvis_chi2(cvis, data, weights = OI_DEFAULT_WEIGHTS) -> chi2

χ² of complex visibilities against an `OIdata` block, **with a hand-written adjoint**.

Thin wrapper over OITOOLS' `cvis_to_chi2_fg`, which returns the value *and* the complex
adjoint source `g_cvis` from a single pass over the data. The `rrule` below simply hands
that back, so reverse-mode AD through this call costs one hand-written sweep rather than a
tape of `abs2`/`angle`/`mod360` broadcasts.

Use this instead of `cvis_to_obs` + three `sum(abs2, …)` whenever the χ² is going to be
differentiated. It is the analytic-component counterpart of
[`interferometric_chi2`](@ref), whose rrule is fused to the polygon-FT image path and so
cannot be reused here.

# Adjoint convention

OITOOLS documents `g_params = real(Jᵀ · g_cvis)` with `J = ∂V/∂p` — a **transpose, not an
adjoint**, i.e. no conjugation. ChainRules/Zygote instead want the cotangent
`∂L/∂Re(cvis) + i ∂L/∂Im(cvis)`, which is the *complex conjugate* of that. Hence the
`conj(g)` in the pullback.

Getting this backwards is not a uniformly wrong answer, which is what makes it dangerous:
it leaves flux-like parameters nearly right (they move `V` along the real direction) while
making phase-like parameters — separation, position angle, T0 — wrong by ~100%. Verified
against finite differences to 1e-11 in both directions; `test_reflection.jl` pins it.
"""
@inline _ascomplexfloat(v::AbstractVector{<:Complex{<:AbstractFloat}}) = v
@inline _ascomplexfloat(v::AbstractVector) = complex.(float.(v))

cvis_chi2(cvis, data) = cvis_chi2(cvis, data, OI_DEFAULT_WEIGHTS)
# Value-only path uses `cvis_to_chi2_f`, NOT `first(cvis_to_chi2_fg(...))`: the latter also
# builds and fills the adjoint source, which cost ~40 % on the primal for nothing whenever
# no gradient was wanted. The rrule below still uses `_fg`, which returns both in one sweep.
cvis_chi2(cvis, data, weights) =
    cvis_to_chi2_f(_ascomplexfloat(cvis), data; weights = weights)

function ChainRulesCore.rrule(::typeof(cvis_chi2), cvis, data, weights)
    chi2, g = cvis_to_chi2_fg(_ascomplexfloat(cvis), data; weights = weights)
    function cvis_chi2_pullback(c̄raw)
        return (NoTangent(), unthunk(c̄raw) .* conj(g), NoTangent(), NoTangent())
    end
    return chi2, cvis_chi2_pullback
end

"""
    observables(x, star, data; intensity_model=:linear, band=nothing) -> (v2, t3amp, t3phi)

Model observables for a single component. See [`poly_to_cvis`](@ref) for the
intensity-model keywords.
"""
function observables(x, star, data; intensity_model::Symbol = :linear, band = nothing)
  # AN EMPTY `polyft` MEANS "no matrix was built", not "an error". `poly_to_cvis` reads
  # `star.polyft`, which only exists after `setup_oi!`; when it does not, the visibilities are
  # computed by the matrix-free route instead of failing.
  #
  # That makes the choice automatic and puts it where it belongs. A caller with a FIXED
  # geometry and a varying map — imaging, hundreds of VMLMB iterations — runs `setup_oi!` once
  # and gets the matrix. A caller whose geometry changes every evaluation, or that just wants
  # one χ² for a table, skips it and pays nothing for a matrix it would read once. MEASURED on
  # polaris at HEALPix 3: `setup_oi!` alone is 741 ms, and the whole fused evaluation is 4 ms.
  cvis_model = isempty(star.polyft) ?
      fused_cvis(x, star, data; intensity_model=intensity_model, band=band) :
      poly_to_cvis(x, star; intensity_model=intensity_model, band=band);
  return cvis_to_obs(cvis_model, data)
end

"""
    chi2_breakdown(x, star, data; intensity_model=:linear, band=nothing) -> NamedTuple
    chi2_breakdown(x, stars, data::AbstractVector; ...) -> NamedTuple

χ² of a map, per observable, as a NamedTuple:

    (; v2, t3amp, t3phi, total,            # absolute
       nv2, nt3amp, nt3phi, ndata,         # counts
       v2r, t3ampr, t3phir, totalr)        # reduced

`totalr` divides by the number of data points, not by the degrees of freedom. An observable
with no data gives `NaN`. The multi-epoch method sums over epochs and adds `epochs`, the
per-epoch breakdowns. [`chi2s`](@ref) returns the absolute values only.
"""
function chi2_breakdown(x, star, data; intensity_model::Symbol = :linear, band = nothing)
  chi2_v2, chi2_t3amp, chi2_t3phi = chi2s(x, star, data; verbose = false,
                                          intensity_model = intensity_model, band = band)
  return _chi2_breakdown(chi2_v2, chi2_t3amp, chi2_t3phi,
                         data.nv2, data.nt3amp, data.nt3phi)
end

function chi2_breakdown(x, stars, data::AbstractVector; intensity_model::Symbol = :linear,
                        band = nothing)
  per = [chi2_breakdown(x, stars[i], data[i]; intensity_model = intensity_model, band = band)
         for i in eachindex(data)]
  tot = _chi2_breakdown(sum(e.v2 for e in per), sum(e.t3amp for e in per),
                        sum(e.t3phi for e in per), sum(e.nv2 for e in per),
                        sum(e.nt3amp for e in per), sum(e.nt3phi for e in per))
  return merge(tot, (; epochs = per))
end

# 0/0 is NaN for floats already, but the guard is explicit: an integer count of 0 would make
# this an error rather than a NaN, and "no T3 in this file" must not look like a failure.
_reduced(chi2, n) = n > 0 ? chi2 / n : oftype(chi2 / oneunit(chi2), NaN)

function _chi2_breakdown(v2, t3amp, t3phi, nv2, nt3amp, nt3phi)
  total = v2 + t3amp + t3phi
  ndata = nv2 + nt3amp + nt3phi
  return (; v2, t3amp, t3phi, total, nv2, nt3amp, nt3phi, ndata,
          v2r = _reduced(v2, nv2), t3ampr = _reduced(t3amp, nt3amp),
          t3phir = _reduced(t3phi, nt3phi), totalr = _reduced(total, ndata))
end

function chi2s(x, star, data; verbose::Bool = true,
               intensity_model::Symbol = :linear, band = nothing)
  v2_model, t3amp_model, t3phi_model = observables(x, star, data;
                                                   intensity_model=intensity_model, band=band);
  chi2_v2 = sum(abs2, (v2_model - data.v2)./data.v2_err);
  chi2_t3amp = sum(abs2, (t3amp_model - data.t3amp)./data.t3amp_err);
  chi2_t3phi = sum(abs2, mod360(t3phi_model - data.t3phi)./data.t3phi_err);
  if verbose == true
    printstyled(@sprintf("V2: %.4f ", chi2_v2/data.nv2), color=:red)
    printstyled(@sprintf("T3A: %.4f ", chi2_t3amp/data.nt3amp), color=:blue)
    printstyled(@sprintf("T3P: %.4f\n", chi2_t3phi/data.nt3phi), color=:green)
  end
  return chi2_v2, chi2_t3amp, chi2_t3phi
end

function chi2s2(x, star, data; verbose::Bool = true)
  v2_model, t3amp_model, t3phi_model = observables(x, star, data);
  chi2_v2 = n2((v2_model - data.v2)./data.v2_err);
  chi2_t3amp = n2((t3amp_model - data.t3amp)./data.t3amp_err);
  chi2_t3phi = n2( mod360(t3phi_model - data.t3phi)./data.t3phi_err);
  if verbose == true
    printstyled(@sprintf("V2: %.4f ", chi2_v2/data.nv2), color=:red)
    printstyled(@sprintf("T3A: %.4f ", chi2_t3amp/data.nt3amp), color=:blue)
    printstyled(@sprintf("T3P: %.4f\n", chi2_t3phi/data.nt3phi), color=:green)
  end
  return chi2_v2, chi2_t3amp, chi2_t3phi
end

function spheroid_chi2_f(x, star, data; verbose::Bool = false,
                         intensity_model::Symbol = :linear, band = nothing,
                         weights = (1, 1, 1))
  v2_model, t3amp_model, t3phi_model = observables(x, star, data;
                                                   intensity_model=intensity_model, band=band);
  chi2_v2 = sum(abs2, (v2_model - data.v2)./data.v2_err);
  chi2_t3amp = sum(abs2, (t3amp_model - data.t3amp)./data.t3amp_err);
  chi2_t3phi = sum(abs2, mod360(t3phi_model - data.t3phi)./data.t3phi_err);
  if verbose == true
    printstyled(@sprintf("V2: %.4f ", chi2_v2/data.nv2), color=:red)
    printstyled(@sprintf("T3A: %.4f ", chi2_t3amp/data.nt3amp), color=:blue)
    printstyled(@sprintf("T3P: %.4f\n", chi2_t3phi/data.nt3phi), color=:green)
  end
  return weights[1]*chi2_v2 + weights[2]*chi2_t3amp + weights[3]*chi2_t3phi
end

# `weights` scales `[V², T3amp, T3φ]` in the criterion AND its gradient, as in
# `parametric_chi2`: a zero drops that observable from the reconstruction (e.g. V² and T3amp
# only while the closure phases are not trusted). The verbose trace still prints all three.
@views function spheroid_chi2_fg(x, g, star, data; verbose::Bool = true,
                                 intensity_model::Symbol = :linear, band = nothing,
                                 weights = (1, 1, 1))
  npix = star.npix;
  T = eltype(x);
  indx = star.index_quads_visible
  w = star.vis_weights[indx] .* star.ldmap[indx]  # soft visibility × limb darkening
  # THE SURFACE VARIABLE IS A TEMPERATURE once `intensity_model !== :linear`, and the forward
  # model needs its INTENSITY. `spheroid_chi2_f` (:363) and `poly_to_cvis` (:18) have taken
  # `intensity_model` since the Planck law was added; this gradient did not, so an imaging run
  # that asked for `:planck` got its VALUE from the Planck model and its GRADIENT from the
  # linear proxy — a descent on one function reported as another. Under `:linear` every line
  # below reduces to what it was.
  I = intensity_model === :linear ? x : intensity(x, intensity_model, band)
  xw = I[indx] .* w           # weighted pixel values
  cvis_model = poly_to_cvis(x, star; intensity_model = intensity_model, band = band);
  v2_model = cvis_to_v2(cvis_model, data.indx_v2);
  t3_model, t3amp_model, t3phi_model = cvis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2 ,data.indx_t3_3);
  chi2_v2 = sum(abs2, (v2_model - data.v2)./data.v2_err);
  chi2_t3amp = sum(abs2, (t3amp_model - data.t3amp)./data.t3amp_err);
  chi2_t3phi = sum(abs2, mod360(t3phi_model - data.t3phi)./data.t3phi_err);
  # Gradient w.r.t. weighted pixels (xw), computed via polyft^T @ adjoint_signal
  g_v2 = real(transpose(star.polyft[data.indx_v2,:])*(4*((v2_model-data.v2)./data.v2_err.^2).*conj(cvis_model[data.indx_v2])));
  g_t3amp = real(transpose(star.polyft[data.indx_t3_1,:]) *(2*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*conj(cvis_model[data.indx_t3_1])./abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_3]) ))+real(transpose(star.polyft[data.indx_t3_2,:])*(2*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*conj(cvis_model[data.indx_t3_2])./abs.(cvis_model[data.indx_t3_2]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_3]) ))+real(transpose(star.polyft[data.indx_t3_3,:])*(2*((t3amp_model-data.t3amp)./data.t3amp_err.^2).*conj(cvis_model[data.indx_t3_3])./abs.(cvis_model[data.indx_t3_3]).*abs.(cvis_model[data.indx_t3_1]).*abs.(cvis_model[data.indx_t3_2]) ));
  g_t3phi = T(360/pi)*imag(transpose(star.polyft[data.indx_t3_1,:])*(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*cvis_model[data.indx_t3_2].*cvis_model[data.indx_t3_3].*conj(t3_model))+transpose(star.polyft[data.indx_t3_2,:])*(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_3].*conj(t3_model))+transpose(star.polyft[data.indx_t3_3,:])*(((mod360(t3phi_model-data.t3phi)./data.t3phi_err.^2)./abs2.(t3_model)).*cvis_model[data.indx_t3_1].*cvis_model[data.indx_t3_2].*conj(t3_model)));
  gsum = T(weights[1])*g_v2 + T(weights[2])*g_t3amp + T(weights[3])*g_t3phi;
  flux = poly_to_flux(x, star; intensity_model = intensity_model, band = band);
  # Gradient w.r.t. x: chain rule through soft visibility weights
  # gsum is ∂χ²/∂(xw), so ∂χ²/∂I = w .* ∂χ²/∂(xw) (after flux normalization)
  g_normalized = (gsum .- dot(xw, gsum) * star.polyflux / flux) / flux;
  # AND THROUGH THE INTENSITY LAW. It is pointwise, so dI/dx is diagonal and `_intensity_dT`
  # returns it; under `:linear` it is identically 1 and this is the old expression.
  if intensity_model === :linear
    g[indx] = w .* g_normalized;
  else
    dI = _intensity_dT(x, intensity_model, band)
    g[indx] = w .* g_normalized .* dI[indx];
  end
  if verbose == true
    printstyled(@sprintf("V2: %.4f ", chi2_v2/data.nv2), color=:red)
    printstyled(@sprintf("T3A: %.4f ", chi2_t3amp/data.nt3amp), color=:blue)
    printstyled(@sprintf("T3P: %.4f ", chi2_t3phi/data.nt3phi), color=:green)
    printstyled(@sprintf("Flux: %.4f\n", flux), color=:normal)
  end
  return weights[1]*chi2_v2 + weights[2]*chi2_t3amp + weights[3]*chi2_t3phi
end

function spheroid_chi2_allepochs_f(x, stars, data; epochs_weights=[], verbose=false,
                                   weights = (1, 1, 1))
nepochs = length(data)
chi2_t = zeros(eltype(x), nepochs);
Threads.@threads for i=1:nepochs # weighted sum -- should probably do the computation in parallel
  chi2_t[i] = spheroid_chi2_f(x, stars[i], data[i], verbose=verbose, weights=weights);
end
# The weights go on the per-epoch terms, BEFORE the sum. `f = f.*w` after summing — which is
# what stood here — multiplies a scalar by a vector and returns a VECTOR from a function whose
# every caller expects a number.
w = _epoch_weights(epochs_weights, nepochs, eltype(x))
f = w === nothing ? sum(chi2_t) : sum(chi2_t .* w)
if verbose == true 
  printstyled(@sprintf("All epochs, total χ²: %.4f\n\n", f), color=:white);
end
return f;
end

"""
    _epoch_weights(w, nepochs, T) -> Vector{T} or nothing

Validate a per-epoch weight vector. `nothing` means "unweighted", which is the common case and
is worth not multiplying by a vector of ones.

Rejected rather than recycled or padded: a weight vector of the wrong length is a mistake about
which epochs are being weighted, and quietly broadcasting it would weight the wrong nights.
"""
function _epoch_weights(w, nepochs::Int, ::Type{T}) where T
    (w === nothing || isempty(w)) && return nothing
    length(w) == nepochs ||
        throw(DimensionMismatch("epochs_weights has $(length(w)) entries for $(nepochs) " *
                                "epochs"))
    any(x -> x < 0 || !isfinite(x), w) &&
        throw(ArgumentError("epochs_weights must be finite and non-negative"))
    return collect(T, w)
end

function spheroid_crit_allepochs_fg(x, g, stars, data; regularizers=[], epochs_weights=[],
                                    verbose=false, T=eltype(x),
                                    intensity_model::Symbol = :linear, band = nothing,
                                    weights = (1, 1, 1))
  #g[:] .= T(0);
  nepochs = length(data)
  chi2_t = zeros(T, nepochs);
  #npix = stars[1].npix
  singleepoch_g = [zeros(T, length(x)) for i=1:nepochs];
  Threads.@threads for i=1:nepochs # weighted sum -- should probably do the computation in parallel
    chi2_t[i] = spheroid_chi2_fg(x, singleepoch_g[i], stars[i], data[i], verbose=verbose,
                                 intensity_model=intensity_model, band=band, weights=weights);
  end
  # PER-EPOCH WEIGHTS, applied to the criterion AND to the gradient. This warned "not
  # implemented" and silently ignored the argument; the commented-out attempt beside it
  # (`f = f.*epochs_weights`) could not have worked either, because `f` is already the scalar
  # sum by that point and the weights have to go on the per-epoch terms.
  #
  # What they are for: a night with a bad calibration should be able to count for less without
  # being thrown away, and down-weighting is the only alternative to dropping it.
  w = _epoch_weights(epochs_weights, nepochs, T)
  f = w === nothing ? sum(chi2_t) : sum(chi2_t .* w)
  # The gradient carries the same weights, or the map is optimised against a different
  # criterion from the one being reported — which converges, to the wrong thing.
  g[:] .= w === nothing ? sum(singleepoch_g) :
          sum(w[i] .* singleepoch_g[i] for i in 1:nepochs)
  if verbose == true
    printstyled(@sprintf("Total χ²: %.4f\n", f), color=:white);
  end
  # Map regularization
  if regularizers!=[]
    reg_g = zeros(T, length(x));
    f += spheroid_regularization(x, reg_g, regularizers=regularizers, verbose = verbose);
    g[:] += reg_g
  end
  return f;
end

function parametric_temperature_map(parameters, star; secondary=false) # dispatches parametric
  if star.surface_type == 3
    return temperature_map_vonZeipel_roche_single(parameters, star, star.t, secondary=secondary);
  elseif star.surface_type == 2
    # Whichever gravity-darkening law `parameters.gravity_law` names; von Zeipel when the
    # field is absent, which is what every model written before the laws became selectable
    # means. See `gravity_law_spec`.
    return temperature_map_rapid_rotator(parameters,star);
  elseif star.surface_type == 1
    return temperature_map_vonZeipel_ellipsoid(parameters,star);
  elseif star.surface_type == 0 # sphere
    T = eltype(star.vertices_spherical)
    return T(parameters.tpole)*ones(T, star.npix);
  else
    println("Unimplemented parametric von Zeipel function")
  end
  return
end

function spheroid_parametric_f(parameters, tessels, data, tepochs) 
  stars = create_star_multiepochs(tessels, parameters, tepochs);
  setup_oi!(data, stars) 
  x = parametric_temperature_map(parameters, stars[1]);
  return spheroid_chi2_allepochs_f(x, stars, data)
end

function proj_positivity(ztilde)
z = copy(ztilde)
z[ztilde.>0]=0
return z
end


# function spheroid_total_variation_fg(x, tv_g, tvinfo; ϵ = 1e-13, T=Float32, verbose = true)
#   # Add total variation regularization
#   #tvinfo: neighbors,south_neighbors,west_neighbors,south_neighbors_reverse,west_neighbors_reverse
#   ϵ = T(ϵ)
#   npix = length(x)
#   xs = x[tvinfo[2]];
#   xw = x[tvinfo[3]];
#   tv_f = sum(sqrt.( (x-xs).^2 + (x-xw).^2) .+ ϵ )
#   tv_g[1:npix] = (2*x-xs-xw)./(sqrt.( (x-xs).^2 + (x-xw).^2) .+ ϵ)
#   for j=1:length(x)
#     k = tvinfo[4][j]
#     l = tvinfo[5][j]
#     if length(k)>0
#       kk=k[1]
#       tv_g[j] += (x[j]-x[kk])/sqrt((x[kk]-x[j])^2+(x[kk]-x[tvinfo[3][kk]])^2 +ϵ)
#     end
#     if length(l)>0
#       ll=l[1]
#       tv_g[j] += (x[j]-x[ll])/sqrt((x[ll]-x[tvinfo[2][ll]])^2+(x[ll]-x[j])^2 +ϵ)
#     end
#    end
# if verbose == true
#       println("TV: ", tv_f);
# end
#   return tv_f
# end


function spheroid_total_variation2_fg(x, tv_g, tvinfo; verbose = true)
  npix = length(x)
  tv_f = norm(tvinfo[6]*x)^2
  tv_g[:] = 2*(tvinfo[7]*x)
  if verbose == true
      printstyled(@sprintf("TV2: %.4f\n", tv_f), color=:yellow);
  end
  return tv_f
end

function spheroid_total_variation_fg(x, tv_g, tvinfo; verbose = true)
  npix = length(x)
  tv_f = norm(tvinfo[6]*x)
  if tv_f>0 
  tv_g[:] = (tvinfo[7]*x)/tv_f
  else
    tv_g[:] .= 0
  end
  if verbose == true
      println("TV: ", tv_f);
  end
  return tv_f
end

function spheroid_mean_fg(x, g; verbose = true)
f = sum(abs.(x.-sum(x)/length(x)))
g[:] = sign.(x.-sum(x)/length(x))
if verbose == true
println(" MeanReg: ", f);
end
return f;
end

function spheroid_harmon_bias_fg(x, g, B::Float64; verbose = true)
n = length(x);
avg_x = mean(x);
bcorr = (B-1.0)*Int.((x.-avg_x).>0).+1.0
#reg_f = sum(bcorr.*(x.-avg_x).^2);
reg_f = sum(bcorr.*(x.-avg_x).^2)/n;
#reg_g = 2*(n-1)/n*bcorr.*(x.-avg_x)
reg_g = 2/n*bcorr.*(x.-avg_x)
g[:] = reg_g .- mean(reg_g); # necessary ?

if verbose == true
println(" Bias: ", reg_f);
end
return reg_f;
end

"""
    max_entropy_fg(x, g; verbose=false, ϵ=1e-9) -> f

Maximum entropy of the map normalised by its mean (`"mem"`):

```
f = Σᵢ xmᵢ log(xmᵢ + ϵ),      xm = x / mean(x)
```

Invariant under `x → c·x`, as χ² is. Pointwise: it imposes no spatial correlation.

# Gradient

With `m = mean(x)`, `n = length(x)` and `dᵢ = log(xmᵢ + ϵ) + xmᵢ/(xmᵢ + ϵ)`,

```
∂f/∂xⱼ = dⱼ/m − (1/(n m²)) Σᵢ dᵢ xᵢ
```
"""
function max_entropy_fg(x, g; verbose=false, ϵ=1e-9)
  n = length(x)
  m = sum(x)/n
  abs(m) > ϵ || (g .= 0; return zero(m))
  xm = x ./ m
  d  = log.(xm .+ ϵ) .+ xm ./ (xm .+ ϵ)
  reg_f = sum(xm .* log.(xm .+ ϵ))
  glob  = sum(d .* x) / (n*m^2)
  g[:] = d ./ m .- glob
  verbose && println(" MEM: ", reg_f)
  return reg_f
end

"""
    radflat_bins(star; nbins=6, subset=nothing) -> NamedTuple

Bin the visible tessels (or `subset`) by projected radius, for `"radflat"` and
`"radialvar"`. Pass the result as the third element of the regulariser entry and its `idx`
as the fourth.

`ρ` is the distance of each tessel's centre (mean of its four projected corners) from the
disk centre, normalised to 1 at the outermost tessel; the `nbins` bins are equal in `ρ`. Each
bin's weight is `w_b = ρ_b² · Σ polyflux` over the bin, with `ρ_b` the bin centre.

Returns `(idx, bin, nbins, w, count, rho)`. Warns when a bin is empty or holds fewer than 3
tessels; the innermost bin empties first, and nside ≥ 8 (`tessellation_healpix(3)`) is
stable to about 1 %.
"""
function radflat_bins(star; nbins::Int = 6, subset = nothing)
    idx = subset === nothing ? star.index_quads_visible : subset
    pw  = @view star.proj_west[idx, :]
    pn  = @view star.proj_north[idx, :]
    n   = length(idx)
    ρ   = [hypot(sum(@view pw[i, :])/4, sum(@view pn[i, :])/4) for i in 1:n]
    ρmax = maximum(ρ)
    ρmax > 0 || error("radflat_bins: all projected radii are zero")
    ρn = ρ ./ ρmax
    # bin edges equispaced in ρ; the outermost patch lands in the last bin, not one past it
    b  = [min(nbins, floor(Int, r*nbins) + 1) for r in ρn]
    pf = length(star.polyflux) == n ? abs.(star.polyflux) : abs.(star.polyflux[idx])
    ρc = [(k - 0.5)/nbins for k in 1:nbins]                 # bin centre radius
    w  = zeros(Float64, nbins)
    cnt = zeros(Int, nbins)
    @inbounds for i in 1:n
        w[b[i]]  += pf[i]
        cnt[b[i]] += 1
    end
    w .= (ρc .^ 2) .* w                                     # ρ_b² · flux_weight_b
    # Under-populated bins are the failure mode at coarse tessellations. Equal-ρ bins on a
    # projected disk hold ∝ ρ dρ patches, so the INNERMOST bin starves first: at HEALPix
    # nside=1 (24 visible patches) one bin comes out empty and the penalty is ~30 % wrong;
    # at nside=2 the innermost bin holds a single patch, so its "mean brightness" is one
    # number. An empty bin contributes zero rather than erroring — quietly, which is exactly
    # the kind of silent wrong answer worth a warning. nside ≥ 3 is stable to ~1 %.
    if any(iszero, cnt)
        @warn "radflat_bins: $(count(iszero, cnt)) of $nbins radial bins are EMPTY " *
              "($(cnt)); RADFLAT will be biased. Use a finer tessellation or fewer bins."
    elseif minimum(cnt) < 3
        @warn "radflat_bins: thinnest radial bin holds only $(minimum(cnt)) patch(es) " *
              "($(cnt)); RADFLAT will be noisy. Consider a finer tessellation or fewer bins."
    end
    return (idx = idx, bin = b, nbins = nbins, w = w, count = cnt, rho = ρn)
end

"""
    spheroid_radflat_fg(x, g, bins; scale=1e5, verbose=false) -> f

RADFLAT (J. Monnier, priv. comm. 2026-07-01): penalise departures of the azimuthally averaged
radial brightness profile from flat.

```
f = scale · Σ_b w_b · (I_b/I_mean − 1)²
```

`I_b` is the mean of `x` in radial bin `b`, `I_mean` the mean over all binned tessels, and
`bins` comes from [`radflat_bins`](@ref). `x` is the map without limb darkening, so the
radial profile is left to the limb-darkening law. An empty bin contributes zero. Raises if
`x` and `bins` have different lengths.

Intended for a non-rotating star observed at one epoch.

# Gradient

With `n_b` tessels in bin `b`, `N` in total and `r_b = I_b/I_mean`,

```
∂f/∂x_j = 2·scale·[ w_{b(j)}(r_{b(j)} − 1)/(n_{b(j)}·I_mean)
                    − (1/(N·I_mean))·Σ_b w_b (r_b − 1) r_b ]
```
"""
function spheroid_radflat_fg(x, g, bins; scale = 1e5, verbose = false, ϵ = 1e-12)
    nb = bins.nbins; b = bins.bin; w = bins.w; cnt = bins.count
    N  = length(x)
    N == length(b) ||
        error("spheroid_radflat_fg: x has $N entries but the binning was built for $(length(b)). " *
              "The regularizer's pixel subset must match the one passed to radflat_bins.")
    Imean = sum(x)/N
    abs(Imean) > ϵ || (g .= 0; return 0.0)
    S = zeros(Float64, nb)
    @inbounds for i in 1:N; S[b[i]] += x[i]; end
    r = [cnt[k] > 0 ? (S[k]/cnt[k])/Imean : 1.0 for k in 1:nb]   # empty bin contributes nothing
    f = scale * sum(w[k]*(r[k]-1)^2 for k in 1:nb)
    glob = sum(w[k]*(r[k]-1)*r[k] for k in 1:nb) / (N*Imean)
    @inbounds for i in 1:N
        k = b[i]
        loc = cnt[k] > 0 ? w[k]*(r[k]-1)/(cnt[k]*Imean) : 0.0
        g[i] = 2*scale*(loc - glob)
    end
    verbose && println(" RADFLAT: ", f)
    return f
end

"""
    spheroid_radialvar_fg(x, g, bins; scale=1e5, normalize=true, verbose=false) -> f

RADIALVAR: penalise the azimuthal scatter of the map within each projected radial bin.

```
f = scale · Σ_b var_b / I_mean²,        var_b = (1/(n_b−1)) Σ_{j∈b} (x_j − x̄_b)²
```

Uses the binning of [`radflat_bins`](@ref) and, like [`spheroid_radflat_fg`](@ref), the map
without limb darkening. A bin with one tessel contributes zero. RADFLAT and RADIALVAR are the
between-bin and within-bin terms of the variance of the map (one-way ANOVA decomposition).

This is the tessellated counterpart of OITOOLS' `"radialvar"` (`radial_variance`), which
works on an image grid with elliptical annuli and sums raw variances. `normalize = true`
divides by `I_mean²`; `normalize = false` gives OITOOLS' scaling.

# Gradient

With `c_b = 1/(n_b−1)` and `T = Σ_b c_b Σ_{j∈b}(x_j − x̄_b)²`,

```
∂f/∂x_j = (2·scale / I_mean²) · [ c_{b(j)}(x_j − x̄_{b(j)}) − T/(N·I_mean) ]
```

The second term vanishes when `normalize = false` (and `I_mean²` is replaced by 1).
"""
function spheroid_radialvar_fg(x, g, bins; scale = 1e5, normalize::Bool = true,
                               verbose = false, ϵ = 1e-12)
    nb = bins.nbins; b = bins.bin; cnt = bins.count
    N  = length(x)
    N == length(b) ||
        error("spheroid_radialvar_fg: x has $N entries but the binning was built for " *
              "$(length(b)). The regularizer's pixel subset must match the one passed to " *
              "radflat_bins.")
    Imean = sum(x)/N
    if normalize && abs(Imean) <= ϵ
        g .= 0; return 0.0
    end
    S = zeros(Float64, nb)
    @inbounds for i in 1:N; S[b[i]] += x[i]; end
    mean_b = [cnt[k] > 0 ? S[k]/cnt[k] : 0.0 for k in 1:nb]
    c = [cnt[k] > 1 ? 1/(cnt[k] - 1) : 0.0 for k in 1:nb]    # a 1-patch bin has no variance
    T = 0.0
    @inbounds for i in 1:N
        d = x[i] - mean_b[b[i]]
        T += c[b[i]] * d * d
    end
    den  = normalize ? Imean^2 : 1.0
    f    = scale * T / den
    glob = normalize ? T / (N * Imean) : 0.0
    @inbounds for i in 1:N
        g[i] = 2 * scale / den * (c[b[i]] * (x[i] - mean_b[b[i]]) - glob)
    end
    verbose && println(" RADIALVAR: ", f)
    return f
end

"""
    orthold_direction(star, x0, star_params; subset=nothing) -> NamedTuple

The map direction that is degenerate with the limb-darkening coefficient `ld1`, for an
`"orthold"` regulariser entry (third element; `.idx` is the fourth).

Each tessel enters the forward model weighted by `vis_weights · ldmap`. A change `δ` in `ld1`
changes the weighted map by `x₀ ∘ vis_weights ∘ ∂ldmap/∂ld1 · δ`, so a map change along

```
v = x₀ ∘ vis_weights ∘ (∂ldmap/∂ld1),   v̂ = v/‖v‖
```

is indistinguishable from it in the data. For the Hestroffer law (`ldmap = μ^α`) this is
`δx ∝ x₀ ln μ`, which grows toward the limb; its weight in the data,
`∝ μ^{2·ld1} |ln μ|`, peaks at `μ = exp(−1/(2·ld1))` (0.73 for ld1 = 1.6).

Returns `(idx, c, x0, m0)` with `c = vis_weights ∘ ldmap ∘ v̂` and `m0 = mean(x₀)`. `x0` is
the reference map, normally the parametric starting map. The direction uses the surface
normals, so it applies to any surface type. Raises if `ld1` has no effect on the map or if
`x0` has zero mean.
"""
function orthold_direction(star, x0, star_params; subset = nothing)
    idx = subset === nothing ? star.index_quads_visible : subset
    ldtype = Int(star_params.ldtype)
    ld1 = Float64(star_params.ld1)
    ld2 = Float64(hasproperty(star_params, :ld2) ? star_params.ld2 : 0.0)
    nz  = Float64.(star.normals[idx, 3])
    ld3 = Float64(hasproperty(star_params, :ld3) ? star_params.ld3 : 0.0)
    ld4 = Float64(hasproperty(star_params, :ld4) ? star_params.ld4 : 0.0)
    ld, _, dld_dld1, _ = ld_and_derivs(nz, ldtype, ld1, ld2, ld3, ld4)
    vis = Float64.(star.vis_weights[idx])
    x0v = Float64.(x0[idx])
    v   = x0v .* vis .* dld_dld1                 # model-space degenerate direction
    nv  = sqrt(sum(abs2, v))
    nv > 0 || error("orthold_direction: the degenerate direction is identically zero " *
                    "(ld1 has no effect on this map — is ldtype/ld1 set?)")
    v ./= nv
    m0 = sum(x0v) / length(x0v)
    abs(m0) > 0 || error("orthold_direction: reference map has zero mean")
    return (idx = idx, c = vis .* ld .* v, x0 = x0v, m0 = m0)
end

"""
    spheroid_orthold_fg(x, g, od; scale=1e5, verbose=false) -> f

orthoLD: penalise the component of `x − x₀` along the limb-darkening-degenerate direction
from [`orthold_direction`](@ref).

```
f = scale · ⟨x − x₀, c⟩² / m₀²,     ∂f/∂x = 2·scale·⟨x − x₀, c⟩/m₀² · c
```

Rank one: it removes one degree of freedom. `m₀` is the mean of the reference map.
"""
function spheroid_orthold_fg(x, g, od; scale = 1e5, verbose = false)
    length(x) == length(od.c) ||
        error("spheroid_orthold_fg: x has $(length(x)) entries but the direction was built " *
              "for $(length(od.c)). The regularizer's pixel subset must match `od.idx`.")
    P = 0.0
    @inbounds for i in eachindex(x); P += (x[i] - od.x0[i]) * od.c[i]; end
    P /= od.m0
    f = scale * P * P
    k = 2 * scale * P / od.m0
    @inbounds for i in eachindex(x); g[i] = k * od.c[i]; end
    verbose && println(" orthoLD: ", f)
    return f
end

"""
    spheroid_sobel2_fg(x, g, S; scale=1.0, verbose=false, ϵ=1e-12) -> f

Quadratic gradient penalty (`"sobel2"`) on the operator `S = sobel_gradient_healpix(n)`:

```
f = scale · (4π/npix) · ‖D x‖² / mean(x)²   ≈   scale · ∫|∇x|²dΩ / x̄²
```

Invariant under `x → c·x`, as χ² is, and independent of `nside` for a smooth map, so a
weight carries over between resolutions. Typical weights are O(1–10²).

Compared with `"tv2"` (`‖Lx‖²`, `L` the graph Laplacian):

| | penalises | Fourier response | scale-invariant | nside-stable (smooth map) |
|---|---|---|---|---|
| `"tv2"` | `∇²x` | `k⁴` | no | no (×0.36 per doubling) |
| `"sobel2"` | `∇x` | `k²` | yes | yes |

# Gradient

With `Q = ‖D x‖²`, `m = mean(x)`, `n = length(x)` and `A = 4π/npix`,

```
∂f/∂xⱼ = scale·A·[ 2(Dᵀ D x)ⱼ / m² − 2Q/(n m³) ]
```
"""
function spheroid_sobel2_fg(x, g, S; scale = 1.0, verbose = false, ϵ = 1e-12)
    n = length(x)
    m = sum(x)/n
    abs(m) > ϵ || (g .= 0; return zero(m))
    # Squared neighbour differences (`S.D`), NOT |Gx x|² + |Gy x|²: the least-squares gradient
    # is blind to checkerboards, the squared differences are not. See `sobel_gradient_healpix`.
    d  = S.D*x
    Q  = sum(abs2, d)
    A  = S.area
    f  = scale * A * Q / m^2
    g[:] = (2*scale*A/m^2) .* (S.D'*d) .- (2*scale*A*Q/(n*m^3))
    verbose && println(" SOBEL2: ", f)
    return f
end

"""
    spheroid_sobel_fg(x, g, S; scale=1.0, verbose=false, ϵ=1e-8) -> f

Isotropic, edge-preserving total variation (`"sobel"`) on the operator
`S = sobel_gradient_healpix(n)`:

```
f = scale · (4π/npix) · Σᵢ √(|∇x|ᵢ² + ϵ) / mean(x)   ≈   scale · ∫|∇x|dΩ / x̄
```

The square root is taken per tessel, so an edge costs in proportion to its height. ROTIR's
`"tv"` is a different quantity, `‖Lx‖` (a global norm of the graph Laplacian, `k⁴` response).
Use `"sobel"` to keep edges such as spot boundaries, `"sobel2"` for a smooth map.

`ϵ` (in units of `|∇x|²` of the mean-normalised map) keeps the square root differentiable
at `∇x = 0`.

# Gradient

With `rᵢ = √(|∇x|ᵢ² + ϵ)`, `m = mean(x)`, `n = length(x)` and `A = 4π/npix`,

```
∂f/∂xⱼ = scale·A·[ (Dᵀ(D x ./ r[owner]))ⱼ / m − (Σᵢrᵢ)/(n m²) ]
```
"""
function spheroid_sobel_fg(x, g, S; scale = 1.0, verbose = false, ϵ = 1e-8)
    n = length(x)
    m = sum(x)/n
    abs(m) > ϵ || (g .= 0; return zero(m))
    # |∇x|²ᵢ from the squared neighbour differences each tessel owns (`S.D`, `S.owner`), then
    # the √ per tessel: isotropic graph total variation, with no checkerboard null space.
    d  = S.D*x
    q  = zeros(eltype(d), n)
    @inbounds for k in eachindex(d); q[S.owner[k]] += d[k]^2; end
    r  = sqrt.(q .+ ϵ)
    R  = sum(r)
    A  = S.area
    f  = scale * A * R / m
    g[:] = (scale*A/m) .* (S.D'*(d ./ r[S.owner])) .- (scale*A*R/(n*m^2))
    verbose && println(" SOBEL: ", f)
    return f
end

"""
    spheroid_laplace_beltrami_fg(x, g, S; scale=1.0, verbose=false, ϵ=1e-12) -> f

Dirichlet energy (`"laplace_beltrami"`) of the cotangent Laplace–Beltrami operator
`S = laplace_beltrami_healpix(n)`, normalised like `"sobel2"`:

```
f = scale · Σₑ wₑ (x_j − x_i)² / mean(x)²   ≈   scale · ∫|∇x|²dΩ / x̄²
```

Weights carry over from `"sobel2"` and between resolutions. Invariant under `x → c·x`; the
only null mode is the constant map.

# Gradient

With `Q = ‖D x‖²`, `m = mean(x)`, `n = length(x)`:
`∂f/∂xⱼ = scale·[ 2(Dᵀ D x)ⱼ / m² − 2Q/(n m³) ]`.
"""
function spheroid_laplace_beltrami_fg(x, g, S; scale = 1.0, verbose = false, ϵ = 1e-12)
    n = length(x)
    m = sum(x)/n
    abs(m) > ϵ || (g .= 0; return zero(m))
    d = S.D*x
    Q = sum(abs2, d)
    f = scale * Q / m^2
    g[:] = (2*scale/m^2) .* (S.D'*d) .- (2*scale*Q/(n*m^3))
    verbose && println(" LAPLACE-BELTRAMI: ", f)
    return f
end

function spheroid_regularization(x,g; printcolor = :black, regularizers=[], verbose=false)
  reg_f = 0.0;
  for ireg =1:length(regularizers)
      x_sub = x[regularizers[ireg][4]] # take the pixel subset if needed (example: only regularize visible pixels)
      temp_g = similar(x_sub);
      if regularizers[ireg][1] == "mem"
          reg_f += regularizers[ireg][2]*max_entropy_fg(x_sub, temp_g, verbose = verbose);
      elseif regularizers[ireg][1] == "tv2"
          reg_f += regularizers[ireg][2]*spheroid_total_variation2_fg(x_sub, temp_g, regularizers[ireg][3], verbose = verbose);
      elseif regularizers[ireg][1] == "tv"
          reg_f += regularizers[ireg][2]*spheroid_total_variation_fg(x_sub, temp_g, regularizers[ireg][3], verbose = verbose);
      elseif regularizers[ireg][1] == "mean"
          reg_f += regularizers[ireg][2]*spheroid_mean_fg(x_sub, temp_g, verbose = verbose);
      elseif regularizers[ireg][1] == "bias"
          reg_f += regularizers[ireg][2]*spheroid_harmon_bias_fg(x_sub, temp_g, regularizers[ireg][3], verbose = verbose );
      elseif regularizers[ireg][1] == "radflat"
          reg_f += regularizers[ireg][2]*spheroid_radflat_fg(x_sub, temp_g, regularizers[ireg][3], verbose = verbose);
      elseif regularizers[ireg][1] == "radialvar"
          reg_f += regularizers[ireg][2]*spheroid_radialvar_fg(x_sub, temp_g, regularizers[ireg][3], verbose = verbose);
      elseif regularizers[ireg][1] == "orthold"
          reg_f += regularizers[ireg][2]*spheroid_orthold_fg(x_sub, temp_g, regularizers[ireg][3], verbose = verbose);
      elseif regularizers[ireg][1] == "sobel2"
          reg_f += regularizers[ireg][2]*spheroid_sobel2_fg(x_sub, temp_g, regularizers[ireg][3], verbose = verbose);
      elseif regularizers[ireg][1] == "sobel"
          reg_f += regularizers[ireg][2]*spheroid_sobel_fg(x_sub, temp_g, regularizers[ireg][3], verbose = verbose);
      elseif regularizers[ireg][1] == "laplace_beltrami"
          reg_f += regularizers[ireg][2]*spheroid_laplace_beltrami_fg(x_sub, temp_g, regularizers[ireg][3], verbose = verbose);
      else
          # An unrecognised name must RAISE. Falling through would contribute exactly
          # zero, so a typo would look like a regularizer that simply does nothing.
          error("spheroid_regularization: unknown regularizer \"$(regularizers[ireg][1])\"; " *
                "known: mem, tv, tv2, mean, bias, radflat, radialvar, orthold, sobel, sobel2, " *
                "laplace_beltrami")
      end
      g[regularizers[ireg][4]] += regularizers[ireg][2]*temp_g
  end
  if verbose == true
      println("\n");
  end
return reg_f
end


using OptimPackNextGen

"""
    image_reconstruct_oi(x_start, data, stars; maxiter=100, regularizers=[], lower=0,
                         upper=Inf, epochs_weights=[], weights=(1, 1, 1),
                         intensity_model=:linear, band=nothing, callback=nothing,
                         callback_every=25, verbose=true) -> x

Reconstruct a surface map by VMLMB (OptimPackNextGen), minimising χ² over all epochs plus
`regularizers` within the bounds `[lower, upper]`. Uses the matrix path: run `setup_oi!`
first.

`weights = (w_V², w_T3amp, w_T3φ)` scales each observable's χ²; `(1, 1, 0)` reconstructs from
V² and T3amp. `epochs_weights` scales each epoch's χ².

`callback(x, nevals, f)` is called on the first evaluation and every `callback_every`
evaluations. `x` is the optimiser's working vector: copy it to keep it. `nevals` counts
criterion evaluations, which exceed iterations when a line search re-evaluates. A callback
that throws is dropped with a warning and the reconstruction continues.
"""
function image_reconstruct_oi(x_start, data, stars; epochs_weights =[], printcolor= [], verbose = true, lower=0, upper=Inf, maxiter = 100, regularizers =[],
                              intensity_model::Symbol = :linear, band = nothing,
                              weights = (1, 1, 1),
                              callback = nothing, callback_every::Int = 25)
  x_sol = [];
  nevals = Ref(0)
  crit_imaging = function (x, g)
      f = spheroid_crit_allepochs_fg(x, g, stars, data, regularizers=regularizers,
                                     epochs_weights=epochs_weights, verbose = verbose,
                                     intensity_model = intensity_model, band = band,
                                     weights = weights)
      nevals[] += 1
      if callback !== nothing && (nevals[] == 1 || nevals[] % max(1, callback_every) == 0)
          # A callback that throws must not take the reconstruction with it: it is a progress
          # report, and the run is the thing worth finishing.
          try
              callback(x, nevals[], f)
          catch err
              @warn "image_reconstruct_oi: callback failed, continuing without it" exception = err
              callback = nothing
          end
      end
      return f
  end
  x_sol = OptimPackNextGen.vmlmb(crit_imaging, x_start, verb=verbose, lower=lower, upper=upper, maxiter=maxiter, blmvm=false, gtol=(0,1e-8));
  dummy = similar(x_sol);
  crit_opt = crit_imaging(x_sol,dummy);
  return x_sol
end

# `verbose` defaults to FALSE, not to `verbose`. These three read `verbose = verbose` until
# 2026-08, which is not a self-reference to the keyword but a lookup of a GLOBAL named
# `verbose` — there is none, so every call that did not pass `verbose=` explicitly raised
# `UndefVarError(:verbose)`. It went unnoticed because every caller in the repo does pass it.
# False is the right default here: unlike `image_reconstruct_oi`, these evaluate a criterion
# at a given map and are typically called in a scan loop (see demos/rwcep_radflat.jl,
# demos/alphacen_ld.jl), where per-call printing is noise.
function image_reconstruct_oi_crit(x, data, stars; regularizers = [], verbose = false)
  g = similar(x);
  crit = spheroid_crit_allepochs_fg(x, g, stars, data, regularizers=regularizers,epochs_weights=[],  verbose = verbose);
  return crit
end

function image_reconstruct_oi_chi2(x, data, stars; verbose = false)
  g = similar(x);
  crit = spheroid_crit_allepochs_fg(x, g, stars, data, regularizers=[], epochs_weights= [], verbose = verbose);
  return crit
end

function image_reconstruct_oi_chi2_fg(x, data, stars; verbose = false)
  g = similar(x);
  crit = spheroid_crit_allepochs_fg(x, g, stars, data, regularizers=[], epochs_weights= [], verbose = verbose);
  return crit,g
end

"""
    multires_reconstruct_oi(data, star_params, tepochs; n_start=2, n_end=4, maxiter=500,
                            reg_weight=10.0, reg_type="sobel2", verbose=true, kwargs...)
                            -> (tmap, stars)

Reconstruct at HEALPix levels `n_start` to `n_end`, starting from the parametric temperature
map and initialising each level by upsampling the previous result. Other keywords go to
[`image_reconstruct_oi`](@ref).

The same `reg_weight` is used at every level. `"sobel"` and `"sobel2"` are normalised by
solid angle and by `mean(x)`, so the weight means the same at each level; their weights are
O(1–10²). A `"tv2"` weight depends on the map's level and is O(10⁻⁵–10⁻⁴); its effect on a
smooth map falls by about ×0.4 per level.
"""
function multires_reconstruct_oi(data, star_params, tepochs; n_start=2, n_end=4, maxiter=500, reg_weight=10.0, reg_type="sobel2", verbose=true, kwargs...)
  tmap = nothing
  stars = nothing
  for n in n_start:n_end
    tessels = tessellation_healpix(n)
    stars = create_star_multiepochs(tessels, star_params, tepochs)
    if tmap === nothing
      tmap = parametric_temperature_map(star_params, stars[1])
    else
      tmap = vec(repeat(tmap, 1, 4)')  # upsample from previous level
    end
    setup_oi!(data, stars)
    reginfo = reg_type in ("sobel", "sobel2") ? sobel_gradient_healpix(n) :
                                                tv_neighbors_healpix(n)
    regularizers = [[reg_type, reg_weight, reginfo, 1:length(tmap)]]
    if verbose
      println("Multi-resolution: HEALPix level n=$n, npix=$(nside2npix(2^n))")
    end
    tmap = image_reconstruct_oi(tmap, data, stars; maxiter=maxiter, regularizers=regularizers, verbose=verbose, kwargs...)
  end
  return tmap, stars
end

"""
    rescale_temperature_tpole(tmap, stars, star_params)

Rescale a temperature map so that the temperature at the pole equals `star_params.tpole`.
`stars` can be a single `stellar_geometry` or a vector of them (multi-epoch).
Warns if the pole pixel is never visible across all epochs.
"""
function rescale_temperature_tpole(tmap, stars, star_params)
    star = stars isa AbstractVector ? stars[1] : stars
    colat = @view star.vertices_spherical[:, 5, 2]
    ipole = argmin(colat)
    # Check visibility
    if stars isa AbstractVector
        visible = any(ipole in s.index_quads_visible for s in stars)
    else
        visible = ipole in star.index_quads_visible
    end
    if !visible
        @warn "Pole pixel (index $ipole) is never visible — rescaling is based on an unconstrained pixel"
    end
    scale = star_params.tpole / tmap[ipole]
    return tmap .* scale
end

"""
    rescale_temperature_teff(tmap, stars, teff)

Rescale a temperature map so that the mean temperature of the visible pixels equals `teff`.
`stars` can be a single `stellar_geometry` or a vector of them (multi-epoch).
Visibility is the union across all epochs.
"""
function rescale_temperature_teff(tmap, stars, teff)
    vis = sometimes_visible(stars)
    if isempty(vis)
        @warn "No visible pixels — returning unmodified map"
        return copy(tmap)
    end
    mean_T = mean(tmap[vis])
    scale = teff / mean_T
    return tmap .* scale
end


# "Multitemporal" = dynamical reconstruction

# function oi_reconstruct_mutitemporal(x_start::Array{Float64,1}, data::Array{OIdata,1}, polyflux::Array{Array{Float64,1},1}, polyft::Array{Array{Complex{Float64},2},1}; epochs_weights =[], printcolor= [], verbose = true, maxiter = 100, regularizers =[])
#   crit_imaging = (x,g)->oi_multitemporal_fg(x, g, polyflux, polyft, data, regularizers=regularizers, epochs_weights=  epochs_weights);
#   x_sol = OptimPackNextGen.vmlmb(crit_imaging, x_start, verbose=verbose, lower=0, maxiter=maxiter, blmvm=false, gtol=(0,1e-8));
#   return reshape(x_sol, size(polyflux[1],1), length(data))
# end

# function oi_multitemporal_fg(x, g, polyflux, polyft, data; regularizers=[], epochs_weights = [], verbose=true)
#   # Explanation of the following: optimpack optimizes a vector, but we want an array of images
#   npix = size(polyflux[1],1);
#   chi2_f = 0.0;
#   g[:] .= 0.0;
#   #temp_g = deepcopy(g[:]);
#   #npix = size(x);
#   nepochs = length(data);
#   if epochs_weights == []
#     epochs_weights = ones(Float64, nepochs)/nepochs;
#   end

#   #singleepoch_g = zeros(Float64, npix);
#   #for i=1:nepochs # weighted sum -- should probably do the computation in parallel
#   #  chi2_f += epochs_weights[i]*spheroid_chi2_fg_alt(x[1:npix], singleepoch_g, polyflux[i], polyft[i], data[i], verbose=verbose);
#   #  g[1:npix] += epochs_weights[i]*singleepoch_g;
#   #end


#   for i=1:nepochs # weighted sum -- in the future, do the computation in parallel
#     tslice = 1+(i-1)*npix:i*npix; # temporal slice
#     subg = zeros(Float64, npix);
#     chi2_f += epochs_weights[i]*spheroid_chi2_fg_alt(x[tslice], subg, polyflux[i], polyft[i], data[i], verbose=verbose);
#     g[tslice] = epochs_weights[i]*subg;
#     #x[tslice] = x[1:npix];
#     #g[tslice] = g[1:npix];
#   end

#   # cross temporal regularization -- weight needs to be defined in the "regularizers" variable
#   if length(regularizers)>nepochs
#     if (regularizers[nepochs+1][1][1] == "temporal_tvsq")  & (regularizers[nepochs+1][1][2] > 0.0) & (nepochs>1)
#       y = reshape(x,(npix,nepochs))
#       temporalf = sum( (y[:,2:end]-y[:,1:end-1]).^2 )
#       tv_g = Array{Float64}(undef, npix,nepochs)
#       if nepochs>2
#          tv_g[:,1] = 2*(y[:,1] - y[:,2])
#          tv_g[:,2:end-1] = 4*y[:,2:end-1]-2*(y[:,1:end-2]+y[:,3:end])
#          tv_g[:,end] = 2*(y[:,end] - y[:,end-1])
#       else
#          tv_g[:,1] = 2*(y[:,1]-y[:,2]);
#          tv_g[:,2] = 2*(y[:,2]-y[:,1]);
#       end
#       chi2_f += regularizers[nepochs+1][1][2]*temporalf
#       g[:] += regularizers[nepochs+1][1][2]*vec(tv_g);
#       if verbose == true
#            printstyled("Temporal regularization: $temporalf\n", color=:yellow)
#       end
#     end
#   end
#   reg_f = spheroid_regularization(x, g, regularizers=regularizers[1:nepochs], verbose = verbose); # Note: adds to g
#   return chi2_f + reg_f;
# end

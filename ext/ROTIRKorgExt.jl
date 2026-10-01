# Korg.jl as a specific-intensity grid builder for ROTIR.
#
#     using ROTIR, Korg
#     grid = build_korg_grid(Teff = (6000, 8000, 9), logg = (3.0, 4.5, 4))
#     prov = TabulatedProvider(grid)
#
# Korg is a weak dependency because it is heavy: a MARCS atmosphere artifact plus a
# ~42 000-line VALD linelist. Making it weak keeps `using ROTIR` fast, and — the point of
# building a grid at all — once the grid is saved with `save_intensity_grid` nothing needs
# Korg again. Evaluation is `TabulatedProvider`'s job, which is pure ROTIR.
#
# Structured like ROTIRNautilusExt: the stubs are declared in src/intensity_provider.jl so
# that `build_korg_grid` can be named, and documented, before Korg is loaded.
#
# WHAT WAS MEASURED, not assumed (Korg 1.3.1, verified against the installed source):
#
#   * `interpolate_marcs` spans Teff ∈ [2800, 8000] K, logg ∈ [-0.5, 5.5]. Outside that it
#     throws `LazyMultilinearInterpError`. beta Cas (7208 K pole) is inside; Vega, Regulus,
#     beta Lyr and Spica are not, and need a tabulated Kurucz/TLUSTY grid instead.
#   * `I_scheme = "linear"` is MANDATORY for per-mu intensity. The default
#     `"linear_flux_only"` ignores `mu_values` entirely and returns a `(1, nlambda)` array.
#   * with it, `result.intensity` is a 3-TENSOR `(n_mu, n_lambda, n_layers)` — not the
#     matrix the published docs imply — and **layer 1 is the emergent surface**. Confirmed
#     physically: at layer 1, I(mu=1)/I(mu=0.001) = 1.86 and rises monotonically with mu,
#     while at the deepest layer the ratio is 1.000, i.e. the radiation field is isotropic
#     at depth, as it must be. Reading the wrong end of that axis yields a limb-darkening
#     law that is flat and plausible-looking.
#   * an INTEGER `mu_values` gives Gauss-Legendre nodes that exclude both mu = 1 and
#     mu -> 0 (for 5 points: 0.047 … 0.953). Interferometry needs the disc centre and the
#     extreme limb, so this passes an explicit mu vector, which Korg honours exactly.
#   * cost is negligible: a 200 A window at 6 mu takes 0.29 s after warm-up, so a
#     12 x 6 (Teff, logg) grid is ~25 s.
#   * `result.intensity`'s FIRST axis is rays, ordered INWARD THEN OUTWARD:
#     `RadiativeTransfer.radiative_transfer` allocates
#     `n_inward_rays + length(mu_surface_grid)` rows and fills the inward ones first. So the
#     emergent intensity is the **LAST nmu rows**, in the requested mu order, and taking the
#     first nmu instead yields mostly zeros (inward rays that never reach the surface layer)
#     that are still monotone in mu and so look plausible. For a PlanarAtmosphere
#     `n_inward_rays = 0` and the two coincide, which is exactly how such a mistake survives
#     testing on dwarfs and corrupts every giant.
#   * `interpolate_marcs` silently switches geometry at logg = 3.5 — `ShellAtmosphere` below,
#     `PlanarAtmosphere` above — and the two give QUALITATIVELY different mu dependence, so
#     a grid must not mix them. See `spherical` in `build_korg_grid` for the choice and why
#     it is the caller's.
module ROTIRKorgExt

using ROTIR
using Korg
using Printf

import ROTIR: build_korg_grid, korg_provider, RectGrid4, TabulatedProvider

# Korg works in Angstrom; ROTIR's `band` is in metres (see `band_of`). Convert once, here,
# so only one convention is ever live inside a RectGrid4.
const _ANGSTROM_M = 1e-10

_axis(spec::Tuple{<:Real,<:Real,<:Integer}) = collect(range(float(spec[1]), float(spec[2]),
                                                            length = spec[3]))
_axis(v::AbstractVector) = collect(float.(v))

"""
    build_korg_grid(; Teff, logg, μ, λ, linelist, A_X, vmic, continuum_only, verbose)

Synthesize a `(Teff, logg, μ, λ)` specific-intensity grid.

`Teff`, `logg` and `μ` are each either `(lo, hi, n)` or an explicit vector. `λ` is
`(λ_start_Å, λ_stop_Å)` plus an optional step in Å.

## `spherical`: plane-parallel or extended, and why it is your choice

`interpolate_marcs` switches geometry at `logg = 3.5` on its own. This never uses that
default, because the two give qualitatively different μ dependence and a grid must not mix
them inside one interpolation:

| logg | spherical `I(μ)/I(1)` at μ = .001, .1, .3, .6, .9, 1 | planar `I(1)/I(.001)` |
|---|---|---|
| 3.4 | 0, 0.75, 0.86, 0.94, 0.99, 1 | 1.9 |
| 1.0 | 0, 0.0002, 0.44, 0.78, 0.96, 1 | 6.4 |
| 0.5 | 0, 0, 0.07, 0.70, 0.94, 1 | 12.2 |

A spherical model's intensity **genuinely vanishes** below a μ_min that moves with
(Teff, logg): those rays pass outside the extended photosphere altogether. That is the sharp
spherical limb, and it is real — it is what makes spherical limb-darkening tables differ
from plane-parallel ones at all.

`spherical = false` (the default) is what is **self-consistent with a tessellated Roche
surface**. ROTIR's mesh already carries the global geometry: each tessel is a local patch
with its own `(Teff, logg, μ)`, and what belongs there is the intensity emerging at angle
`acos(μ)` from a locally plane-parallel atmosphere. A 1-D *spherical* structure cannot be
pasted onto a rotationally distorted surface without double-counting extension — its
extension is spherical, the surface is not. This is also, in effect, what PMOIRED does:
`rotastar.py:271-273` rescales `newMu = (μ − μ₀)/(1 − μ₀)` with
`μ₀ = √(1 − (d_inner/d_outer)²)`, which strips the spherical extension back out of a SATLAS
model to recover a plane-parallel-like law.

`spherical = true` is the right choice for a single, genuinely extended star — a ρ Cas or
RW Cep near `logg = 0.5`, where plane-parallel misses the extension that interferometry can
actually see. Two things then matter. The μ zero-region is data, not padding, so keep the
`logg` axis fine enough that interpolation does not smear the edge (it moves with `logg`);
and ROTIR's μ is measured at the Roche photosphere while the model's is at the outermost
shell radius, so the two reference radii differ by the atmospheric extension.

`continuum_only = true` runs with an empty linelist and `hydrogen_lines = false`, giving
the continuum intensity on the same axes — that is Korg's only route to a per-μ continuum,
since `result.cntm` is a disk-integrated flux rather than an intensity. Useful for
normalising a line profile; not needed for interferometry, which wants the absolute
intensity.

Returns a [`RectGrid4`](@ref) with `λ` in **metres** and intensity in Korg's units
(erg/s/cm²/Å/sr).
"""
function build_korg_grid(; Teff = (4000.0, 7500.0, 8),
                           logg = (2.5, 4.5, 5),
                           μ = [0.001, 0.05, 0.1, 0.2, 0.3, 0.45, 0.6, 0.75, 0.9, 1.0],
                           λ = (6540.0, 6580.0),
                           λ_step::Real = 0.05,
                           linelist = nothing,
                           A_X = Korg.format_A_X(0),
                           vmic::Real = 1.0,
                           spherical::Bool = false,
                           continuum_only::Bool = false,
                           verbose::Bool = true)
    Ta = _axis(Teff); ga = _axis(logg); ma = _axis(μ)

    # A mu = 0 node makes the outermost interpolation cell blend toward a value no
    # atmosphere produced — the defect in src/di.jl:260-262. RectGrid4 refuses it, but
    # failing here names the actual cause.
    first(ma) <= 0 && throw(ArgumentError(
        "build_korg_grid: the μ axis starts at $(first(ma)); use ~1e-3, not 0. A μ = 0 ray " *
        "is tangent to the surface and carries no emergent intensity, so the node would be " *
        "an artificial zero. (For `spherical = true` the intensity is legitimately zero " *
        "over a RANGE of small μ — that is the extended limb and it is data; this guard is " *
        "about the exactly-tangent ray, which is not.)"))
    issorted(ma) || throw(ArgumentError("build_korg_grid: μ must be increasing"))

    lo, hi = float(λ[1]), float(λ[2])
    ll = continuum_only ? Korg.Line[] :
         (linelist === nothing ? Korg.get_VALD_solar_linelist() : linelist)

    # One synthesis fixes the λ sampling for every node, so take the axis from the first.
    verbose && @printf("build_korg_grid: %d Teff x %d logg x %d μ over %.1f-%.1f Å (%s)\n",
                       length(Ta), length(ga), length(ma), lo, hi,
                       spherical ? "spherical" : "plane-parallel")
    t0 = time()
    first_res = _synth(ll, A_X, Ta[1], ga[1], ma, lo, hi, λ_step, vmic, continuum_only, spherical)
    λÅ = collect(first_res.wavelengths)
    nλ = length(λÅ)
    values = Array{Float64,4}(undef, length(Ta), length(ga), length(ma), nλ)
    values[1, 1, :, :] = _surface_intensity(first_res, length(ma), nλ)

    nnode = length(Ta) * length(ga)
    done = 1
    for (i, tk) in enumerate(Ta), (j, lg) in enumerate(ga)
        (i == 1 && j == 1) && continue
        res = _synth(ll, A_X, tk, lg, ma, lo, hi, λ_step, vmic, continuum_only, spherical)
        length(res.wavelengths) == nλ || error(
            "build_korg_grid: node (Teff=$tk, logg=$lg) returned $(length(res.wavelengths)) " *
            "wavelengths, the first node gave $nλ — the λ sampling must be identical " *
            "across nodes for a rectilinear grid")
        values[i, j, :, :] = _surface_intensity(res, length(ma), nλ)
        done += 1
        verbose && @printf("\r  %d/%d nodes (%.1f s)", done, nnode, time() - t0)
    end
    verbose && @printf("\r  %d/%d nodes in %.1f s\n", nnode, nnode, time() - t0)

    return RectGrid4(Ta, ga, ma, λÅ .* _ANGSTROM_M, values)
end

function _synth(ll, A_X, Teff, logg, mu, lo, hi, step, vmic, continuum_only, spherical)
    # ONE geometry for the whole grid, chosen by the caller — never Korg's logg-dependent
    # default, which would mix planar and shell nodes inside a single interpolation.
    atm = Korg.interpolate_marcs(Teff, logg, A_X; spherical = spherical)
    # Wavelengths as ONE TUPLE, and the step is its THIRD ELEMENT — `(start, stop, step)`,
    # per Korg.Wavelengths' constructor. There is no `wavelength_step` keyword. Passing the
    # bounds as separate positional arguments is also deprecated in 1.3 and would warn on
    # every one of the (potentially hundreds of) nodes.
    return Korg.synthesize(atm, ll, A_X, (lo, hi, step);
                           mu_values = mu, I_scheme = "linear",
                           hydrogen_lines = !continuum_only, return_cntm = false,
                           vmic = vmic)
end

# The emergent intensity is layer 1 of the third axis; see the header for how that was
# established and why reading the wrong end is a silent, physical-looking error.
function _surface_intensity(res, nμ, nλ)
    I = res.intensity
    ndims(I) == 3 || error(
        "ROTIRKorgExt: expected a (n_mu, n_lambda, n_layers) intensity tensor, got " *
        "$(size(I)). `I_scheme = \"linear\"` is required — the default " *
        "\"linear_flux_only\" ignores `mu_values` and returns one row.")
    size(I, 1) >= nμ || error(
        "ROTIRKorgExt: Korg returned only $(size(I,1)) ray rows for $nμ requested μ values")
    size(I, 2) == nλ || error("ROTIRKorgExt: Korg returned $(size(I,2)) λ, asked $nλ")
    # The LAST nmu rows are the outward (emergent) rays; see the header. `n_inward` is 0 for
    # a planar atmosphere and grows with atmospheric extension for a shell.
    return I[end-nμ+1:end, :, 1]
end

korg_provider(; strict::Bool = false, kwargs...) =
    TabulatedProvider(build_korg_grid(; kwargs...); strict = strict, name = "Korg/MARCS")

end # module

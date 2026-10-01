#!/usr/bin/env julia
# =======================================================================================
# A hot-star intensity grid from Kurucz ATLAS9 — the stars Korg cannot reach
# =======================================================================================
# Korg's `interpolate_marcs` stops at Teff = 8000 K. That covers beta Cas (7208 K) and nothing
# hotter, so Spica (25300 / 20585 K), Vega (~9500 K), Regulus (~15000 K) and beta Lyr all need
# a different source. Kurucz's intensity packs run **3500 to 50000 K at 17 emergent angles**,
# which covers every one of them with a single file.
#
# WHERE TO GET IT. http://kurucz.harvard.edu/grids/gridp00/ holds two solar-metallicity packs
# and the choice matters more than the names suggest:
#
#     ip00k2.pck19    67.5 MB   Teff  3500-50000 K   <- the one you want
#     ip00k2new.pck   31.5 MB   Teff  3500- 8750 K   <- newer opacities, cool stars only
#
# `ip00k2new.pck` is the newer ODFNEW grid but tops out at 8750 K, i.e. no better than Korg
# for this purpose. Other metallicities are in the sibling `gridm05`, `gridm10`, … directories
# under the same naming (`im05k2.pck19`). Nothing is bundled with ROTIR: these are third-party
# data, and `read_kurucz_intensity` converts whichever you fetch.
#
#     curl -O http://kurucz.harvard.edu/grids/gridp00/ip00k2.pck19
#     julia --project=demos demos/spica_kurucz_grid.jl ip00k2.pck19
#
# ON TLUSTY. The plan called for TLUSTY OSTAR2002/BSTAR2006 as well, for NLTE above ~30000 K.
# Its canonical host is GONE: nova.astro.umd.edu now 301-redirects to tlusty.oca.eu, which
# serves an unrelated conference page, and the Arizona mirror 404s. The Spanish VO
# (svo2.cab.inta-csic.es) and CDS mirror TLUSTY *fluxes*, not per-mu intensities, which are
# what a limb-darkened surface needs. So TLUSTY is not currently fetchable and this demo does
# not pretend otherwise. It matters less than it sounds: Kurucz is LTE, and LTE is adequate
# through Spica's 25300 K — NLTE becomes important for O stars, which ROTIR has no target for
# yet. If you have the intensity files locally, a reader for them belongs beside
# `read_kurucz_intensity` and `load_intensity_grid` is the seam.
# =======================================================================================

using ROTIR, Printf

const PACK = length(ARGS) >= 1 ? ARGS[1] :
             get(ENV, "KURUCZ_PACK", joinpath(@__DIR__, "data", "ip00k2.pck19"))

if !isfile(PACK)
    println("""
    Kurucz pack not found: $(PACK)

    Fetch one (67.5 MB) and pass it as an argument:
        curl -O http://kurucz.harvard.edu/grids/gridp00/ip00k2.pck19
        julia --project=demos demos/spica_kurucz_grid.jl ip00k2.pck19
    """)
    exit(0)
end

println("="^88)
println("Kurucz ATLAS9 intensity pack -> a ROTIR grid")
println("="^88)

# H band, where MIRC-X and GRAVITY observe. Restricting lambda at READ time matters: a full
# pack is 1221 wavelengths x 17 angles x ~476 models, and an H-band grid wants ~35 of them.
const HBAND = (1.50e-6, 1.75e-6)

all_models = read_kurucz_models(PACK; λrange = HBAND, verbose = false)
Ts = sort(unique(m.Teff for m in all_models)); gs = sort(unique(m.logg for m in all_models))
@printf("  %d models; Teff %.0f-%.0f K, logg %.1f-%.1f, %d wavelengths in H\n",
        length(all_models), minimum(Ts), maximum(Ts), minimum(gs), maximum(gs),
        length(first(all_models).λ))

# Spica A: B1 III-IV, Teff ~ 25300 K. Take a window around it rather than the whole pack —
# the grid only needs to span the SURFACE, and a gravity-darkened surface spans a few hundred
# K and a few tenths of a dex.
println("\n--- a grid for Spica's primary (Teff ~ 25300 K) ---")
g = read_kurucz_intensity(PACK; λrange = HBAND, Teff = (17000.0, 30000.0),
                          logg = (3.0, 4.5), verbose = false)
@printf("  grid %s   mu ascending: %s   lambda %.4f-%.4f um\n",
        string(size(g)), issorted(g.μ), g.λ[1]*1e6, g.λ[end]*1e6)

mid = size(g, 4) ÷ 2
hot  = g.values[findfirst(>=(25000.0), g.Teff), 2, :, mid]
@printf("\n  Teff = 25000 K, logg = %.1f:  I(mu)/I(1) = %s\n", g.logg[2],
        string(round.(hot ./ hot[end], digits = 3)))
@printf("  I(1)/I(0.01) = %.3f\n", hot[end] / hot[1])
# Computed, not hardcoded: the number depends on which (Teff, logg) node this lands on, and a
# literal in the prose would drift from it the first time the window changed.
@printf("""
  A HOT star is LESS limb-darkened in H than a warm one: %.2f here against 1.95 for beta Cas
  at 7000 K (and 1.95 is itself what the Korg H-band grid predicts for it).
""", hot[end] / hot[1])
println("""  That is the Rayleigh-Jeans limit asserting itself — at 1.65 um a 25000 K
  photosphere is far down the tail of the Planck function, where the source function varies
  slowly with depth, so the emergent intensity depends weakly on the viewing angle. It is also
  why fitting limb darkening on a hot star is even less well constrained than on beta Cas, and
  why PREDICTING it matters more there, not less.""")

out = joinpath(@__DIR__, "data", "spica_H_kurucz.fits")
save_intensity_grid(out, g; comment = "Kurucz ATLAS9 $(basename(PACK)), H band, Spica range")
@printf("\n  saved %s (%.2f MB) — readable with `load_intensity_grid`, no Kurucz file needed\n",
        basename(out), filesize(out)/1e6)

# And check it against a real gravity-darkened Spica-like surface.
println("\n--- does it cover a Spica-like surface? ---")
sp = default_star_params(2; rpole = 0.45, d = 80.0, frac_escapevel = 0.5,
                         rotation_period = 4.0, tpole = 25300.0, inclination = 60.0,
                         position_angle = 0.0, beta = 0.25, gravity_law = 2, ldtype = 0)
tess = tessellation_healpix(3, T = Float64)
st = create_star(tess, sp, 0.0)
Tm = parametric_temperature_map(sp, st)
θ  = tess.unit_spherical[:, 5, 2]
lg = logg_map(sp.rpole, sp.d, sp.frac_escapevel, sp.rotation_period, sin.(θ), cos.(θ))
@printf("  surface spans Teff %.0f-%.0f K, logg %.2f-%.2f\n", extrema(Tm)..., extrema(lg)...)
@printf("  inside the grid: Teff %s, logg %s\n",
        all(g.Teff[1] .<= Tm .<= g.Teff[end]), all(g.logg[1] .<= lg .<= g.logg[end]))
prov = TabulatedProvider(g; name = "Kurucz ATLAS9 H")
println("  provider consistency (ldtype = 0): ",
        isempty(check_provider_consistency(prov, sp)) ? "ok" :
        check_provider_consistency(prov, sp))
q = derived_quantities(sp)
@printf("  derived: M = %.2f Msun, R_p = %.2f Rsun, logg_p = %.2f, vsini = %.0f km/s\n",
        q.mass, q.rpole_rsun, q.logg_pole, q.vsini)
println("\n  (those Spica parameters are illustrative, not a fit — Spica is a BINARY and")
println("   belongs on surface_type 3, which this grid serves equally well.)")

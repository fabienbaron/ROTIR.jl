#!/usr/bin/env julia
# The Korg.jl grid builder (ext/ROTIRKorgExt.jl) against a REAL synthesised grid.
#
# OPT-IN: needs Korg, which pulls a MARCS atmosphere artifact and a ~42 000-line VALD
# linelist, and ~30 s of synthesis. Korg is deliberately NOT in `[targets] test`.
#
#     ROTIR_TEST_KORG=1 julia --project=<env-with-Korg> test/test_korg_grid.jl
#
# WHAT THIS GUARDS, and none of it is hypothetical - every one was a live bug during
# development:
#
#  * `result.intensity`'s first axis is rays ordered INWARD then OUTWARD, so the emergent
#    intensity is the LAST nmu rows. Taking the first nmu gives mostly zeros that are still
#    monotone in mu, hence plausible; and for a PlanarAtmosphere n_inward = 0 so the two
#    coincide, which is how the mistake survives testing on dwarfs and corrupts every giant.
#  * `I_scheme = "linear"` is mandatory; the default "linear_flux_only" ignores `mu_values`
#    and returns one row.
#  * Korg returns VACUUM wavelengths. Halpha is 6562.79 A in air and 6564.60 in vacuum - a
#    1.8 A offset, which is 3e-4, i.e. 83 km/s. The line-position check below asserts the
#    vacuum value on purpose.
#  * EXTNAME must be set with `write_key` AFTER the HDU exists; passing it inside a
#    FITSHeader to `write` is silently dropped, and the file then has no axis labels at all.

using ROTIR, Printf, LinearAlgebra

if get(ENV, "ROTIR_TEST_KORG", "0") != "1"
    @info "test_korg_grid.jl skipped (set ROTIR_TEST_KORG=1 to run)"
    npass = Ref(0); nfail = Ref(0); korg_ran = Ref(false)
else
    korg_ran = Ref(true)
    using Korg
    npass = Ref(0); nfail = Ref(0)
    function cb(label, ok)
        ok ? (npass[] += 1) : (nfail[] += 1)
        @printf("  %-58s %s\n", label, ok ? "\u2713" : "\u2717")
        return ok
    end

    println("\n[1] build a small real grid (Halpha region, beta Cas-ish)")
    t0=time()
    g = build_korg_grid(Teff=(6600.0,7400.0,3), logg=(3.4,4.0,2),
                        μ=[0.001,0.1,0.3,0.6,0.9,1.0], λ=(6555.0,6570.0), λ_step=0.05)
    @printf("  built in %.1f s: %s\n", time()-t0, size(g))
    cb("axes match values", size(g.values)==(3,2,6,length(g.λ)))
    cb("lambda axis in METRES", 6.5e-7 < g.λ[1] < 6.6e-7)
    cb("all intensities finite and positive", all(isfinite,g.values) && all(>(0),g.values))

    println("\n[2] physics of the grid")
    mid = length(g.λ)÷2
    # Limb darkening: I must increase with mu at every (Teff, logg).
    mono = all(issorted(g.values[i,j,:,mid]) for i in 1:3, j in 1:2)
    cb("I increases with mu everywhere (limb darkening)", mono)
    r = g.values[2,1,end,mid]/g.values[2,1,1,mid]
    @printf("      I(mu=1)/I(mu=0.001) = %.3f at Teff=%.0f logg=%.1f\n", r, g.Teff[2], g.logg[1])
    cb("limb-darkening contrast is physical (1.2 < r < 4)", 1.2 < r < 4.0)
    # Hotter -> brighter, at fixed mu and lambda.
    cb("I increases with Teff", all(issorted(g.values[:,j,k,mid]) for j in 1:2, k in 1:6))
    # Halpha absorption: the line core must be fainter than the window edge.
    core = argmin(g.values[2,1,end,:])
    # Korg returns VACUUM wavelengths (air_wavelengths=false is its default). Halpha is
    # 6562.79 A in air, 6564.60 A in vacuum - a 1.8 A offset, which at 6563 A is 3e-4, i.e.
    # 83 km/s. Asserting the air value here is what first made this check "fail".
    @printf("      deepest lambda = %.3f A (Halpha: 6562.79 air, 6564.60 vacuum)\n", g.λ[core]*1e10)
    cb("absorption line at Halpha VACUUM (6564.60)", abs(g.λ[core]*1e10 - 6564.60) < 0.3)
    cb("line core is fainter than continuum", g.values[2,1,end,core] < g.values[2,1,end,1])

    println("\n[3] FITS round trip (no JLD2, no Korg needed to read it back)")
    f = tempname()*".fits"
    save_intensity_grid(f, g; comment="test grid")
    h = load_intensity_grid(f)
    cb("values identical", h.values == g.values)
    cb("all four axes identical", h.Teff==g.Teff && h.logg==g.logg && h.μ==g.μ && h.λ==g.λ)
    @printf("      file size %.2f MB\n", filesize(f)/1e6)
    rm(f)

    println("\n[4] as a provider, on a real gravity-darkened surface")
    sp = default_star_params(2; rpole=0.849, d=16.8, frac_escapevel=0.92, rotation_period=1/1.12,
                             tpole=7208.0, inclination=19.9, position_angle=-7.09,
                             beta=0.25, ldtype=0, gravity_law=2)
    tess = tessellation_healpix(3, T=Float64)
    st = create_star(tess, sp, 0.0)
    Tmap = parametric_temperature_map(sp, st)
    lgm  = logg_map(sp.rpole, sp.d, sp.frac_escapevel, sp.rotation_period,
                    sin.(tess.unit_spherical[:,5,2]), cos.(tess.unit_spherical[:,5,2]))
    mu   = max.(st.normals[:,3], 0.0)
    @printf("      Teff %.0f-%.0f K, logg %.2f-%.2f, both inside the grid? ",
            extrema(Tmap)..., extrema(lgm)...)
    inside = all(g.Teff[1] .<= Tmap .<= g.Teff[end]) && all(g.logg[1] .<= lgm .<= g.logg[end])
    println(inside)
    prov = TabulatedProvider(g)
    I = provider_intensity(prov, Tmap, lgm, mu, g.λ[mid])
    cb("intensity evaluated on the mesh", length(I)==tess.npix && all(isfinite,I))
    vis = st.index_quads_visible
    cb("visible-side intensity all positive", all(>(0), I[vis]))
    cb("ldtype=0 + provider is consistent", isempty(check_provider_consistency(prov,sp)))
    # The whole point: the pole is hotter AND nearer disc-centre, so it must be brightest.
    @printf("      I spans %.3g - %.3g over the visible disc (ratio %.2f)\n",
            extrema(I[vis])..., maximum(I[vis])/minimum(I[vis]))
    cb("brightness varies across the disc", maximum(I[vis])/minimum(I[vis]) > 1.1)

    println("\n[5] spherical vs plane-parallel: both geometries, and the row extraction")
    # THE check the design turns on. `interpolate_marcs` switches geometry at logg = 3.5 on
    # its own, and a ShellAtmosphere emits extra INWARD ray rows ahead of the outward ones.
    # If the builder took the first nmu rows instead of the last, planar would still look
    # right (n_inward = 0) and spherical would be mostly zeros. So: build the same node both
    # ways and require each to be a sane limb-darkening law.
    gp = build_korg_grid(Teff=(5000.0,5200.0,2), logg=(1.0,1.5,2), μ=[0.001,0.1,0.3,0.6,0.9,1.0],
                         λ=(6560.0,6566.0), λ_step=0.2, spherical=false, verbose=false)
    gs = build_korg_grid(Teff=(5000.0,5200.0,2), logg=(1.0,1.5,2), μ=[0.001,0.1,0.3,0.6,0.9,1.0],
                         λ=(6560.0,6566.0), λ_step=0.2, spherical=true,  verbose=false)
    mid5 = length(gp.λ)÷2
    cb("planar and spherical grids have identical axes",
       gp.Teff==gs.Teff && gp.logg==gs.logg && gp.μ==gs.μ && gp.λ==gs.λ)
    cb("planar: I increases with mu at every node",
       all(issorted(gp.values[i,j,:,mid5]) for i in 1:2, j in 1:2))
    cb("spherical: I increases with mu at every node",
       all(issorted(gs.values[i,j,:,mid5]) for i in 1:2, j in 1:2))
    # If the FIRST rows had been taken, the disc centre (mu=1) would not be the brightest.
    cb("planar: disc centre is brightest", all(argmax(gp.values[i,j,:,mid5])==length(gp.μ) for i in 1:2, j in 1:2))
    cb("spherical: disc centre is brightest", all(argmax(gs.values[i,j,:,mid5])==length(gs.μ) for i in 1:2, j in 1:2))
    cb("both have a nonzero disc centre",
       all(>(0), gp.values[:,:,end,mid5]) && all(>(0), gs.values[:,:,end,mid5]))
    # The physics that makes the choice matter: a spherical model's limb is SHARP - the
    # intensity falls to zero over a range of small mu, because those rays miss the extended
    # photosphere entirely. Planar has no such edge. This is the difference the `spherical`
    # keyword exists to expose, not a defect.
    pedge = gp.values[1,1,1,mid5] / gp.values[1,1,end,mid5]
    sedge = gs.values[1,1,1,mid5] / gs.values[1,1,end,mid5]
    @printf("      I(mu=.001)/I(mu=1) at Teff=5000 logg=1.0:  planar %.4f   spherical %.4f\n",
            pedge, sedge)
    cb("planar limb fades gradually (I(.001)/I(1) > 0.01)", pedge > 0.01)
    cb("spherical limb is sharp (I(.001)/I(1) < planar)", sedge < pedge)

    @printf("\n=== %d passed, %d failed ===\n", npass[], nfail[])
end

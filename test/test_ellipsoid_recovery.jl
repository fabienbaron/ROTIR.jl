# Does `fit_ellipsoid_ld` recover a KNOWN triaxial figure?
#
# WHY THIS EXISTS. Until 2026-10-03 the function imposed `radius_x = radius_y = req` and
# `radius_z = req*(1-f)` — a 5-parameter OBLATE figure, a two-dimensional slice of the model's
# three-radius shape space, with the symmetry axis pinned to the spin axis. A tidally distorted
# star is elongated towards its companion (`radius_x > radius_y ≈ radius_z`), which that slice
# cannot represent at any parameter value, so on ζ And it drove `f` to its bound at 0.3997 —
# past the 1/3 a rigid rotator reaches at break-up — and the number was read as a shape.
#
# Nothing caught it, because nothing had ever checked that this function recovers a figure it
# was given. The gradient tests (`test_ellipsoid_logpi.jl`) check the log-posterior's
# derivatives; the sampler test checks `isfinite`. Neither can see a wrong shape family.
#
# The generator and the fit use the SAME mesh on purpose: discretisation then cancels exactly
# and what is left is the fitter. Mesh bias is a separate question, measured in
# `demos/zetand_triaxial.jl` by generating at a finer mesh than it fits.

using Test
using ROTIR

if !@isdefined(PKGDIR)
    const PKGDIR = pkgdir(ROTIR)
end

@testset "the triaxial ellipsoid recovers a known figure" begin

ZETDIR = joinpath(PKGDIR, "demos", "zetAnd")
PROT   = 17.769426        # d — ζ And is tidally locked (Roettenbacher et al. 2016)

"Real ζ And uv coverage and epoch times; the observables get overwritten with the model."
function zetdata()
    files = sort(filter(f -> endswith(f, ".oifits"), readdir(ZETDIR)))
    d = [readoifits(joinpath(ZETDIR, f); verbose = false)[1] for f in files]
    m = [minimum(x.v2_mjd) for x in d]
    (d, m .- minimum(m))
end

"Overwrite `data`'s observables with this figure's own forward model, uniform surface."
function synth!(data, tess, truth, ts)
    stars = create_star_multiepochs(tess, truth, ts)
    x = fill(truth.tpole, stars[1].npix)
    for (i, s) in enumerate(stars)
        v2m, t3am, t3pm = ROTIR.cvis_to_obs(ROTIR.fused_cvis(x, s, data[i]), data[i])
        data[i].v2 .= v2m; data[i].t3amp .= t3am; data[i].t3phi .= t3pm
    end
end

# A ζ And-like figure: elongated towards the companion (x), flattened at the pole.
TRUTH = (surface_type = 1, radius_x = 1.350, radius_y = 1.310, radius_z = 1.255,
         tpole = 5000.0, ldtype = 3, ld1 = 0.20, ld2 = 0.0, beta = 0.25,
         inclination = 65.0, position_angle = 40.0, rotation_period = PROT)

@testset "cold start, three rotation phases" begin
    tess = tessellation_healpix(2)
    data, ts = zetdata()
    keep = [3, 5, 4]                      # phases 0.00, 0.50, 0.17 — well separated
    data = data[keep]; ts = ts[keep]
    synth!(data, tess, TRUTH, ts)
    # COLD: the fit is told the period and nothing else. `seed_from_sphere` finds the radii.
    r = fit_ellipsoid_ld(data, tess; tepochs = ts, rotation_period = PROT, ldtype = 3,
                         r_bounds = (0.5, 3.0), restarts = 4)
    for (name, got, want) in (("radius_x", r.radius_x, TRUTH.radius_x),
                              ("radius_y", r.radius_y, TRUTH.radius_y),
                              ("radius_z", r.radius_z, TRUTH.radius_z))
        @test abs(got/want - 1) < 1e-3
    end
    @test abs(r.inclination - TRUTH.inclination) < 0.5
    @test abs(r.position_angle - TRUTH.position_angle) < 0.5
    @test abs(r.ld1 - TRUTH.ld1) < 5e-3
    # THE NUMBER THE OLD PARAMETERIZATION COULD NOT PRODUCE: an equatorial elongation. The
    # oblate slice forced this to be exactly 1 by construction, whatever the data said.
    @test abs(r.equatorial_elongation - TRUTH.radius_x/TRUTH.radius_y) < 1e-3
    @test r.equatorial_elongation > 1.01        # and it is resolved, not a rounding artefact
end

@testset "one rotation phase cannot determine three radii" begin
    # The projected outline is ONE ellipse — three numbers — against five shape and
    # orientation parameters. The fit must still run and must say so.
    tess = tessellation_healpix(2)
    data, _ = zetdata()
    data = data[1:3]; ts = zeros(3)
    synth!(data, tess, TRUTH, ts)
    run1() = fit_ellipsoid_ld(data, tess; tepochs = ts, rotation_period = PROT, ldtype = 3,
                              r_bounds = (0.5, 3.0), restarts = 1)
    r = @test_logs (:warn, r"same rotation phase") match_mode=:any run1()
    @test length(r.phases) == 1
end

@testset "the oblate keywords are refused, not reinterpreted" begin
    tess = tessellation_healpix(2)
    data, ts = zetdata()
    for kw in (:req0 => 1.3, :f0 => 0.05, :req_bounds => (0.5, 3.0), :f_bounds => (0.0, 0.4))
        @test_throws ErrorException fit_ellipsoid_ld(data[1:1], tess; tepochs = ts[1:1],
                                                     rotation_period = PROT, (kw,)...)
    end
end

@testset "the period is required when the epochs differ" begin
    tess = tessellation_healpix(2)
    data, ts = zetdata()
    # Multi-epoch without a period: a triaxial figure's projection depends on the phase, so
    # there is no safe placeholder — and the old code's hardcoded 1.0 would have turned a
    # 20-day baseline into 20 turns.
    @test_throws ErrorException fit_ellipsoid_ld(data[1:2], tess; tepochs = ts[1:2],
                                                 ldtype = 3)
    # One instant: the phase cannot matter, so it is allowed.
    r = fit_ellipsoid_ld(data[1:1], tess; tepochs = [0.0], ldtype = 3,
                         r_bounds = (0.5, 3.0), restarts = 1)
    @test r.rotation_period == 1.0
end

end

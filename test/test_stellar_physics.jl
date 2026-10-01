#!/usr/bin/env julia
# Absolute mass, radius and surface gravity for the rapid rotator (src/stellar_physics.jl).
#
# Standalone script in the style of test_parametric_gradient.jl: it prints its own table,
# never throws on a numerical mismatch, and exposes `nfail[]` for runtests.jl to assert on.
#
#     julia --project=demos test/test_stellar_physics.jl
#
# WHAT THIS IS GUARDING. `logg` is what lets a model atmosphere be indexed by (Teff, logg),
# which is what makes limb darkening PREDICTED instead of fitted. If the gravity or the
# derived mass is wrong, every intensity that follows is wrong in a way that looks like a
# plausible fit — so the checks here are against an INDEPENDENT quantity (beta Cas's
# published mass, radius, logg and vsini) and against finite differences, not against
# values this code produced.

using ROTIR, LinearAlgebra, Printf

const HAVE_FD = try
    @eval using FiniteDifferences
    true
catch
    false
end

const FDM = HAVE_FD ? central_fdm(5, 1) : nothing

relerr(a, b; scale = 0.0) = norm(a .- b) / max(norm(a), norm(b), scale, eps())
status(r) = r < 1e-6 ? "✓" : (r < 1e-3 ? "~" : "✗")

npass = Ref(0); nfail = Ref(0); fd_ran = Ref(false)

function check(label, r; tol = 1e-6)
    ok = r < tol
    ok ? (npass[] += 1) : (nfail[] += 1)
    @printf("  %-52s relerr=%.3e  %s\n", label, r, status(r))
    return ok
end
checkbool(label, ok) = check(label, ok ? 0.0 : 1.0)

# beta Cas as fitted in demos/rapid_rotator_betCas_param_fit.jl, plus the Gaia distance.
# Only (rpole, d, fev, P) are inputs; mass, radius, logg and vsini are all OUTPUTS, which
# is what makes the literature comparison below a real test rather than a restatement.
const RP, DPC, FEV, PROT, INC = 0.849, 16.8, 0.92, 1/1.12, 19.9

println("\n[1] beta Cas — derived quantities against the literature")
let q = derived_quantities((rpole = RP, d = DPC, frac_escapevel = FEV,
                            rotation_period = PROT, inclination = INC))
    @printf("      M = %.3f Msun   R_p = %.3f   R_eq = %.3f Rsun   logg_p = %.3f   vsini = %.1f km/s\n",
            q.mass, q.rpole_rsun, q.req_rsun, q.logg_pole, q.vsini)
    # Generous windows: these are order-of-magnitude-and-then-some sanity gates on an
    # F2III-IV star, not a measurement. They catch a wrong constant or a dropped 1e-3,
    # which is the failure mode that matters.
    checkbool("mass in [1.5, 2.5] Msun (lit. ~1.9)",       1.5 <= q.mass <= 2.5)
    checkbool("R_p in [2.5, 4.0] Rsun (lit. ~3.5)",        2.5 <= q.rpole_rsun <= 4.0)
    checkbool("R_eq > R_p (the star is oblate)",           q.req_rsun > q.rpole_rsun)
    checkbool("logg_pole in [3.3, 4.1] (lit. 3.5-3.8)",    3.3 <= q.logg_pole <= 4.1)
    checkbool("vsini in [55, 90] km/s (lit. ~70)",         55 <= q.vsini <= 90)
end

println("\n[2] the mass/period relation round-trips")
let M = derive_mass(RP, DPC, FEV, PROT)
    # Recover the period from the derived mass through Omega_crit — the inverse of the
    # relation derive_mass inverts. Agreement to round-off means the two agree on what
    # `frac_escapevel` MEANS, which is the invariant that used to be unenforced: nothing
    # tied rotation_period (read by rotate_star) to frac_escapevel (read by the maps).
    Rp = polar_radius_m(RP, DPC)
    P_back = 2π / (FEV * sqrt(8 * M * 1.3271244e20 / (27 * Rp^3))) / 86400
    check("period -> mass -> period", relerr(P_back, PROT); tol = 1e-12)
end

println("\n[3] logg map: structure")
let θ = collect(range(1e-3, π/2, length = 33))
    lg = logg_map(RP, DPC, FEV, PROT, sin.(θ), cos.(θ))
    checkbool("logg decreases monotonically pole -> equator", all(diff(lg) .< 0))
    checkbool("pole-to-equator span > 0.2 dex (a scalar-logg grid cannot serve this)",
              lg[1] - lg[end] > 0.2)
    check("logg_pole() == map at the pole", relerr(logg_pole(RP, DPC, FEV, PROT), lg[1]);
          tol = 1e-7)
end

println("\n[4] equator symmetry: logg(θ) == logg(π−θ)")
let θ = collect(range(1e-3, π - 1e-3, length = 41))
    lg = logg_map(RP, DPC, FEV, PROT, sin.(θ), cos.(θ))
    # The southern hemisphere is the northern one mirrored for a rigidly rotating star.
    # This is the check that caught the hemisphere-folding bug in the ELR flux factor
    # (see src/gravity_darkening.jl), where the analytic derivative made the SAME mistake
    # self-consistently and so a finite-difference check could not see it.
    check("north/south mirror", relerr(lg[1:20], reverse(lg[22:41])); tol = 1e-12)
end

println("\n[5] the non-rotating limit")
let θ = collect(range(1e-3, π/2, length = 17))
    # As fev -> 0 the surface becomes a sphere and the gravity uniform, so logg must go
    # flat. It cannot be evaluated AT fev = 0: a sphere spinning at a finite rate needs
    # infinite mass for that rate to be a vanishing fraction of critical, which is why
    # derive_mass diverges there rather than returning something plausible.
    lg = logg_map(RP, DPC, 1e-6, PROT, sin.(θ), cos.(θ))
    check("logg uniform at fev = 1e-6", (maximum(lg) - minimum(lg)) / abs(lg[1]); tol = 1e-9)
    checkbool("derive_mass diverges as fev -> 0", isinf(derive_mass(RP, DPC, 0.0, PROT)))
end

println("\n[6] closed forms for the mass log-derivatives")
let (_, l1, l2, l3, l4) = derive_mass_and_dlog(RP, DPC, FEV, PROT)
    # ln M = const + 3 ln R_p + 2 ln Omega - 2 ln fev, and R_p ∝ rpole·d, Omega ∝ 1/P.
    # Every derivative is therefore a bare power, exact to the last bit.
    check("dlnM/drpole == 3/rpole", relerr(l1,  3/RP);   tol = 1e-15)
    check("dlnM/dd     == 3/d",     relerr(l2,  3/DPC);  tol = 1e-15)
    check("dlnM/dfev   == -2/fev",  relerr(l3, -2/FEV);  tol = 1e-15)
    check("dlnM/dP     == -2/P",    relerr(l4, -2/PROT); tol = 1e-15)
end

println("\n[7] type stability: a Float32 model stays Float32")
let θ = collect(range(1e-3, π/2, length = 17))
    lg32 = logg_map(Float32(RP), Float32(DPC), Float32(FEV), Float32(PROT),
                    Float32.(sin.(θ)), Float32.(cos.(θ)))
    lg64 = logg_map(RP, DPC, FEV, PROT, sin.(θ), cos.(θ))
    checkbool("eltype(logg_map) === Float32", eltype(lg32) === Float32)
    check("Float32 agrees with Float64", relerr(Float64.(lg32), lg64); tol = 1e-5)
end

if HAVE_FD
    println("\n[8] analytic derivatives vs finite differences  —  requires FiniteDifferences")
    let θ = collect(range(1e-3, π - 1e-3, length = 41))
        s, c = sin.(θ), cos.(θ)
        base = [RP, DPC, FEV, PROT]
        ana = logg_map_and_derivs(base..., s, c)[2:5]
        for (k, nm) in enumerate(("rpole", "d", "fev", "rotation_period"))
            fd = FiniteDifferences.jacobian(
                     FDM, p -> logg_map(p[1], p[2], p[3], p[4], s, c), base)[1][:, k]
            check("dlogg/d$(nm)  (41 colatitudes)", relerr(ana[k], fd); tol = 1e-8)
        end
        for (k, nm) in enumerate(("rpole", "d", "fev", "rotation_period"))
            an = derive_mass_and_dlog(base...)[k+1]
            fd = FiniteDifferences.grad(
                     FDM, p -> log(derive_mass(p[1], p[2], p[3], p[4])), base)[1][k]
            check("dlnM/d$(nm)", relerr(an, fd); tol = 1e-9)
        end
        fd_ran[] = true
    end
else
    @warn "FiniteDifferences unavailable — section [8] skipped"
end

println("\n[9] ldtype = 0 is the identity, on BOTH the forward and the gradient path")
let nz = collect(range(-0.5, 1.0, length = 25))
    # The trap: ld_and_derivs' fallback branch is ldtype 3, so without an explicit ldtype
    # 0 branch the gradient path silently applies the Hestroffer law mu^ld1 while
    # compute_ldmap returns ones. The fit would then optimise a model nobody reports.
    fwd = ROTIR.compute_ldmap(max.(nz, 0.0), default_star_params(2; ldtype = 0, ld1 = 0.37))
    ad, dnz, dld1, dld2 = ROTIR.ld_and_derivs(nz, 0, 0.37, 0.1)
    checkbool("compute_ldmap(ldtype=0) is all ones", all(==(1.0), fwd))
    # Exact equality, expressed through checkbool rather than `check(..., tol = 0.0)`:
    # the criterion is `relerr < tol`, so a zero tolerance can never be met even by two
    # identical vectors.
    checkbool("ld_and_derivs matches it (not mu^ld1)", ad == fwd)
    checkbool("every LD derivative is zero",
              all(iszero, dnz) && all(iszero, dld1) && all(iszero, dld2))
    check("ldtype=3 still is mu^ld1",
          relerr(first(ROTIR.ld_and_derivs(nz, 3, 0.37, 0.1)), max.(nz, 0.0) .^ 0.37))
end

println("\n[10] the schema carries the distance, optionally")
let sp = surface_params(2)
    checkbool("`d` is in surface_params(2)", :d in [x.name for x in sp])
    checkbool("`d` is OPTIONAL, so older models still validate",
              :d in [x.name for x in surface_spec(2).optional])
    checkbool("`d` is grouped :geometry here, :orbit on the Roche surface",
              sp[findfirst(x -> x.name === :d, sp)].group === :geometry &&
              surface_params(3)[findfirst(x -> x.name === :d,
                                          surface_params(3))].group === :orbit)
    checkbool("ldtype offers 0", 0 in first.(sp[findfirst(x -> x.name === :ldtype, sp)].choices))
    checkbool("ld_coefficients_used(0) is empty", ld_coefficients_used(0) == ())
end

println("\n[11] an inconsistent rotation_period is reported")
let good = default_star_params(2; rpole = RP, d = DPC, frac_escapevel = FEV,
                               rotation_period = PROT)
    bad = merge(good, (rotation_period = 50.0,))
    # ADVISORY, not blocking: the GUI gates on `validate_star_params` being empty, so a
    # plausibility message there would disable `epoch_chi2` rather than inform anyone.
    checkbool("beta Cas draws no mass complaint",
              isempty(filter(m -> occursin("Msun", m), advise_star_params(good))))
    checkbool("a 50-day period at fev = 0.92 does",
              length(filter(m -> occursin("Msun", m), advise_star_params(bad))) == 1)
    checkbool("the mass advice does NOT block validate_star_params",
              isempty(filter(m -> occursin("Msun", m), validate_star_params(bad))))
end

println("\n[12] derived_summary_text — the GUI readout, tested without the GUI")
# Lives in the core precisely so it can be tested here: the GUI extension needs GLMakie and
# QML, which a headless test run does not have.
let good = default_star_params(2; rpole = RP, d = DPC, frac_escapevel = FEV,
                               rotation_period = PROT)
    t = derived_summary_text(good)
    checkbool("reports a mass", occursin("M = ", t))
    checkbool("reports both radii", occursin("R_p = ", t) && occursin("R_eq = ", t))
    checkbool("reports logg and both velocities",
              occursin("logg_p", t) && occursin("v_eq", t) && occursin("vsini", t))
    # The numbers must be the ones `derived_quantities` gives, not a second computation.
    q = derived_quantities(good)
    checkbool("the mass shown matches derived_quantities",
              occursin(string(round(q.mass, digits = 3)), t))
    # Empty, not an error, wherever there is nothing to say — that is what hides the label.
    checkbool("empty without a distance",
              isempty(derived_summary_text(Base.structdiff(good, NamedTuple{(:d,)}))))
    checkbool("empty for a sphere", isempty(derived_summary_text(default_star_params(0))))
    checkbool("empty for a Roche surface", isempty(derived_summary_text(default_star_params(3))))
    # The case it exists for: a fit that wandered to an unphysical mass is visible at a glance.
    runaway = merge(good, (rpole = 0.98, frac_escapevel = 0.59))
    checkbool("an implausible mass is plainly shown",
              occursin("M = 7.2", derived_summary_text(runaway)))
end

println("\n[13] a parameter that cannot describe a surface is BLOCKED, not advised")
# THE REGRESSION THIS PINS. `rpole = 0` takes a Roche equipotential through a zero radius and
# the temperature map comes back entirely NaN. Nothing downstream raises: `_map_range`'s
# `minimum` propagates NaN, its `pmax - pmin < 1.0` widening never fires because `NaN < 1.0`
# is false, and the NaN reaches a Makie Colorbar — whose tick machinery then throws
# `InexactError: convert(UInt64, …)` out of Ryu's `writefixed`, inside a QML callback, which
# freezes the window with nothing printed.
#
# What stops it is the RANGE pass being in `validate_star_params`, because that is the list
# the GUI's build gate reads (`build_epoch_star`: `isempty(validate_star_params(p)) || return
# nothing`). Moving that pass into `advise_star_params` — to keep the derived-mass message
# from disabling `epoch_chi2` — took the guard with it and cost 46 GUI callback tests. Only
# the mass check needed to move. These assertions are the difference between the two.
let good = default_star_params(2; rpole = RP, d = DPC, frac_escapevel = FEV,
                               rotation_period = PROT)
    checkbool("a sound model is still clean", isempty(validate_star_params(good)))
    for v in (0.0, -5.0)
        checkbool("rpole = $v BLOCKS the build",
                  !isempty(validate_star_params(merge(good, (rpole = v,)))))
    end
    # NaN fails every comparison, so the range test alone reports it as "outside the plausible
    # range" or skips it; it needs its own branch and its own message.
    bad = validate_star_params(merge(good, (rpole = NaN,)))
    checkbool("rpole = NaN blocks too", !isempty(bad))
    checkbool("and says so as NaN, not as out of range",
              any(m -> occursin("cannot describe a surface", m), bad))
    # The complement, which is the half that had to stay advisory.
    checkbool("an implausible MASS still does not block",
              isempty(validate_star_params(merge(good, (rotation_period = 50.0,)))))
end

println("\n[14] _map_range is total — no colour range can throw or come back degenerate")
# Every canvas funnels its colour range through this, from inside QML callbacks where a throw
# freezes the window rather than printing. The two inputs a plotting helper would not expect
# are an ALL-NaN map (above) and an EMPTY one (a total eclipse, or an orthographic view with
# nothing visible), on which `minimum` throws "reducing over an empty collection".
let mr = ROTIR._map_range
    for (what, v) in (("all NaN", fill(NaN, 7)), ("all NaN (Float32)", fill(NaN32, 7)),
                      ("empty", Float64[]), ("mixed finite/NaN", [1.0, NaN, 9.0, Inf]),
                      ("all Inf", [Inf, -Inf]), ("uniform", fill(5750.0, 4)),
                      ("one element", [3.0]), ("integers", [1, 2, 3]))
        lo, hi = mr(v)
        checkbool("$what -> finite and ordered", isfinite(lo) && isfinite(hi) && hi > lo)
        # And it must survive the Float32 narrowing the Colorbar applies, which a width of
        # 1e-30 at 5750 would not.
        checkbool("$what -> survives Float32", Float32(hi) > Float32(lo))
    end
    # The good values set the range; the bad ones are dropped, not mapped to zero.
    lo, hi = mr([1.0, NaN, 9.0, Inf])
    checkbool("non-finite values are dropped, not counted", lo == 1.0 && hi == 9.0)
    # A uniform map widens SYMMETRICALLY, or it renders solid black (see the docstring).
    lo, hi = mr(fill(5750.0, 4))
    checkbool("a uniform map widens about its value", lo < 5750.0 < hi)
    # An explicit range is honoured, and a non-finite or inverted one falls back rather than
    # being passed through to Makie.
    checkbool("explicit vmin/vmax honoured", mr([1.0, 2.0]; vmin = 0.0, vmax = 100.0) == (0.0, 100.0))
    lo, hi = mr([1.0, 9.0]; vmin = NaN)
    checkbool("a NaN vmin falls back to the data", isfinite(lo) && hi > lo)
    lo, hi = mr([1.0, 9.0]; vmin = 100.0, vmax = 0.0)
    checkbool("an inverted explicit range is swapped", hi > lo)
end

@printf("\n=== %d passed, %d failed ===\n", npass[], nfail[])

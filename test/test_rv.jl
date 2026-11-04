# Radial velocities (src/radial_velocity.jl) and the RV term of fit_orbit (src/orbit_fit.jl).
#
# Three kinds of check:
#   * the CSV format is OITOOLS' — both packages read the same files, and `binary_RV` equals
#     OITOOLS' `orbit_rv`, which pins the ω convention and the a/b ↔ primary/secondary mapping;
#   * synthetic round trips for the fit, the instrument offsets and the jitter;
#   * real data: HD 210763 (Gallenne et al. 2023) against PMOIRED's published fit, and Spica's
#     CHARA visibilities with its RVs.
using Test
using ROTIR
import OITOOLS
using NLopt
using Random

const RVDATA = joinpath(dirname(@__DIR__), "demos", "data")

# β Lyr uv coverage for the synthetic interferometry (the same files test_orbit_fit.jl uses)
const BLYR = joinpath(dirname(@__DIR__), "demos", "betlyr", "older_data")

@testset "radial velocities" begin
    @testset "CSV: the OITOOLS format" begin
        sp = read_rv(joinpath(RVDATA, "spica_rv.csv"))
        @test length(sp) == 66 && count(==(:a), sp.star) == 33
        @test sp.jd[1] ≈ 2451624.685 atol = 1e-8                 # MJD 51624.185 → JD
        @test sp.rv[1] == -75.4 && sp.sigma[1] == 1.0
        @test rv_instruments(sp) == [""]
        # the legacy two-file format gives the same table
        lg = read_rv(joinpath(RVDATA, "all_rv_1_ORIG.txt"), joinpath(RVDATA, "all_rv_2_ORIG.txt"))
        @test lg.jd ≈ sp.jd atol = 1e-8
        @test lg.rv == sp.rv && lg.star == sp.star
        hd = read_rv(joinpath(RVDATA, "hd210763_rv_gallenne2023.csv"))
        @test length(hd) == 34 && count(==(:a), hd.star) == 17

        # OITOOLS reads ROTIR's files to the same numbers, and the other way round
        for f in ("spica_rv.csv", "hd210763_rv_gallenne2023.csv")
            r = read_rv(joinpath(RVDATA, f))
            o = OITOOLS.read_rv(joinpath(RVDATA, f))
            @test r.jd .- 2400000.5 ≈ o.mjd atol = 1e-8
            @test r.rv == o.rv && r.sigma == o.sigma && r.star == o.star
        end

        mktempdir() do dir
            # write → read round trip, with and without instruments
            p = write_rv(joinpath(dir, "a.csv"), hd; comment = "two lines\nof provenance")
            @test startswith(readline(p), "# two lines")
            r = read_rv(p)
            @test r.jd ≈ hd.jd atol = 1e-7
            @test r.rv == hd.rv && r.sigma == hd.sigma && r.star == hd.star
            @test OITOOLS.read_rv(p).rv == hd.rv                  # OITOOLS reads what we write
            inst = [isodd(k) ? "UVES" : "SOPHIE" for k in 1:length(hd)]
            hi = RVData(hd.jd, hd.rv, hd.sigma, hd.star, inst, "")
            r2 = read_rv(write_rv(joinpath(dir, "b.csv"), hi))
            @test r2.instrument == inst && rv_instruments(r2) == ["UVES", "SOPHIE"]
            @test length(OITOOLS.read_rv(joinpath(dir, "b.csv"))) == length(hd)   # extra column ignored

            # a JD column instead of MJD, any column order, case-insensitive
            write(joinpath(dir, "c.csv"), "# c\nstar,Sigma,JD,rv\nA,0.5,2459510.5,3.0\nb,0.7,2459511.5,-2.0\n")
            c = read_rv(joinpath(dir, "c.csv"))
            @test c.jd == [2459510.5, 2459511.5] && c.star == [:a, :b] && c.sigma == [0.5, 0.7]

            # refusals
            bad(s) = (f = joinpath(dir, "bad.csv"); write(f, s); f)
            @test_throws ArgumentError read_rv(bad("MJD,RV,sigma\n1,2,3\n"))            # no star
            @test_throws ArgumentError read_rv(bad("RV,sigma,star\n2,3,a\n"))            # no time
            @test_throws ArgumentError read_rv(bad("MJD,JD,RV,sigma,star\n1,1,2,3,a\n")) # both
            @test_throws ArgumentError read_rv(bad("MJD,RV,sigma,star\n1,2,3,c\n"))      # star c
            @test_throws ArgumentError read_rv(bad("MJD,RV,sigma,star\n1,2,0,a\n"))      # σ = 0
            @test_throws ArgumentError read_rv(bad("1 2 3\n"))                           # not CSV
            @test_throws ArgumentError read_rv(bad("MJD,RV,sigma,star\n"))              # no rows
        end
    end

    @testset "binary_RV equals OITOOLS' orbit_rv" begin
        # OITOOLS: MJD, MJD0 = periastron, ω of the secondary relative to the primary,
        # VA = γ − Ka[cos(ν+ω) + e cos ω]. ROTIR: JD, T0 = periastron, ω of the relative orbit.
        # Same quantities, so the velocities must agree to rounding.
        for (e, ω) in ((0.0, 90.0), (0.123, 255.0), (0.62, 114.0), (0.9, 10.0))
            bp = (i = 70.0, Ω = 30.0, ω = ω, P = 42.38, a = 3.0, e = e, T0 = 2459785.98,
                  q = 1.0, dP = 0.0)
            t = 2459785.98 .+ collect(range(-50.0, 120.0, length = 37))
            v1, v2 = binary_RV(bp, t; K1 = 50.3, K2 = 58.8, γ = 14.9)
            oa = [OITOOLS.orbit_rv(x - 2400000.5, :a, bp.P, bp.T0 - 2400000.5, e, ω, 50.3, 14.9) for x in t]
            ob = [OITOOLS.orbit_rv(x - 2400000.5, :b, bp.P, bp.T0 - 2400000.5, e, ω, 58.8, 14.9) for x in t]
            @test maximum(abs.(v1 .- oa)) < 1e-8
            @test maximum(abs.(v2 .- ob)) < 1e-8
        end
        # the widened signature: ranges, Ints, scalars, and the same numbers as Float64
        bp = (ω = 30.0, e = 0.2, P = 5.0, T0 = 0.0, dP = 0.0)
        @test binary_RV(bp, 0:4; K1 = 10, K2 = 20, γ = 1)[1] ≈
              binary_RV(bp, collect(0.0:4.0); K1 = 10.0, K2 = 20.0, γ = 1.0)[1]
        @test binary_RV(bp, 1.5; K1 = 10.0, K2 = 20.0, γ = 0.0)[2] isa Float64
    end

    @testset "model, residuals, offsets, jitter" begin
        bp = (i = 95.0, Ω = 70.0, ω = 40.0, P = 12.94, a = 0.9, e = 0.25, T0 = 2454283.0,
              q = 1.0, dP = 0.0, dω = 0.0)
        t = 2454283.0 .+ collect(range(0.3, 40.0, length = 30))
        star = [isodd(k) ? :a : :b for k in 1:30]
        inst = [k <= 15 ? "UVES" : "SOPHIE" for k in 1:30]
        m = rv_model(bp, RVData(t, zeros(30), ones(30), star, inst, "");
                     K1 = 60.0, K2 = 90.0, gamma = -12.0, offsets = Dict("SOPHIE" => 0.8))
        rv = RVData(t, m, fill(0.5, 30), star, inst, "")
        kw = (K1 = 60.0, K2 = 90.0, gamma = -12.0, offsets = Dict("SOPHIE" => 0.8))
        @test rv_chi2(bp, rv; kw...) < 1e-20
        # residual = (rv − model)/σ: shifting γ by 1 gives −1/σ everywhere
        @test rv_residuals(bp, rv; kw..., gamma = -11.0) ≈ fill(-2.0, 30)
        # an offset moves only its own instrument's rows
        r = rv_residuals(bp, rv; kw..., offsets = Dict("SOPHIE" => 0.3))
        @test all(iszero, r[1:15]) && r[16:end] ≈ fill(1.0, 15)
        # jitter inflates σ in quadrature, per star
        r = rv_residuals(bp, rv; kw..., gamma = -11.0, jitter = (a = 1.2, b = 0.0))
        @test r[1] ≈ -1 / sqrt(0.25 + 1.44) && r[2] ≈ -2.0
        # rv_nll2 is χ² at zero jitter, and adds Σ ln((σ²+s²)/σ²) otherwise
        @test rv_nll2(bp, rv; kw..., gamma = -11.0) ≈ rv_chi2(bp, rv; kw..., gamma = -11.0)
        j = (a = 1.0, b = 2.0)
        extra = sum(log((0.25 + (s === :a ? 1.0 : 4.0)) / 0.25) for s in star)
        @test rv_nll2(bp, rv; kw..., gamma = -11.0, jitter = j) ≈
              rv_chi2(bp, rv; kw..., gamma = -11.0, jitter = j) + extra
        st = rv_start(rv)
        @test st.gamma ≈ sum(rv.rv) / 30 && st.K1 > 0 && st.K2 > st.K1
        @test ROTIR.rv_offset_name("SOPHIE-2") == :dgamma_SOPHIE_2
        @test ROTIR.rv_offset_name("2m") == :dgamma_i2m
    end

    @testset "fit_orbit with velocities (synthetic)" begin
        files = sort(filter(f -> endswith(f, ".oifits"), readdir(BLYR, join = true)))
        if isempty(files)
            @info "skipping the synthetic RV fit: no OIFITS under $BLYR"
        else
            el = (a = 0.9, i = 95.0, Omega = 250.0, omega = 40.0, e = 0.25, P = 12.9414,
                  T0 = 2454283.0, dP = 0.0, domega = 0.0)
            bp = (i = el.i, Ω = el.Omega, ω = el.omega, P = el.P, a = el.a, e = el.e,
                  T0 = el.T0, q = 1.0, dP = 0.0, dω = 0.0)
            K1, K2, γ = 60.0, 90.0, -12.0
            t = el.T0 .+ collect(range(0.3, 40.0, length = 40))
            star = [isodd(k) ? :a : :b for k in 1:40]
            inst = [k <= 20 ? "UVES" : "SOPHIE" for k in 1:40]
            v = rv_model(bp, RVData(t, zeros(40), ones(40), star, inst, "");
                         K1 = K1, K2 = K2, gamma = γ, offsets = Dict("SOPHIE" => 0.8))
            rv = RVData(t, v, fill(0.5, 40), star, inst, "")
            c1, c2 = UniformDisk(diameter = 0.55), GaussianDisk(fwhm = 0.80)
            ds = [readoifits(f; T = Float64)[1, 1] for f in files]
            # The amplitudes follow from the orbit: plx and q are the parameters.
            plx0, q0 = rv_physical(el.a, el.i, el.P, el.e, K1, K2)
            @test all(rv_amplitudes(el.a, el.i, el.P, el.e, plx0, q0) .≈ (K1, K2))
            @test q0 ≈ K1 / K2
            spec = orbit_fit_spec(c1, c2; elements = el, flux_ratio = 0.8, rv = rv,
                                  rv_params = (plx = plx0, q = q0, gamma = γ, dgamma_SOPHIE = 0.8))
            @test spec.names[end-5:end] == [:plx, :q, :gamma, :dgamma_SOPHIE, :jit_a, :jit_b]
            @test spec.hi[3] == 360.0                         # Omega: the mirror is broken
            @test :plx in spec.names[spec.free] && :dgamma_SOPHIE in spec.names[spec.free]
            # the K1/K2 parametrisation gives the same model
            speck = orbit_fit_spec(c1, c2; elements = el, flux_ratio = 0.8, rv = rv,
                                   rv_parametrization = :amplitudes,
                                   rv_params = (K1 = K1, K2 = K2, gamma = γ, dgamma_SOPHIE = 0.8))
            @test speck.names[end-5:end-4] == [:K1, :K2]
            @test !(:jit_a in spec.names[spec.free])          # jitters only when asked
            fd = orbit_fit_data(ds)
            for (k, d) in enumerate(ds)                       # noiseless OI from the same orbit
                v2m, t3am, t3pm = cvis_to_obs(orbit_model_cvis(spec, spec.values, fd, k), d)
                d.v2 .= v2m; d.t3amp .= t3am; d.t3phi .= t3pm
            end
            fd = orbit_fit_data(ds)
            @test orbit_chi2(spec, spec.values, fd; rv = rv) < 1e-8
            @test orbit_chi2(speck, speck.values, fd; rv = rv) < 1e-8
            @test length(orbit_chi2(spec, spec.values, fd; rv = rv, split = true)) == 4

            # the residual vector is the χ² parts
            θ = copy(spec.values); θ[1] += 0.01; θ[end-3] += 2.0
            r = orbit_residuals(spec, θ, fd; rv = rv)
            @test length(r) == fd.ntot + length(rv)
            @test sum(abs2, r) ≈ sum(orbit_chi2(spec, θ, fd; rv = rv, split = true)) rtol = 1e-10

            # the mirror: identical interferometry, very different velocities
            θm = copy(spec.values); θm[3] -= 180; θm[4] += 180
            cm = orbit_chi2(spec, θm, fd; rv = rv, split = true)
            @test sum(cm[1:3]) < 1e-8 && cm[4] > 1e4

            # refusals: velocities without an RV block, mismatched instruments, rv_params alone
            # (without velocities Omega is boxed on [0, 180), so take the mirrored branch)
            el0 = merge(el, (Omega = el.Omega - 180, omega = el.omega + 180))
            sp0 = orbit_fit_spec(c1, c2; elements = el0, flux_ratio = 0.8)
            @test_throws ErrorException orbit_chi2(sp0, sp0.values, fd; rv = rv)
            rv1 = RVData(t, v, fill(0.5, 40), star, fill("", 40), "")
            @test_throws ErrorException orbit_chi2(spec, spec.values, fd; rv = rv1)
            @test_throws ErrorException orbit_fit_spec(c1, c2; elements = el0, rv_params = (plx = 1.0,))
            @test_throws ErrorException orbit_fit_spec(c1, c2; elements = el, rv = rv,
                                                       rv_params = (K3 = 1.0,))

            # recover ω, e, T0, plx, q, γ and the instrument offset from a perturbed start
            start = merge(el, (omega = 50.0, e = 0.2, T0 = el.T0 + 0.3))
            res = fit_orbit(ds, c1, c2; elements = start, flux_ratio = 0.8, rv = rv,
                            rv_params = (plx = 1.1plx0, q = 0.9q0, gamma = -10.0),
                            free = [:omega, :e, :T0, :plx, :q, :gamma, :dgamma_SOPHIE],
                            maxeval = 40_000, verbose = false)
            @test res.chi2_red < 1e-5                       # NM stopping tolerance; the parameters below are the test
            @test res.elements.omega ≈ el.omega atol = 0.05
            @test res.elements.e ≈ el.e atol = 1e-3
            @test res.rvpar.plx ≈ plx0 rtol = 1e-3
            @test res.rvpar.q ≈ q0 rtol = 1e-3
            @test res.rvpar.K1 ≈ K1 atol = 0.05                # derived, reported
            @test res.rvpar.K2 ≈ K2 atol = 0.05
            dm = dynamical_masses(el.a, el.P, res.rvpar.plx, res.rvpar.q)
            @test res.rvpar.M ≈ dm.M && res.rvpar.M2 / res.rvpar.M1 ≈ res.rvpar.q
            @test res.rvpar.d ≈ 1000 / res.rvpar.plx
            @test res.rvpar.gamma ≈ γ atol = 0.05
            @test res.rvpar.dgamma_SOPHIE ≈ 0.8 atol = 0.05
            @test res.rvpar.jit_a == 0.0                      # fixed, untouched

            # jitter: add 1.5 km/s of scatter to the secondary and fit it (RVs only)
            rng = MersenneTwister(20261023)
            noisy = copy(v)
            for k in eachindex(noisy); star[k] === :b && (noisy[k] += 1.5randn(rng)); end
            rvn = RVData(t, noisy, fill(0.5, 40), star, inst, "")
            rj = fit_orbit(ds, c1, c2; elements = el, flux_ratio = 0.8, rv = rvn,
                           rv_parametrization = :amplitudes,
                           rv_params = (K1 = K1, K2 = K2, gamma = γ, dgamma_SOPHIE = 0.8,
                                        jit_b = 0.5),
                           free = [:K1, :K2, :gamma, :jit_a, :jit_b], weights = [0.0, 0.0, 0.0],
                           maxeval = 20_000, verbose = false)
            @test 0.8 < rj.rvpar.jit_b < 2.5          # 20 points: s is known to ~±0.3
            @test rj.rvpar.jit_a < 0.2                # no extra scatter on the primary
        end
    end

    @testset "HD 210763: RVs against the published fit" begin
        # Gallenne et al. 2023 table C.1 (UVES, SOPHIE), 34 velocities. Elements start from
        # their table 6 with ω converted to the relative orbit (293.965° − 180°), as OITOOLS'
        # EX6 walk does; P is held at 42.3811 d. PMOIRED's joint GRAVITY + RV fit gives
        # Ka 50.309 ± 0.054, Kb 58.76 ± 0.12, γ 14.929 ± 0.028, e 0.62295 ± 0.00052.
        # These are RVs alone, so agreement to ~2σ is what to expect.
        hd = read_rv(joinpath(RVDATA, "hd210763_rv_gallenne2023.csv"))
        P = 42.3811; T0 = 54275.93 + 42.3811 * 130 + 2400000.5
        f(p) = rv_chi2((ω = p[1], e = p[2], P = P, T0 = p[3], dP = 0.0), hd;
                       K1 = p[4], K2 = p[5], gamma = p[6])
        opt = NLopt.Opt(:LN_NELDERMEAD, 6)
        NLopt.min_objective!(opt, (p, g) -> f(p))
        NLopt.lower_bounds!(opt, [0.0, 0.0, T0 - 5, 0.0, 0.0, -100.0])
        NLopt.upper_bounds!(opt, [360.0, 0.95, T0 + 5, 200.0, 200.0, 100.0])
        NLopt.maxeval!(opt, 40_000); NLopt.xtol_rel!(opt, 1e-12)
        (χ2, p, _) = NLopt.optimize(opt, [113.965, 0.6228, T0, 50.28, 58.76, 14.94])
        @test abs(p[4] - 50.309) < 3 * 0.054
        @test abs(p[5] - 58.76)  < 3 * 0.12
        @test abs(p[6] - 14.929) < 3 * 0.028
        @test abs(p[2] - 0.62295) < 3 * 0.00052
        @test χ2 / (length(hd) - 6) < 1.5             # the tabulated errors fit the data
    end

    @testset "Spica: CHARA visibilities and RVs" begin
        raw = readoifits(joinpath(RVDATA, "2007_2012_2015.Spica.oifits"); T = Float64)[1, 1]
        mj = sort(raw.v2_mjd); jm = findall(diff(mj) .> 0.5)
        st, sp = mj[[1; jm .+ 1]], mj[[jm; length(mj)]]
        data = [d for d in [filter_data(raw, set_data_filter(raw; mjd_range = [st[i] - 0.01, sp[i] + 0.01]))
                            for i in eachindex(st)] if d.nv2 > 0]
        rv0 = read_rv(joinpath(RVDATA, "spica_rv.csv"))
        # Aufdenberg et al. 2015 elements, apsidal period 139 yr (Robinette & Aufdenberg 2015)
        lit = (a = 1.54, i = 116.0, Omega = 309.938, omega = 255.0, e = 0.123, P = 4.0145,
               T0 = 2454189.40, dP = 0.0, domega = 360 / (139 * 365.25))
        c1 = LimbDarkenedDisk(diameter = 0.894, law = :power, ld1 = 0.15)
        c2 = LimbDarkenedDisk(diameter = 0.452, law = :power, ld1 = 0.15)
        # plx, q giving the literature amplitudes 123.9 / 198.8 km/s at these elements
        Kl = (K1 = 123.9, K2 = 198.8)
        p0, q0 = rv_physical(lit.a, lit.i, lit.P, lit.e, Kl.K1, Kl.K2)
        K = (plx = p0, q = q0, gamma = 0.0)
        @test 10 < p0 < 16                                   # Spica: ϖ ≈ 12.4-13.1 mas
        @test q0 ≈ 123.9 / 198.8
        fd = orbit_fit_data(data)

        # The velocities choose between the two orbits interferometry cannot tell apart.
        s = orbit_fit_spec(c1, c2; elements = lit, flux_ratio = 0.196, rv = rv0, rv_params = K)
        θm = copy(s.values); θm[3] -= 180; θm[4] -= 180
        c, cm = orbit_chi2(s, s.values, fd; rv = rv0, split = true),
                orbit_chi2(s, θm, fd; rv = rv0, split = true)
        @test collect(c[1:3]) ≈ collect(cm[1:3]) rtol = 1e-10
        @test cm[4] > 10 * c[4]

        # The files' σ = 1 km/s is a placeholder; the blended-line error model of
        # demos/spica_orbit_fit_ultranest.jl, σ = 7 + 30 exp(−|V2 − V1|/100) km/s, evaluated
        # once at the literature orbit, makes it usable. RVs alone then give e = 0.11 —
        # Robinette & Aufdenberg 2015 find 0.116 ± 0.004 from RVs spanning 1889–2000. The
        # same minimum is reached from several starts; ω is boxed on [180°, 540°) because
        # this solution sits near ω = 0°, on the edge of the default [0°, 360°] box.
        bp = (ω = lit.omega, e = lit.e, P = lit.P, T0 = lit.T0, dP = 0.0, dω = lit.domega)
        v1, v2 = binary_RV(bp, rv0.jd; K1 = Kl.K1, K2 = Kl.K2, γ = 0.0)
        rv = RVData(rv0.jd, rv0.rv, 7.0 .+ 30.0 .* exp.(-abs.(v2 .- v1) ./ 100.0), rv0.star,
                    rv0.instrument, "")
        es = Float64[]
        for (ΔT0, ω0) in ((1.0, 255.0), (0.5, 300.0))
            # RVs alone constrain the AMPLITUDES; plx and q need the interferometry's a and i,
            # which this fit holds fixed — so the amplitudes are the parameters here.
            r = fit_orbit(data, c1, c2; elements = merge(lit, (T0 = lit.T0 + ΔT0, omega = ω0)),
                          flux_ratio = 0.196, rv = rv, rv_parametrization = :amplitudes,
                          rv_params = (K1 = Kl.K1, K2 = Kl.K2, gamma = 0.0),
                          free = [:omega, :e, :T0, :K1, :K2, :gamma],
                          bounds = Dict(:omega => (180.0, 540.0)), weights = [0.0, 0.0, 0.0],
                          maxeval = 20_000, verbose = false)
            push!(es, r.elements.e)
            @test 0.09 < r.elements.e < 0.13
            @test 115 < r.rvpar.K1 < 130 && 180 < r.rvpar.K2 < 200
            @test r.chi2_split[4] / length(rv) < 2
            @test sum(r.chi2_split[1:3]) == 0                 # OI weights 0: not evaluated
        end
        @test abs(es[1] - es[2]) < 1e-3

        # A joint fit from the literature lowers the total and stays on the RV branch of Ω.
        r = fit_orbit(data, c1, c2; elements = lit, flux_ratio = 0.196, rv = rv, rv_params = K,
                      free = [:a, :i, :Omega, :omega, :e, :T0, :f, :plx, :q, :gamma],
                      maxeval = 1500, verbose = false)
        s2 = orbit_fit_spec(c1, c2; elements = lit, flux_ratio = 0.196, rv = rv, rv_params = K)
        @test r.chi2 < sum(orbit_chi2(s2, s2.values, fd; rv = rv, split = true))
        @test 180 < r.elements.Omega < 360
        @test r.chi2_red ≈ r.chi2 / (fd.ntot + length(rv))
    end
end

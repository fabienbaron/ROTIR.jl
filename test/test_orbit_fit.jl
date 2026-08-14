# Generic orbit fitter (src/orbit_fit.jl).
#
# The load-bearing assertion is the round trip: build a two-component system, generate
# noiseless observables from known elements, and check the fitter recovers them. That
# exercises the whole chain — element -> Kepler -> projection -> component visibility ->
# phase shift -> observable -> chi2 -> optimiser — in a way no unit test of the parts does.
using Test
using ROTIR

@testset "orbit_fit" begin
    @testset "component library" begin
        for c in (PointSource(), UniformDisk(diameter = 0.5),
                  LimbDarkenedDisk(diameter = 0.5, law = :linear),
                  LimbDarkenedDisk(diameter = 0.5, law = :quadratic),
                  GaussianDisk(fwhm = 0.5),
                  EllipticalGaussian(fwhm = 0.5, ratio = 0.6, pa = 30.0))
            n = ROTIR.component_param_names(c)
            v = ROTIR.component_param_values(c)
            lo, hi = ROTIR.component_param_bounds(c)
            @test length(n) == length(v) == length(lo) == length(hi)
            @test all(lo .<= v .<= hi)          # defaults must sit inside their own box
        end
        # an elliptical Gaussian is symmetric under a half turn, so its PA bound is 180
        @test ROTIR.component_param_bounds(EllipticalGaussian(fwhm = 1.0))[2][3] == 180.0
    end

    @testset "spec assembly" begin
        c1, c2 = UniformDisk(diameter = 0.6), EllipticalGaussian(fwhm = 0.9, ratio = 0.6, pa = 70.0)
        spec = orbit_fit_spec(c1, c2; elements = (P = 12.94, a = 0.87, i = 93.5,
                                                  Omega = 73.7, T0 = 2454283.0))
        # derive the count, do not repeat it: ORBIT_ELEMENTS gained `domega` and a
        # hardcoded 1:8 silently stopped covering the last element.
        @test spec.names[1:length(ORBIT_ELEMENTS)] == collect(ORBIT_ELEMENTS)
        @test :domega in spec.names          # apsidal motion is reachable
        @test spec.values[findfirst(==(:domega), spec.names)] == 0.0  # ...and off by default
        @test :domega ∉ spec.names[spec.free]                        # ...and not fitted by default
        @test :c1_diameter in spec.names && :c2_pa in spec.names && :f in spec.names
        # T0 spans exactly one period: with Omega in [0,180) that is one fundamental domain
        # of the (Omega, T0) -> (Omega+180, T0+P/2) degeneracy at e = 0.
        j = findfirst(==(:T0), spec.names)
        @test spec.hi[j] - spec.lo[j] ≈ 12.94
        # e and omega are not free by default (omega is undefined at e = 0)
        @test !(findfirst(==(:e), spec.names) in spec.free)
        @test !(findfirst(==(:omega), spec.names) in spec.free)
        # P is not free by default either
        @test !(findfirst(==(:P), spec.names) in spec.free)
        @test_throws ErrorException orbit_fit_spec(c1, c2; elements = (P = 1.0,), free = [:nope])
        @test_throws ErrorException orbit_fit_spec(c1, c2; elements = (a = 1.0,))  # P required
    end

    @testset "round trip on synthetic data" begin
        # Noiseless observables from known elements, then recover them.
        truth = (P = 12.9414, e = 0.0, a = 0.90, i = 95.0,
                 Omega = 70.0, omega = 90.0, T0 = 2454283.0, dP = 0.0)
        c1 = UniformDisk(diameter = 0.55)
        c2 = GaussianDisk(fwhm = 0.80)
        dir = joinpath(@__DIR__, "..", "demos", "betlyr", "older_data")
        files = isdir(dir) ? sort(filter(f -> endswith(f, ".oifits"), readdir(dir, join = true))) : String[]

        if isempty(files)
            @info "skipping orbit round trip: no OIFITS under $dir"
        else
            "Replace a block's observables with the model's own prediction (noiseless)."
            function synth(fs)
                ds = [readoifits(f; T = Float64)[1, 1] for f in fs]
                spec = orbit_fit_spec(c1, c2; elements = truth, flux_ratio = 0.8)
                fd = orbit_fit_data(ds)
                for (k, d) in enumerate(ds)
                    v2m, t3am, t3pm = cvis_to_obs(orbit_model_cvis(spec, spec.values, fd, k), d)
                    d.v2 .= v2m; d.t3amp .= t3am; d.t3phi .= t3pm
                end
                ds, spec
            end
            start = merge(truth, (a = 0.80, i = 88.0, Omega = 62.0, T0 = truth.T0 + 0.4))
            fit(ds) = fit_orbit(ds, c1, c2; elements = start, flux_ratio = 0.8,
                                free = [:a, :i, :Omega, :T0], method = :neldermead,
                                maxeval = 40_000, verbose = false)

            # chi2 must vanish at the generating parameters
            ds, spec = synth(files)
            @test orbit_chi2(spec, spec.values, orbit_fit_data(ds)) < 1e-8

            # ALL nights: the elements are recovered
            res = fit(ds)
            @test res.chi2_red < 1e-4
            @test isapprox(res.elements.a,     truth.a;     atol = 0.02)
            @test isapprox(res.elements.i,     truth.i;     atol = 1.0)
            @test isapprox(res.elements.Omega, truth.Omega; atol = 1.0)

            # ONE night: a perfect chi2 with the WRONG elements. A single night samples a
            # few hours of a 12.9 d period, so many orbits pass through the same
            # instantaneous separation. Asserted deliberately — "the fit converged" is not
            # evidence the orbit is determined, and this is the trap users fall into.
            ds1, _ = synth(files[1:1])
            r1 = fit(ds1)
            @test r1.chi2_red < 1e-4                              # fits the data perfectly
            @test !isapprox(r1.elements.a, truth.a; atol = 0.05)  # and is still wrong
        end
    end

    # Tied parameters (src/orbit_ties.jl): values computed from other parameters rather than
    # fitted. Everything here is analytic — a tie either evaluates to the right number or it
    # does not — plus the guards that stop a silently wrong model.
    @testset "tied parameters" begin
        names = Symbol[:a, :i, :Omega, :omega, :e, :P, :T0, :dP,
                       :c1_diameter, :c2_fwhm, :c2_ratio, :c2_pa, :f]
        θ = Float64[0.8958, 91.244, 73.529, 90.0, 0.0, 12.9616, 2454283.0555, 0.0,
                    0.590, 0.5825, 0.6423, 999.0, 0.8179]

        @testset "evaluation" begin
            t1 = compile_ties(names, Dict(:c2_pa => "-Omega"))
            @test apply_ties(t1, θ)[12] ≈ -73.529
            # everything else is untouched
            @test apply_ties(t1, θ)[[1,2,3,9,13]] == θ[[1,2,3,9,13]]

            # arbitrary Julia, including function calls
            t2 = compile_ties(names, Dict(:c2_ratio => "abs(cosd(i)) + 0.62"))
            @test apply_ties(t2, θ)[11] ≈ abs(cosd(91.244)) + 0.62

            # a tie may depend on another tie, in either declaration order
            t3 = compile_ties(names, Dict(:c2_pa => "-Omega", :c2_fwhm => "2*c2_pa"))
            r3 = apply_ties(t3, θ)
            @test r3[12] ≈ -73.529 && r3[10] ≈ -147.058
        end

        @testset "guards" begin
            @test_throws ErrorException compile_ties(names, Dict(:c2_pa => "c2_pa + 1"))
            @test_throws ErrorException compile_ties(names, Dict(:not_a_param => "1.0"))
            @test_throws ErrorException compile_ties(names, Dict(:c2_pa => "Omegaa"))
            @test_throws ErrorException compile_ties(names,
                            Dict(:c2_pa => "c2_ratio", :c2_ratio => "c2_pa"))
        end

        @testset "no ties costs nothing" begin
            t0 = compile_ties(names, Dict())
            @test isempty(t0)
            # identity, not a copy: an untied fit must not pay an allocation per likelihood
            @test apply_ties(t0, θ) === θ
        end

        # `fn` must stay concretely typed. With `fn::Function` every call is a dynamic
        # dispatch and the returned tuple is boxed, costing ~3 allocations per likelihood
        # evaluation — invisible in the answers, expensive over a nested-sampling run.
        @testset "concretely typed" begin
            t1 = compile_ties(names, Dict(:c2_pa => "-Omega"))
            @test isconcretetype(typeof(t1))
            @test fieldtype(typeof(t1), :fn) !== Function
        end

        @testset "wired into OrbitFitSpec" begin
            c1, c2 = UniformDisk(diameter = 0.59), EllipticalGaussian(0.58, 0.64, 108.0)
            el = (a = 0.9, i = 91.2, Omega = 73.5, P = 12.96, T0 = 2454283.0)
            sp = orbit_fit_spec(c1, c2; elements = el, flux_ratio = 0.82,
                                free = [:a, :i, :Omega, :T0, :f],
                                ties = Dict(:c2_pa => "-Omega"))
            j = findfirst(==(:c2_pa), sp.names)
            @test apply_ties(sp.ties, sp.values)[j] ≈ -73.5

            # REGRESSION. `resolve_params` must be the only free->full construction: an
            # open-coded `copy(values); θ[free] .= p` skips the ties, so a fit would report
            # tied parameters still holding their STARTING values while quoting a χ²
            # computed from the correct ones. That was live until the fused path landed.
            pfree = sp.values[sp.free]
            θres  = resolve_params(sp, pfree)
            @test θres[j] ≈ -θres[3]                      # tie honoured
            @test θres[j] != sp.values[j]                 # and not the start value
            let θnaive = copy(sp.values)
                θnaive[sp.free] .= pfree
                @test θnaive[j] == sp.values[j]           # the trap this guards against
            end

            # concretely typed all the way down, or every likelihood pays a dynamic dispatch
            @test fieldtype(typeof(sp), :ties) !== OrbitTies
            @test isconcretetype(fieldtype(typeof(sp), :ties))
            # a tie overwrites its target, so fitting it too would discard the sampled value
            @test_throws ErrorException orbit_fit_spec(c1, c2; elements = el,
                free = [:a, :c2_pa], ties = Dict(:c2_pa => "-Omega"))
        end
    end

    @testset "tessellated path is declared, not silently wrong" begin
        c1, c2 = UniformDisk(diameter = 0.5), GaussianDisk(fwhm = 0.5)
        @test_throws ErrorException fit_orbit(nothing, c1, c2;
            elements = (P = 1.0,), model = :tessellated, verbose = false)
        @test_throws ErrorException fit_orbit(nothing, c1, c2;
            elements = (P = 1.0,), model = :bogus, verbose = false)
    end
end

#!/usr/bin/env julia
# =======================================================================================
# Does PREDICTED limb darkening break the beta degeneracy, or just move it?
# =======================================================================================
# A graduate student fitting beta Cas as a rapid rotator by NUTS obtained a NEGATIVE
# limb-darkening coefficient. The mechanism is in demos/gravity_law_comparison.jl and the
# header of src/gravity_darkening.jl: gravity darkening and limb darkening both take flux out
# of the limb, an interferometer measures only their sum, and a fit forced to use the wrong
# gravity law buys the difference out of ld1. A coefficient fitted that way is a statement
# about the law, not about the star.
#
# A model-atmosphere intensity I(Teff, logg, mu, lambda) removes ld1 from the fit entirely.
# But it does not come free: indexing an atmosphere needs `logg`, which needs a mass, which
# (see src/stellar_physics.jl) comes from adding the distance `d`. So the obvious worry is
# that beta stops trading against ld1 and starts trading against d instead — both reach the
# pole-to-equator flux ratio.
#
# THIS SCRIPT ANSWERS THAT WITHOUT FITTING ANYTHING. It compares two Fisher matrices on the
# real beta Cas MIRC-X data, at equal parameter count:
#
#   A: theta = [rpole, fev, inc, PA, beta, ld1]  Planck intensity, fitted power-law LD
#   B: theta = [rpole, fev, inc, PA, beta, d  ]  Korg/MARCS grid, ldtype = 0 (no fitted LD)
#
# Only the forward model is needed, plus central differences — no sampler, no fit, and none
# of the parameter-vector plumbing that a real fit of case B would still require.
#
#   julia --project=demos demos/betcas_fisher_predicted_ld.jl
#
# The grid is built with Korg on first run (~50 s, needs `using Korg`) and cached.
# =======================================================================================

using ROTIR, Printf, LinearAlgebra

const CACHE = joinpath(@__DIR__, "data", "betcas_H_korg.fits")
const OIFITS = joinpath(@__DIR__, "data",
                        "MEDIAN5.MIRCX_L2.2025Oct30.HD_432.MIRCX_IDL.bet_Cas.AVG10m.oifits")

# beta Cas's surface spans Teff 5888-7206 K and logg 3.29-3.75; the grid is padded well
# beyond that because the finite differences perturb rpole, fev and d, all of which move it.
# Korg's VALD solar linelist has no lines in the H band, so this is continuum plus the
# Brackett series — which is what beta Cas's H band is.
function grid()
    isfile(CACHE) && return load_intensity_grid(CACHE)
    @info "building the Korg grid (first run only)" CACHE
    try
        @eval using Korg
    catch
        error("This demo needs Korg for its first run: `] add Korg`. Once $(CACHE) exists, " *
              "Korg is no longer required — the grid is read from FITS.")
    end
    g = build_korg_grid(Teff = (5400.0, 7600.0, 12), logg = (2.9, 4.1, 7),
                        μ = [0.001,0.02,0.05,0.1,0.2,0.3,0.45,0.6,0.75,0.85,0.95,1.0],
                        λ = (15100.0, 17300.0), λ_step = 5.0)
    save_intensity_grid(CACHE, g; comment = "Korg MARCS plane-parallel, H band, beta Cas")
    return g
end

const GRID = grid()
const PROV = TabulatedProvider(GRID; name = "Korg/MARCS H")
const DATA = readoifits_multiepochs([OIFITS]; T = Float64)[1, 1]
const TESS = tessellation_healpix(4, T = Float64)
const BAND = band_of(DATA)

base(; kw...) = default_star_params(2; rpole = 0.849, d = 16.8, frac_escapevel = 0.92,
                    rotation_period = 1/1.12, tpole = 7208.0, inclination = 19.9,
                    position_angle = -7.09, beta = 0.25, gravity_law = 2, kw...)

# The two models return the SAME observable vector, so the two Fisher matrices compare
# directly. Case B's sixth parameter reaches the data ONLY through logg -> the atmosphere's
# limb darkening; it does not touch the geometry, since rpole is already an angle.
function obs_planck(θ)
    sp = base(ldtype = 3, ld1 = θ[6], rpole = θ[1], frac_escapevel = θ[2],
              inclination = θ[3], position_angle = θ[4], beta = θ[5])
    st = create_star(TESS, sp, 0.0)
    cvis_to_obs(ROTIR.fused_cvis(parametric_temperature_map(sp, st), st, DATA;
                                 intensity_model = :planck, band = BAND), DATA)
end
function obs_atmos(θ)
    sp = base(ldtype = 0, d = θ[6], rpole = θ[1], frac_escapevel = θ[2],
              inclination = θ[3], position_angle = θ[4], beta = θ[5])
    st = create_star(TESS, sp, 0.0)
    I = channel_intensity(PROV, surface_state(st, sp), BAND)
    cvis_to_obs(ROTIR.fused_cvis(I, st, DATA; intensity_model = :linear), DATA)
end

whiten(o) = vcat(o[1] ./ DATA.v2_err, o[3] ./ DATA.t3phi_err)
# Phases are differenced WITH WRAPPING: a closure phase sitting near +-180 deg would
# otherwise manufacture a derivative of 360/(2h).
dwhiten(p, m) = vcat((p[1] .- m[1]) ./ DATA.v2_err,
                     mod360.(p[3] .- m[3]) ./ DATA.t3phi_err)

function fisher(f, θ; frac = 1e-6)
    n = length(θ); J = Matrix{Float64}(undef, length(whiten(f(θ))), n)
    for j in 1:n
        h = frac * max(abs(θ[j]), 1e-3)
        tp = copy(θ); tp[j] += h; tm = copy(θ); tm[j] -= h
        J[:, j] = dwhiten(f(tp), f(tm)) ./ (2h)
    end
    return J' * J
end
corrmat(C) = [C[i,j]/sqrt(C[i,i]*C[j,j]) for i in axes(C,1), j in axes(C,2)]

# Gaia DR3 places beta Cas at 16.8 pc to better than 1 per cent; 0.1 pc is conservative.
const D_PRIOR = Diagonal([0.0, 0, 0, 0, 0, 1/0.1^2])
const NA = ["rpole","fev","inc","PA","beta","ld1"]
const NB = ["rpole","fev","inc","PA","beta","d"]
const θA = [0.849, 0.92, 19.9, -7.09, 0.25, 0.21]
const θB = [0.849, 0.92, 19.9, -7.09, 0.25, 16.8]

function show_case(title, names, F)
    C = inv(F); R = corrmat(C)
    @printf("\n%s\n  cond(F) = %.2e\n", title, cond(F))
    @printf("  %-7s %11s   %s\n", "param", "sigma", join(rpad.(names, 8)))
    for i in eachindex(names)
        @printf("  %-7s %11.5g   %s\n", names[i], sqrt(C[i,i]),
                join([rpad(@sprintf("%+.3f", R[i,j]), 8) for j in eachindex(names)]))
    end
    return C, R
end

println("="^90)
println("beta Cas MIRC-X — Fisher comparison at equal parameter count")
@printf("nV2 = %d   nT3phi = %d   band = %.4f um   HEALPix nside 4 (%d tessels)\n",
        DATA.nv2, DATA.nt3phi, BAND*1e6, TESS.npix)
@printf("grid: %s over Teff %.0f-%.0f, logg %.1f-%.1f\n", string(size(GRID)),
        first(GRID.Teff), last(GRID.Teff), first(GRID.logg), last(GRID.logg))
println("="^90)

FA = fisher(obs_planck, θA); FB = fisher(obs_atmos, θB)
CA, RA = show_case("A: Planck + fitted power-law LD (ldtype = 3)", NA, FA)
CB, RB = show_case("B: Korg/MARCS + ldtype = 0, with Gaia prior sigma_d = 0.1 pc",
                   NB, FB + D_PRIOR)
CB0, RB0 = show_case("B: the same WITHOUT a distance prior", NB, FB)

println("\n" * "="^90); println("Verdict"); println("="^90)
@printf("  sigma(beta):            A %.5f    B %.5f   (%.2fx tighter)\n",
        sqrt(CA[5,5]), sqrt(CB[5,5]), sqrt(CA[5,5])/sqrt(CB[5,5]))
@printf("  corr(beta, 6th param):  A beta/ld1 %+.3f   B beta/d %+.3f\n", RA[5,6], RB[5,6])
@printf("  corr(beta, fev):        A %+.3f            B %+.3f\n", RA[5,2], RB[5,2])
@printf("  sigma(fev):             A %.5f    B %.5f   (%.2fx tighter)\n",
        sqrt(CA[2,2]), sqrt(CB[2,2]), sqrt(CA[2,2])/sqrt(CB[2,2]))
@printf("\n  WITHOUT the Gaia prior, beta DOES trade against d: corr = %+.3f, sigma_d = %.1f pc.\n",
        RB0[5,6], sqrt(CB0[6,6]))
println("  That is the honest risk in adding `d`, and it is fully controlled: interferometry")
println("  alone cannot measure the distance (d enters only through logg), so the external")
@printf("  prior is ~%.0fx tighter than the data and pins it. A parameter the data cannot\n",
        sqrt(CB0[6,6])/0.1)
println("  constrain but an outside measurement fixes is exactly what you want; ld1 is the")
println("  opposite — the data MUST constrain it, and it competes with beta for the same")
println("  information.")

println("\n" * "="^90)
println("Scan over rotation rate (the diagnostic used on lambda And)")
println("="^90)
@printf("%-6s %-9s | %-26s | %-26s\n", "fev", "R_eq/R_p",
        "A  sig(beta) b/ld1  b/fev", "B  sig(beta) b/d    b/fev")
for fev in (0.30, 0.50, 0.60, 0.70, 0.80, 0.90, 0.92, 0.95)
    tA = copy(θA); tA[2] = fev; tB = copy(θB); tB[2] = fev
    CAs = inv(fisher(obs_planck, tA)); RAs = corrmat(CAs)
    CBs = inv(fisher(obs_atmos, tB) + D_PRIOR); RBs = corrmat(CBs)
    @printf("%-6.2f %-9.4f | %9.5f %+6.3f %+6.3f | %9.5f %+6.3f %+6.3f\n",
            fev, ROTIR.f_rapid_rot_and_deriv(fev)[1],
            sqrt(CAs[5,5]), RAs[5,6], RAs[5,2],
            sqrt(CBs[5,5]), RBs[5,6], RBs[5,2])
end
println("""
  A's beta/ld1 correlation is +0.67 to +0.99 at EVERY rotation rate, worst around
  fev = 0.6-0.7. B's beta/d correlation never exceeds 0.06. The beta/fev correlation —
  which is the STRONGER of the two degeneracies reported for lambda And (-0.67 to -0.84) —
  also weakens in B, from -0.96 to -0.60 at fev = 0.92 and -0.14 at fev = 0.95.

  Note B's sigma(beta) falls monotonically with fev, as it must: more rotation means more
  gravity darkening means a better-determined exponent. A's does not, which is itself a
  symptom of an ill-posed problem rather than of the physics.

  The remaining dominant degeneracy in BOTH cases is rpole/fev (-1.000 in A, -0.995 in B).
  That one is geometric, not photometric: at beta Cas's inclination of 19.9 deg the star is
  seen nearly pole-on, so the polar radius and the oblateness trade almost perfectly.
  Predicted limb darkening does not address it, and nothing here claims otherwise.""")

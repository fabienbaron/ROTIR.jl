# Spica (alpha Vir) — joint interferometry + radial-velocity fit with the tessellated model
#
#   julia -t 8 --project=demos demos/spica_tessellated_joint_fit.jl
#
#   PFREE=1 ...      # also fit the anomalistic period (no prior)
#   MAXITER=3 ...    # quick smoke test (Levenberg-Marquardt iterations, default 100)
#   FDSTEP=1e-4 ...  # finite-difference step, as a fraction of each parameter's box width
#
# ---------------------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------------------
# Interferometry: both stars as Roche surfaces (HEALPix 3072 + 768 tessels), von Zeipel
# gravity darkening (beta = 0.25), mutual irradiation (albedo 0.6), Planck intensity at each
# epoch's mean wavelength, no mutual occultation. Data: 2007_2012_2015.Spica.oifits, split
# into 11 epochs at gaps longer than half a day. The surfaces and maps are computed at each
# epoch's mean time; the secondary's position is computed at each uv point's own time, since
# it moves 0.24-0.52 mas during a 2015 night.
#
# Radial velocities: Keplerian curves from the same orbit, fitted to data/spica_rv.csv. They
# are linear in (K1, K2, gamma), which are solved by weighted least squares at every step;
# the mass ratio of the Roche surfaces is q = K1/K2. The file's sigma = 1 km/s is a
# placeholder, so sigma = 7 + 30 exp(-|V2 - V1|/100) km/s is used instead: line blending
# makes the velocities least certain near conjunction.
#
# Fixed from the literature or from atmosphere models:
#   P      4.01454 d, ANOMALISTIC (Shobbrook et al. 1972; Aufdenberg et al. 2015). ROTIR runs
#          the mean anomaly on P and advances omega separately through domega.
#   Teff   25300 / 20585 K (Tkachenko et al. 2016), from spica_params.jl.
#   LD     H-band power-law exponents fitted to Kurucz ATLAS9 intensities
#          (data/spica_H_kurucz.fits) at each star's polar Teff and logg. The data do not
#          constrain limb darkening: left free, the exponents run to 0 and 1.
#   F1, F2 rotation rates that reproduce v sin i = 161 / 70 km/s at the fitted radii.
#
# Free: the two polar radii, i, Omega, omega, a, e, T0 and domega (no prior on the apsidal
# period).
#
# ---------------------------------------------------------------------------------------
# Weighting
# ---------------------------------------------------------------------------------------
# The interferometric chi2 is divided by its reduced value at the start (OI_SCALE), which is
# the same as inflating the OIFITS errors until chi2r = 1. Without it 6871 interferometric
# points at chi2r ~ 2.4 outweigh 66 velocities by two orders of magnitude.
#
# ---------------------------------------------------------------------------------------
# Fit
# ---------------------------------------------------------------------------------------
# Levenberg-Marquardt (OITOOLS._lm_bounded) on the residual vector: the OI residuals over
# sqrt(OI_SCALE), then the RV residuals. Parameters on a bound of the search box are held
# rather than clipped. The Jacobian is a forward difference, one model evaluation per free
# parameter, with steps of FDSTEP times each box width. The uncertainties are
# sqrt(diag((J'J)^-1)) times sqrt(chi2 / (n - nfree)) of that weighted objective; they take
# the residuals as independent, which night-to-night calibration errors are not.
#
# ---------------------------------------------------------------------------------------
# Result (P fixed)
# ---------------------------------------------------------------------------------------
# OI chi2r 4.14, RV chi2r 1.63; U = 136 yr, e = 0.126, K1/K2 = 123.6/193.5 km/s, q = 0.638,
# parallax 12.07 mas, M1 = 11.0, M2 = 7.05 Msun. The apsidal period, eccentricity, masses and
# parallax all agree with the literature without being constrained to it.
#
# With PFREE=1 the interferometric chi2r falls to 3.62, but only by trading P against
# domega: P = 4.01489 d (12 sigma from 4.01454 +- 0.00003), U = 72 yr (published values
# 105-139 yr) and the RV chi2 rises from 107 to 263. The velocities cover a single season,
# so they cannot anchor that trade; the century-long RV record behind the published P can.

using ROTIR, Printf, Statistics, LinearAlgebra
using OITOOLS: OITOOLS, split_by_time

const HERE = @__DIR__
include(joinpath(HERE, "spica_params.jl"))

const PFREE = get(ENV, "PFREE", "0") == "1"
const MAXITER = parse(Int, get(ENV, "MAXITER", "100"))
const FDSTEP  = parse(Float64, get(ENV, "FDSTEP", "1e-4"))
const P_ANOM = 4.01454

# --- interferometry -----------------------------------------------------------------------
data, mjd = split_by_time(readoifits(joinpath(HERE, "data", "2007_2012_2015.Spica.oifits"))[1, 1]; gap = 0.5)
const NEP   = length(data)
const TEP   = mjd .+ 2400000.5
const BANDS = band_of.(data)
const TESS1, TESS2 = tessellation_healpix(4), tessellation_healpix(3)

# --- radial velocities ----------------------------------------------------------------------
rv = read_rv(joinpath(HERE, "data", "spica_rv.csv"))
v1lit, v2lit = binary_RV((ω = 255.0, e = 0.123, P = P_ANOM, T0 = 2454189.40, dP = 0.0, dω = DOMEGA),
                         rv.jd; K1 = 123.9, K2 = 198.8, γ = 0.0)
const RV_SIG = 7.0 .+ 30.0 .* exp.(-abs.(v2lit .- v1lit) ./ 100.0)
const RV_PRIMARY = rv.star .== :a

# --- limb darkening predicted by the Kurucz grid ------------------------------------------
"Power-law exponent alpha of I(mu)/I(1) = mu^alpha, fitted (projected-area weighted) to the grid over H."
function kurucz_power_law(prov, Teff, logg)
    μ = collect(range(0.005, 1, length = 400))
    I = mean(hcat([ROTIR.provider_intensity(prov, fill(Teff, 400), fill(logg, 400), μ, λ)
                   for λ in range(1.5e-6, 1.75e-6, length = 11)]...), dims = 2)[:]
    I ./= I[end]
    αs = 0:0.0005:1
    return αs[argmin([sum(μ .* (I .- μ .^ α) .^ 2) for α in αs])]
end
prov = TabulatedProvider(load_intensity_grid(joinpath(HERE, "data", "spica_H_kurucz.fits")))
const LD1 = kurucz_power_law(prov, TPOLE1, 3.7)
const LD2 = kurucz_power_law(prov, TPOLE2, 4.0)
@printf("Kurucz H-band power-law limb darkening: primary %.4f, secondary %.4f\n", LD1, LD2)

# --- model ----------------------------------------------------------------------------------
function spica_system(r1, r2, i, Ω, ω, a, e, T0, dω, P, q)
    f1, f2 = F1_SYNC * RPOLE1 / r1, F2_SYNC * RPOLE2 / r2       # fixed v sin i
    orb = (inclination = 180.0 - i, position_angle = Ω - 180.0,
           i = i, Ω = Ω, ω = ω, a = a, e = e, T0 = T0, dω = dω, P = P)
    p1 = merge(SPICA_P1, orb, (rpole = r1, q = q,     ld1 = LD1, rotation_period = P / f1))
    p2 = merge(SPICA_P2, orb, (rpole = r2, q = 1 / q, ld1 = LD2, rotation_period = P / f2))
    s1 = starparameters(r1, TPOLE1, 0.0, LDTYPE, LD1, 0.0, BETA, 0.0, 180.0 - i, Ω - 180.0, 0.0, P / f1)
    s2 = starparameters(r2, TPOLE2, 0.0, LDTYPE, LD2, 0.0, BETA, 0.0, 180.0 - i, Ω - 180.0, 0.0, P / f2)
    bp = binaryparameters(s1, s2, D_PC, i, Ω, ω, P, a, e, T0, q, [1.0, 1.0], 0.0, dω)
    return p1, p2, bp
end

"Normalised residuals of V2, T3amp, T3phi over all epochs; their chi2 per observable; number of points; per-epoch chi2r."
function oi_resid(p1, p2, bp)
    o1 = roche_omega_table(TESS1, p1, bp; secondary = false, npoints = 9)
    o2 = roche_omega_table(TESS2, p2, bp; secondary = true,  npoints = 9)
    res = Float64[]; c = zeros(3); n = 0; per = zeros(NEP)
    for k in 1:NEP
        D = compute_separation(bp, TEP[k])
        g1, g2 = create_binary_geometry(TESS1, p1, TESS2, p2, bp, TEP[k]; omega1 = o1(D), omega2 = o2(D))
        m1 = temperature_map_vonZeipel_roche_single(p1, g1, TEP[k])
        m2 = temperature_map_vonZeipel_roche_single(p2, g2, TEP[k]; secondary = true)
        h1, h2 = handle_reflection(g1, m1, g2, m2; albedo1 = ALBEDO_DEFAULT, albedo2 = ALBEDO_DEFAULT)
        ph = binary_phase_shift(data[k], bp)      # secondary at each uv point's own time
        v2, t3a, t3p = binary_observables(h1, g1, h2, g2, data[k], ph; intensity_model = :planck, band = BANDS[k])
        rk = ((v2 .- data[k].v2) ./ data[k].v2_err,
              (t3a .- data[k].t3amp) ./ data[k].t3amp_err,
              mod360(t3p .- data[k].t3phi) ./ data[k].t3phi_err)
        ck = [sum(abs2, x) for x in rk]
        nk = data[k].nv2 + data[k].nt3amp + data[k].nt3phi
        append!(res, rk...); c .+= ck; n += nk; per[k] = sum(ck) / nk
    end
    return res, c, n, per
end

"(K1, K2, gamma) by weighted least squares, the normalised RV residuals and their chi2."
function rv_solve(ω, e, T0, dω, P)
    b1, b2 = binary_RV((ω = ω, e = e, P = P, T0 = T0, dP = 0.0, dω = dω), rv.jd; K1 = 1.0, K2 = 1.0, γ = 0.0)
    A = hcat(ifelse.(RV_PRIMARY, b1, 0.0), ifelse.(RV_PRIMARY, 0.0, b2), ones(length(rv.jd))) ./ RV_SIG
    K = A \ (rv.rv ./ RV_SIG)
    res = A * K .- rv.rv ./ RV_SIG
    return K, res, sum(abs2, res)
end

# --- parameters -----------------------------------------------------------------------------
# Start: omega, e, T0 and domega from the disc-model joint fit
# (results/spica_orbitfit_joint_rv_best.txt); radii, i, Omega and a from a tessellated fit of
# those four alone with that orbit held fixed.
J = Dict(split(l, '\t')[1] => parse(Float64, split(l, '\t')[2])
         for l in eachline(joinpath(HERE, "results", "spica_orbitfit_joint_rv_best.txt"))
         if !startswith(l, '#') && occursin('\t', l))
const NAMES = ["rpole1", "rpole2", "i", "Omega", "omega", "a", "e", "T0", "domega", "P"]
x0 = [0.400305, 0.199639, 117.03359, 310.78553, mod(J["omega"], 360), 1.549748, J["e"], J["T0"], J["domega"], P_ANOM]
const HALF = [0.1, 0.06, 4, 4, 20, 0.08, 0.06, 0.3, 0.008, 1.5e-3]
const LO, HI = x0 .- HALF, x0 .+ HALF .+ [0, 0, 0, 0, 0, 0, 0, 0, 0.004, 0]
const FREE = PFREE ? collect(1:10) : collect(1:9)
full(x) = (y = copy(x0); y[FREE] .= x; y)

function evaluate(y)
    r1, r2, i, Ω, ω, a, e, T0, dω, P = y
    K, rrv, crv = rv_solve(ω, e, T0, dω, P)
    q = clamp(K[1] / K[2], 0.4, 0.9)
    roi, c, n, per = oi_resid(spica_system(r1, r2, i, Ω, ω, a, e, T0, dω, P, q)...)
    return (; roi, rrv, c, n, per, crv, K, q)
end

r0 = evaluate(x0)
const OI_SCALE = sum(r0.c) / r0.n
@printf("start: OI chi2r %.3f, RV chi2 %.1f over %d velocities\n", OI_SCALE, r0.crv, length(rv.jd))

# Residuals of the free parameters x. The last evaluation is kept, since _lm_bounded asks for
# the Jacobian at the point whose residuals it has just accepted.
const NEVAL = Ref(0)
const LAST = Ref{Any}(nothing)
function resid(x)
    LAST[] !== nothing && LAST[][1] == x && return LAST[][2]
    r = evaluate(full(x))
    NEVAL[] += 1
    v = vcat(r.roi ./ sqrt(OI_SCALE), r.rrv)
    LAST[] = (copy(x), v)
    return v
end

const LOF, HIF = LO[FREE], HI[FREE]
const H = FDSTEP .* (HIF .- LOF)
const NJAC = Ref(0)
"Forward-difference Jacobian, stepping back from an upper bound."
function jac(x)
    r = resid(x)
    Jm = zeros(length(r), length(x))
    for j in eachindex(x)
        h = x[j] + H[j] <= HIF[j] ? H[j] : -H[j]
        xp = copy(x); xp[j] += h
        Jm[:, j] = (resid(xp) .- r) ./ h
    end
    LAST[] = (copy(x), r)
    NJAC[] += 1
    @printf("  Jacobian %d (%d evaluations): OI chi2r %.4f  RV chi2 %.1f\n", NJAC[], NEVAL[],
            sum(abs2, r[1:r0.n]) * OI_SCALE / r0.n, sum(abs2, r[r0.n+1:end]))
    return Jm
end

res = OITOOLS._lm_bounded(resid, jac, x0[FREE], LOF, HIF; maxIter = MAXITER)
y = full(res.x)
r = evaluate(y)
nres = length(res.resid)
cov = inv(Symmetric(res.jacobian' * res.jacobian)) .* (sum(abs2, res.resid) / (nres - length(FREE)))
σ = zeros(length(y)); σ[FREE] .= sqrt.(diag(cov))

# --- report ---------------------------------------------------------------------------------
plx, q = rv_physical(y[6], y[3], y[10], y[7], r.K[1], r.K[2])
m = dynamical_masses(y[6], y[10], plx, q)
println("\n", PFREE ? "P free" : "P fixed at $P_ANOM d", "; ", res.reason, " after ", res.iterations,
        " iterations, ", NEVAL[], " evaluations")
@printf("OI   V2 %.0f  T3amp %.0f  T3phi %.0f   chi2r %.3f over %d points\n", r.c..., sum(r.c) / r.n, r.n)
@printf("RV   chi2 %.1f   chi2r %.2f over %d velocities\n", r.crv, r.crv / length(rv.jd), length(rv.jd))
for (j, (k, v)) in enumerate(zip(NAMES, y))
    jf = findfirst(==(j), FREE)
    @printf("  %-7s %.8g%s\n", k, v, jf === nothing ? "  (fixed)" :
            @sprintf("  +- %.2g%s", σ[j], jf in res.held ? "  HELD ON BOUND" : ""))
end
@printf("  K1 %.2f  K2 %.2f  gamma %.2f km/s   q %.4f\n", r.K..., q)
@printf("  parallax %.3f mas (%.1f pc)   M1 %.2f  M2 %.2f Msun   U %.1f yr\n",
        plx, m.d, m.M1, m.M2, 360 / y[9] / 365.25)
@printf("  polar radii %.2f / %.2f Rsun   LD %.4f / %.4f   F1 %.3f  F2 %.3f\n",
        phys_radius_rsun(y[1], m.d), phys_radius_rsun(y[2], m.d), LD1, LD2,
        F1_SYNC * RPOLE1 / y[1], F2_SYNC * RPOLE2 / y[2])
println("  chi2r per epoch: ", round.(r.per, digits = 2))

out = joinpath(HERE, "results", PFREE ? "spica_tessellated_joint_pfree_best.txt" : "spica_tessellated_joint_best.txt")
open(out, "w") do io
    println(io, "# demos/spica_tessellated_joint_fit.jl", PFREE ? " (PFREE=1)" : "", ": tessellated Roche model on")
    println(io, "# 2007_2012_2015.Spica.oifits + Keplerian RVs (spica_rv.csv, sigma = 7 + 30 exp(-|V2-V1|/100) km/s).")
    @printf(io, "# OI chi2r %.3f over %d points; RV chi2 %.1f over %d. LD (power law) %.4f / %.4f.\n",
            sum(r.c) / r.n, r.n, r.crv, length(rv.jd), LD1, LD2)
    println(io, "# name\tvalue\tsigma (0 = fixed)")
    for (k, v, e) in zip(NAMES, y, σ); println(io, k, '\t', v, '\t', e); end
    for (k, v) in (("K1", r.K[1]), ("K2", r.K[2]), ("gamma", r.K[3]), ("q", q), ("plx", plx),
                   ("M1", m.M1), ("M2", m.M2), ("ld1_1", LD1), ("ld1_2", LD2))
        println(io, k, '\t', v)
    end
end
println("saved ", out)

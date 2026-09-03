# spica_irradiation_detect.jl
#
# Can CHARA interferometry of Spica detect mutual irradiation?
#
# Run:  julia --project=demos demos/spica_irradiation_detect.jl
#
# Four steps:
#   1. Forward model at the TRUE epoch times (Roche shape from D(t), von Zeipel, Planck
#      intensities), with an audit of the assumptions that could invalidate the answer.
#   2. Δχ² against the real data over a grid of bolometric albedos, split by observable,
#      plus the model-difference amplitude |Δmodel|/σ — which says whether the data
#      *could* see the effect independently of how well the rest of the model fits.
#   3. Injection–recovery: synthesize noiseless data at A = 0.6, then evaluate at A = 0,
#      (a) with everything else fixed — the optimistic bound — and (b) re-optimising the
#      nuisance parameters, which is the honest number, because a wrong albedo is largely
#      absorbable by slightly different radii and a slightly different flux ratio.
#   4. The (β, albedo) degeneracy map: gravity darkening COOLS the sub-companion point
#      (it is the tidal bulge tip, lowest local gravity) while irradiation HEATS it, so a
#      significance quoted at fixed β is not a detection of irradiation.
#
# Knobs: NSIDE1=4 NSIDE2=3 NSIDE_FIT1=3 NSIDE_FIT2=2 VCONS=1 MAXFIT=400 OUTDIR=...
#        OCC=exact|soft|false   mutual occultation
#        ASYNC=1|0  APSIDAL=1|0  (see spica_params.jl) — both ON by default, since Spica
#        is neither synchronous (F = 1.92/1.65) nor free of apsidal motion (U = 139 yr).

ENV["MPLBACKEND"] = get(ENV, "MPLBACKEND", "Agg")
using ROTIR, PythonPlot, Printf, Statistics, LinearAlgebra
using Optim
include(joinpath(@__DIR__, "spica_params.jl"))

nside1  = parse(Int, get(ENV, "NSIDE1", "4"))
nside2  = parse(Int, get(ENV, "NSIDE2", "3"))
nsidef1 = parse(Int, get(ENV, "NSIDE_FIT1", "3"))   # coarser mesh for the nuisance re-fit
nsidef2 = parse(Int, get(ENV, "NSIDE_FIT2", "2"))
vcons   = get(ENV, "VCONS", "1") == "1"
# Mutual occultation: :exact area clipping, :soft centre test, or false.
const OCC = let o = get(ENV, "OCC", "exact"); o == "false" ? false : Symbol(o) end
maxfit  = parse(Int, get(ENV, "MAXFIT", "400"))
outdir  = get(ENV, "OUTDIR", joinpath(@__DIR__, "results", "spica_detect"))
mkpath(outdir)

# =======================================================================================
# 1. Data and geometry
# =======================================================================================
spica_audit()

oifitsfile = joinpath(@__DIR__, "data", "2007_2012_2015.Spica.oifits")
data_all = readoifits(oifitsfile)[1, 1]

# Split into epochs by MJD gaps (same rule as demos/spica_binary_roche.jl).
mjds = sort(data_all.v2_mjd)
jumps = findall(diff(mjds) .> 0.5)
starts = mjds[[1; jumps .+ 1]]
stops  = mjds[[jumps; length(mjds)]]
nepochs = length(starts)
data = Vector{typeof(data_all)}(undef, nepochs)
tepochs = zeros(nepochs)          # JD
for i in 1:nepochs
    idx = set_data_filter(data_all; mjd_range=[starts[i] - 0.01, stops[i] + 0.01])
    data[i] = filter_data(data_all, idx)
    tepochs[i] = mean(data[i].v2_mjd) + 2400000.5
end
bands = [band_of(d) for d in data]

# ---------------------------------------------------------------------------------------
# BAND: which combiner are we pretending to be?
# ---------------------------------------------------------------------------------------
# BAND=native (default) uses the wavelengths in the file — CHARA/MIRC-X in H, ~1.61 um.
# BAND=<microns> re-flies the SAME PHYSICAL BASELINES at a different wavelength, which is
# what CHARA/SPICA is: the same array, a visible-light combiner (~0.6-0.9 um).
#
# `uv` is stored as metres/wavelength, i.e. spatial frequency, so a shorter wavelength on
# the same baseline is simply uv * (lambda_native / lambda_new). Errors are left at their
# measured FRACTIONAL values, so this isolates the effect of RESOLUTION and of the Planck
# lever arm — it is not a prediction of SPICA's photometric performance.
#
# Why this matters here: the proximity effects scale with wavelength through two channels
# that pull in opposite directions in importance.
#
#   1. THERMAL. Irradiation and gravity darkening are temperature perturbations, and a
#      given dT/T produces a flux contrast dlnB/dlnT = x e^x/(e^x - 1) with x = hc/(lambda k T).
#      At 20585 K this is 1.23 in H and 1.58 at 0.7 um — only a 1.28x gain. Visible light
#      is NOT deep in the Wien regime for a 20-25 kK star.
#   2. SPATIAL. Resolution is lambda/2B, so 1.61 -> 0.65 um is 2.5x finer. Decisive, because
#      the primary's first visibility null moves from 449 m (beyond CHARA's 331 m, so the
#      curve never turns over) to 181 m (well inside). Only past the null do diameter, limb
#      darkening and tidal elongation stop trading against one another — and that
#      degeneracy, not sensitivity, is what limited the H-band result to 1.3 sigma.
#
# So the expected gain is much larger than the 1.28x thermal factor alone.
const BAND = get(ENV, "BAND", "native")
const BANDLABEL = if BAND == "native"
    @sprintf("H band, MIRC-X (native, %.3f um)", 1e6mean(bands))
else
    lam_new = parse(Float64, BAND) * 1e-6
    for i in 1:nepochs
        data[i].uv .*= bands[i] / lam_new       # same baselines, new wavelength
        bands[i] = lam_new
    end
    @sprintf("R band, SPICA-like (%.3f um, same CHARA baselines)", 1e6lam_new)
end
println("\nBAND: ", BANDLABEL)

ndata = [d.nv2 + d.nt3amp + d.nt3phi for d in data]

println("\n$nepochs epochs, $(sum(ndata)) data points, " *
        "λ = $(round(1e6*minimum(bands), digits=3))–$(round(1e6*maximum(bands), digits=3)) μm")
println("NOTE: ROTIR collapses the H-band channels of each epoch into one block, so the")
println("      model uses one representative λ per epoch (the spread across the band is")
println("      1.49–1.77 μm, a ±8% lever arm on the Planck ratio that is not exploited).")

tes1 = tessellation_healpix(nside1)
tes2 = tessellation_healpix(nside2)

"""
    build(p1, p2, tes1, tes2, ts; vcons) -> (stars1, stars2, kernels, phases)

Geometry, reflection kernels and Fourier phase shifts for every epoch. Everything that
does not depend on the albedo or on β lives here so the scans below only redo the cheap
part.
"""
function build(p1, p2, t1::tessellation, t2::tessellation, ts; vcons=true, npoints=9)
    o1t = vcons ? roche_omega_table(t1, p1, SPICA_BP; secondary=false, npoints=npoints) : nothing
    o2t = vcons ? roche_omega_table(t2, p2, SPICA_BP; secondary=true,  npoints=npoints) : nothing
    n = length(ts)
    # Concretely typed (as create_star_multiepochs does): a Vector{Any} here makes the
    # threaded loop in setup_oi! type-unstable.
    s1 = Array{stellar_geometry}(undef, n); s2 = Array{stellar_geometry}(undef, n)
    for i in 1:n
        D = compute_separation(SPICA_BP, ts[i])
        s1[i], s2[i] = create_binary_geometry(t1, p1, t2, p2, SPICA_BP, ts[i];
                                              omega1 = o1t === nothing ? nothing : o1t(D),
                                              omega2 = o2t === nothing ? nothing : o2t(D))
    end
    return s1, s2
end

@info "building geometry at the true epoch times (volume-conserving = $vcons)…"
stars1, stars2 = build(SPICA_P1, SPICA_P2, tes1, tes2, tepochs; vcons=vcons)
setup_oi!(data, stars1)
setup_oi!(data, stars2)

phases = [binary_phase_shift(data[i].uv, orbit_to_rotir_offset(SPICA_BP, tepochs[i])...)
          for i in 1:nepochs]

println()
overlap = check_binary_overlap(stars1, stars2, SPICA_BP, tepochs; verbose=true)
if any(overlap)
    println("\n*** $(count(overlap))/$nepochs epochs have overlapping disks (Spica is a known")
    println("    grazing eclipser, Desmet+2009). Mutual occultation is $(OCC === false ? "OFF" : "ON ($(OCC))");")
    println("    with it off, the far component's hidden face would keep contributing flux —")
    println("    ~10% of the secondary's light, ~1.7% of the total, at closest approach.\n")
end

# Intrinsic (gravity-darkened) maps, and the reflection kernels, per epoch.
tm1 = [parametric_temperature_map(SPICA_P1, stars1[i]) for i in 1:nepochs]
tm2 = [parametric_temperature_map(SPICA_P2, stars2[i]; secondary=true) for i in 1:nepochs]
kern = [reflection_kernels(stars1[i], stars2[i]) for i in 1:nepochs]

@printf("intrinsic Tmap ranges: primary %.0f–%.0f K, secondary %.0f–%.0f K\n",
        minimum(minimum.(tm1)), maximum(maximum.(tm1)),
        minimum(minimum.(tm2)), maximum(maximum.(tm2)))

# =======================================================================================
# 2. Δχ² against the real data, as a function of albedo
# =======================================================================================
"Heated maps at albedo `A`, reusing the per-epoch kernels."
function heat(i, A; t1=tm1, t2=tm2, k=kern)
    A == 0 && return t1[i], t2[i]
    handle_reflection(stars1[i], t1[i], stars2[i], t2[i];
                      albedo1=A, albedo2=A, kernels=k[i])
end

"Per-observable χ² of the whole dataset at albedo `A`."
function chi2_split(A; imodel=:planck)
    cv2 = ca = cp = 0.0
    for i in 1:nepochs
        h1, h2 = heat(i, A)
        v2, t3a, t3p = binary_observables(h1, stars1[i], h2, stars2[i], data[i], phases[i];
                                          intensity_model=imodel, band=bands[i],
                                          occultation=OCC)
        cv2 += sum(abs2, (v2 .- data[i].v2) ./ data[i].v2_err)
        ca  += sum(abs2, (t3a .- data[i].t3amp) ./ data[i].t3amp_err)
        cp  += sum(abs2, mod360(t3p .- data[i].t3phi) ./ data[i].t3phi_err)
    end
    return cv2, ca, cp
end

# =======================================================================================
# 1b. Ablation: what does each newly-available systematic actually do?
# =======================================================================================
# Each of these is larger than the irradiation signature, so a detectability number
# computed without them is not meaningful. Rebuild the whole model under each
# configuration and report both the absolute fit and the irradiation Δχ².
println("\n" * "="^78)
println("1b. Ablation over the systematics")
println("="^78)
println("Each row rebuilds the geometry from scratch. `Δχ²(A)` is χ²(A=0) − χ²(A=0.6),")
println("i.e. how much the data would prefer irradiation under that model.\n")

function full_chi2(p1, p2, bp; occ = false, A = 0.0)
    o1t = vcons ? roche_omega_table(tes1, p1, bp; secondary=false, npoints=9) : nothing
    o2t = vcons ? roche_omega_table(tes2, p2, bp; secondary=true,  npoints=9) : nothing
    c = 0.0
    for i in 1:nepochs
        D = compute_separation(bp, tepochs[i])
        a, b = create_binary_geometry(tes1, p1, tes2, p2, bp, tepochs[i];
                                      omega1 = o1t === nothing ? nothing : o1t(D),
                                      omega2 = o2t === nothing ? nothing : o2t(D))
        setup_oi!([data[i]], [a]); setup_oi!([data[i]], [b])
        m1 = parametric_temperature_map(p1, a)
        m2 = parametric_temperature_map(p2, b; secondary=true)
        if A > 0
            m1, m2 = handle_reflection(a, m1, b, m2; albedo1=A, albedo2=A)
        end
        pshift = binary_phase_shift(data[i].uv, orbit_to_rotir_offset(bp, tepochs[i])...)
        c += binary_chi2_f(m1, a, m2, b, data[i], pshift;
                           intensity_model=:planck, band=bands[i],
                           occultation = occ ? :exact : false)
    end
    return c
end

P1sync = merge(SPICA_P1, (rotation_period = P_ORB,))
P2sync = merge(SPICA_P2, (rotation_period = P_ORB,))
BPnoap = merge(SPICA_BP, (dω = 0.0,))
P1noap = merge(SPICA_P1, (dω = 0.0,)); P2noap = merge(SPICA_P2, (dω = 0.0,))
P1both = merge(P1sync,   (dω = 0.0,)); P2both = merge(P2sync,   (dω = 0.0,))

configs = (("baseline (as originally run)", P1both, P2both, BPnoap, false),
           ("+ apsidal motion dω",          P1sync, P2sync, SPICA_BP, false),
           ("+ asynchronous rotation",      P1noap, P2noap, BPnoap,  false),
           ("+ mutual occultation",         P1both, P2both, BPnoap,  true),
           ("ALL THREE",                    SPICA_P1, SPICA_P2, SPICA_BP, true))
@printf("  %-30s %12s %10s %12s\n", "configuration", "χ²/n", "vs base", "Δχ²(A=0.6)")
base_c = NaN
for (lab, p1, p2, bp, occ) in configs
    c0 = full_chi2(p1, p2, bp; occ=occ, A=0.0)
    cA = full_chi2(p1, p2, bp; occ=occ, A=0.6)
    isnan(base_c) && (global base_c = c0)
    @printf("  %-30s %12.2f %+10.0f %12.1f\n",
            lab, c0/sum(ndata), c0 - base_c, c0 - cA)
end
println("\n(a POSITIVE Δχ² means the data prefer irradiation; negative means they disprefer it)")

println("\n" * "="^78)
println("2. χ² against the real Spica data vs bolometric albedo")
println("="^78)
albedos = [0.0, 0.3, 0.6, 1.0]
println("  A      χ²ᵥ₂/n     χ²ₜ₃ₐ/n    χ²ₜ₃ₚ/n    χ²/n      Δχ² vs A=0")
chi2_0 = sum(chi2_split(0.0))
for A in albedos
    c = chi2_split(A)
    @printf("  %.1f   %8.3f   %8.3f   %8.3f   %8.3f   %+10.1f\n",
            A, c[1]/sum(d.nv2 for d in data), c[2]/sum(d.nt3amp for d in data),
            c[3]/sum(d.nt3phi for d in data), sum(c)/sum(ndata), sum(c) - chi2_0)
end
println("\n(The absolute χ²/n is large because the orbital elements, radii and temperatures")
println(" are literature values, not a fit to this dataset. What matters below is Δχ².)")

# How big is the irradiation signature compared with the error bars?
println("\nModel-difference amplitude at A = 0.6 (rms of |Δmodel|/σ, and max):")
sv = Float64[]; sa = Float64[]; sp = Float64[]
for i in 1:nepochs
    h1, h2 = heat(i, 0.6)
    v2a, taa, tpa = binary_observables(h1, stars1[i], h2, stars2[i], data[i], phases[i];
                                       intensity_model=:planck, band=bands[i])
    v2b, tab, tpb = binary_observables(tm1[i], stars1[i], tm2[i], stars2[i], data[i], phases[i];
                                       intensity_model=:planck, band=bands[i])
    append!(sv, abs.(v2a .- v2b) ./ data[i].v2_err)
    append!(sa, abs.(taa .- tab) ./ data[i].t3amp_err)
    append!(sp, abs.(mod360(tpa .- tpb)) ./ data[i].t3phi_err)
end
@printf("  V²    rms %.3f σ   max %.3f σ\n", sqrt(mean(abs2, sv)), maximum(sv))
@printf("  T3amp rms %.3f σ   max %.3f σ\n", sqrt(mean(abs2, sa)), maximum(sa))
@printf("  T3phi rms %.3f σ   max %.3f σ\n", sqrt(mean(abs2, sp)), maximum(sp))

# =======================================================================================
# 3. Injection–recovery
# =======================================================================================
# Noiseless injection: Σ((model_A − model_0)/σ)² IS the expected Δχ² of the noisy
# statistic (its mean), so it measures the detectability without the scatter of one
# particular realisation.
println("\n" * "="^78)
println("3. Injection–recovery (inject A = 0.6, fit with A = 0)")
println("="^78)

A_inj = 0.6
synth = [deepcopy(d) for d in data]
for i in 1:nepochs
    h1, h2 = heat(i, A_inj)
    v2, t3a, t3p = binary_observables(h1, stars1[i], h2, stars2[i], data[i], phases[i];
                                      intensity_model=:planck, band=bands[i])
    synth[i].v2    .= v2
    synth[i].t3amp .= t3a
    synth[i].t3phi .= t3p
end

function chi2_vs(dat, A; s1=stars1, s2=stars2, t1=tm1, t2=tm2, k=kern, ph=phases, bd=bands)
    c = 0.0
    for i in eachindex(dat)
        h1, h2 = A == 0 ? (t1[i], t2[i]) :
                 handle_reflection(s1[i], t1[i], s2[i], t2[i]; albedo1=A, albedo2=A, kernels=k[i])
        v2, t3a, t3p = binary_observables(h1, s1[i], h2, s2[i], dat[i], ph[i];
                                          intensity_model=:planck, band=bd[i], occultation=OCC)
        c += sum(abs2, (v2 .- dat[i].v2) ./ dat[i].v2_err)
        c += sum(abs2, (t3a .- dat[i].t3amp) ./ dat[i].t3amp_err)
        c += sum(abs2, mod360(t3p .- dat[i].t3phi) ./ dat[i].t3phi_err)
    end
    return c
end

dchi2_fixed = chi2_vs(synth, 0.0) - chi2_vs(synth, A_inj)
@printf("\n(a) everything else fixed:   Δχ² = %.1f   →  %.1f σ (1 dof)\n",
        dchi2_fixed, sqrt(max(dchi2_fixed, 0)))

# Is that above the mesh-discretisation floor? Injection and recovery share a mesh, so the
# discretisation error largely cancels — but the VALUE of Δχ² still has to converge, or
# the significance is an artefact of the tessellation. Repeat the whole thing at two
# coarser resolutions and check the trend.
println("\n    mesh convergence of Δχ² (inject and recover on the same mesh each time):")
@printf("      nside %d/%d  (%5d/%5d tessels)   Δχ² = %7.2f\n",
        nside1, nside2, tes1.npix, tes2.npix, dchi2_fixed)
for (m1, m2) in ((nside1 - 1, nside2 - 1), (nside1 - 2, nside2 - 2))
    (m1 >= 1 && m2 >= 0) || continue
    u1 = tessellation_healpix(m1); u2 = tessellation_healpix(m2)
    v1, v2s = build(SPICA_P1, SPICA_P2, u1, u2, tepochs; vcons=vcons)
    setup_oi!(data, v1); setup_oi!(data, v2s)
    q1 = [parametric_temperature_map(SPICA_P1, v1[i]) for i in 1:nepochs]
    q2 = [parametric_temperature_map(SPICA_P2, v2s[i]; secondary=true) for i in 1:nepochs]
    kq = [reflection_kernels(v1[i], v2s[i]) for i in 1:nepochs]
    sy = [deepcopy(d) for d in data]
    for i in 1:nepochs
        h1, h2 = handle_reflection(v1[i], q1[i], v2s[i], q2[i];
                                   albedo1=A_inj, albedo2=A_inj, kernels=kq[i])
        a, b, c = binary_observables(h1, v1[i], h2, v2s[i], data[i], phases[i];
                                     intensity_model=:planck, band=bands[i])
        sy[i].v2 .= a; sy[i].t3amp .= b; sy[i].t3phi .= c
    end
    d0 = chi2_vs(sy, 0.0;     s1=v1, s2=v2s, t1=q1, t2=q2, k=kq)
    dA = chi2_vs(sy, A_inj;   s1=v1, s2=v2s, t1=q1, t2=q2, k=kq)
    @printf("      nside %d/%d  (%5d/%5d tessels)   Δχ² = %7.2f\n",
            m1, m2, u1.npix, u2.npix, d0 - dA)
end
println("      ↑ if these are not settling, raise NSIDE1/NSIDE2 before quoting a significance.")

# (b) re-optimise the nuisance parameters at A = 0. Coarser mesh: this needs a full
# geometry + polygon-FT rebuild per objective evaluation.
println("\n(b) re-optimising nuisance parameters at A = 0 (rpole₁, rpole₂, tpole₂, ld1₁, ld1₂)")
println("    on a coarser mesh (nside $nsidef1/$nsidef2) — this is the honest number, since a")
println("    wrong albedo is partly absorbable by a different radius and flux ratio.")

tf1 = tessellation_healpix(nsidef1)
tf2 = tessellation_healpix(nsidef2)

# Re-inject on the fit mesh so the comparison is mesh-consistent (otherwise the
# discretisation difference between the two meshes leaks into Δχ²).
sf1, sf2 = build(SPICA_P1, SPICA_P2, tf1, tf2, tepochs; vcons=vcons)
setup_oi!(data, sf1); setup_oi!(data, sf2)
tf1m = [parametric_temperature_map(SPICA_P1, sf1[i]) for i in 1:nepochs]
tf2m = [parametric_temperature_map(SPICA_P2, sf2[i]; secondary=true) for i in 1:nepochs]
kf   = [reflection_kernels(sf1[i], sf2[i]) for i in 1:nepochs]
synth_f = [deepcopy(d) for d in data]
for i in 1:nepochs
    h1, h2 = handle_reflection(sf1[i], tf1m[i], sf2[i], tf2m[i];
                               albedo1=A_inj, albedo2=A_inj, kernels=kf[i])
    v2, t3a, t3p = binary_observables(h1, sf1[i], h2, sf2[i], data[i], phases[i];
                                      intensity_model=:planck, band=bands[i])
    synth_f[i].v2 .= v2; synth_f[i].t3amp .= t3a; synth_f[i].t3phi .= t3p
end

θ0 = [RPOLE1, RPOLE2, TPOLE2, LD1_H, LD1_H]
"χ² of the injected data for nuisance vector θ at albedo A."
function chi2_nuisance(θ, A)
    all(isfinite, θ) || return 1e30
    (θ[1] > 0 && θ[2] > 0 && θ[3] > 3000 && 0 <= θ[4] <= 1 && 0 <= θ[5] <= 1) || return 1e30
    p1 = merge(SPICA_P1, (rpole = θ[1], ld1 = θ[4]))
    p2 = merge(SPICA_P2, (rpole = θ[2], tpole = θ[3], ld1 = θ[5]))
    s1, s2 = build(p1, p2, tf1, tf2, tepochs; vcons=vcons, npoints=7)
    setup_oi!(data, s1); setup_oi!(data, s2)
    t1 = [parametric_temperature_map(p1, s1[i]) for i in 1:nepochs]
    t2 = [parametric_temperature_map(p2, s2[i]; secondary=true) for i in 1:nepochs]
    k  = A == 0 ? nothing : [reflection_kernels(s1[i], s2[i]) for i in 1:nepochs]
    c = 0.0
    for i in 1:nepochs
        h1, h2 = A == 0 ? (t1[i], t2[i]) :
                 handle_reflection(s1[i], t1[i], s2[i], t2[i]; albedo1=A, albedo2=A, kernels=k[i])
        v2, t3a, t3p = binary_observables(h1, s1[i], h2, s2[i], synth_f[i], phases[i];
                                          intensity_model=:planck, band=bands[i], occultation=OCC)
        c += sum(abs2, (v2 .- synth_f[i].v2) ./ synth_f[i].v2_err)
        c += sum(abs2, (t3a .- synth_f[i].t3amp) ./ synth_f[i].t3amp_err)
        c += sum(abs2, mod360(t3p .- synth_f[i].t3phi) ./ synth_f[i].t3phi_err)
    end
    return c
end

# Scale the simplex so all five parameters move comparably.
scale = [0.02 * RPOLE1, 0.02 * RPOLE2, 400.0, 0.05, 0.05]
obj0 = u -> chi2_nuisance(θ0 .+ scale .* u, 0.0)
@info "running Nelder–Mead (up to $maxfit evaluations)…"
res0 = optimize(obj0, zeros(5), NelderMead(),
                Optim.Options(iterations=maxfit, show_trace=false))
θfit = θ0 .+ scale .* Optim.minimizer(res0)
chi2_A0_refit = Optim.minimum(res0)
chi2_truth    = chi2_nuisance(θ0, A_inj)     # should be ≈ 0 (noiseless self-consistency)

@printf("\n    χ² at the injected truth (A=0.6, θ=θ₀)      %12.3f   [≈0 confirms self-consistency]\n", chi2_truth)
@printf("    χ² at A=0, θ=θ₀ (no re-fit)                 %12.1f\n", chi2_nuisance(θ0, 0.0))
@printf("    χ² at A=0 after re-fitting nuisance params  %12.1f\n", chi2_A0_refit)
dchi2_refit = chi2_A0_refit - chi2_truth
@printf("    Δχ² (honest)                                %12.1f   →  %.1f σ (1 dof)\n",
        dchi2_refit, sqrt(max(dchi2_refit, 0)))
println("\n    best-fit nuisance parameters when irradiation is (wrongly) omitted:")
@printf("      rpole₁ %.5f → %.5f mas (%+.2f%%,  %.3f → %.3f R☉)\n",
        RPOLE1, θfit[1], 100*(θfit[1]/RPOLE1 - 1), R1_RSUN, phys_radius_rsun(θfit[1], D_PC))
@printf("      rpole₂ %.5f → %.5f mas (%+.2f%%,  %.3f → %.3f R☉)\n",
        RPOLE2, θfit[2], 100*(θfit[2]/RPOLE2 - 1), R2_RSUN, phys_radius_rsun(θfit[2], D_PC))
@printf("      tpole₂ %.0f → %.0f K (%+.0f K)\n", TPOLE2, θfit[3], θfit[3] - TPOLE2)
@printf("      ld1₁   %.3f → %.3f      ld1₂ %.3f → %.3f\n", LD1_H, θfit[4], LD1_H, θfit[5])
println("    ↑ how much of the irradiation signature these absorb is exactly the difference")
println("      between the optimistic and the honest Δχ².")

# =======================================================================================
# 4. (β, albedo) degeneracy
# =======================================================================================
println("\n" * "="^78)
println("4. (β, albedo) Δχ² map on the injected data")
println("="^78)
println("Gravity darkening cools the sub-companion point (tidal bulge tip = lowest g)")
println("while irradiation heats it, so a significance quoted at fixed β has to be checked")
println("against the joint surface. β only changes the temperature maps, not the geometry,")
println("so this scan reuses the kernels.")

βs = collect(range(0.05, 0.35, length=7))
As = collect(range(0.0, 1.0, length=9))
Z = fill(NaN, length(βs), length(As))
for (ib, β) in enumerate(βs)
    p1β = merge(SPICA_P1, (beta = β,))
    p2β = merge(SPICA_P2, (beta = β,))
    t1β = [parametric_temperature_map(p1β, stars1[i]) for i in 1:nepochs]
    t2β = [parametric_temperature_map(p2β, stars2[i]; secondary=true) for i in 1:nepochs]
    for (ia, A) in enumerate(As)
        Z[ib, ia] = chi2_vs(synth, A; t1=t1β, t2=t2β)
    end
end
Zmin = minimum(Z)
@printf("\nminimum Δχ² = 0 at β = %.3f, A = %.2f  (injected: β = %.2f, A = %.2f)\n",
        βs[argmin(Z)[1]], As[argmin(Z)[2]], BETA, A_inj)
print("\n   β \\ A "); for A in As; @printf("%8.2f", A); end; println()
for (ib, β) in enumerate(βs)
    @printf("  %.3f  ", β); for ia in eachindex(As); @printf("%8.1f", Z[ib, ia] - Zmin); end; println()
end

fig, ax = subplots(figsize=(8, 6))
lev = [1.0, 4.0, 9.0, 25.0, 100.0, 400.0]
cs = ax.contour(As, βs, Z .- Zmin, levels=lev, colors="k")
ax.clabel(cs, inline=true, fmt="%.0f")
im = ax.pcolormesh(As, βs, Z .- Zmin, shading="auto", cmap="viridis")
ax.plot([A_inj], [BETA], "r*", markersize=16, label="injected")
ax.set_xlabel("bolometric albedo A"); ax.set_ylabel("von Zeipel β")
ax.set_title("Δχ² on injected Spica data — $BANDLABEL\n(contours at 1, 4, 9, 25, 100, 400)")
ax.legend(); fig.colorbar(im, ax=ax, label="Δχ²")
_tag = BAND == "native" ? "H" : "R" * replace(BAND, "." => "p")
fig.savefig(joinpath(outdir, "beta_albedo_dchi2_$(_tag).png"), dpi=130, bbox_inches="tight")
pyplot.close(fig)

# ---- summary --------------------------------------------------------------------------
println("\n" * "="^78)
println("SUMMARY")
println("="^78)
@printf("  irradiation signature vs errors : V² %.2fσ rms, T3φ %.2fσ rms\n",
        sqrt(mean(abs2, sv)), sqrt(mean(abs2, sp)))
@printf("  Δχ² for A=0.6 vs A=0, all else fixed : %.1f  (%.1fσ)\n",
        dchi2_fixed, sqrt(max(dchi2_fixed, 0)))
@printf("  Δχ² after absorbing into nuisance params : %.1f  (%.1fσ)\n",
        dchi2_refit, sqrt(max(dchi2_refit, 0)))
println("  band: ", BANDLABEL)
println("  caveats: mutual occultation is $(OCC === false ? "OFF" : "ON ($OCC)"); " *
        "$(count(overlap))/$nepochs epochs have overlapping disks (grazing eclipser);")
println("           one λ per epoch; β and A partly degenerate (see the map).")
@info "figures in $outdir"

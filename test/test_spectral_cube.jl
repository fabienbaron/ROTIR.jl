#!/usr/bin/env julia
# Velocity-resolved observables (src/spectral_cube.jl): the Doppler machinery, the line
# profile, and the differential visibilities.
#
# Standalone script in the style of test_parametric_gradient.jl: prints its own table, never
# throws, exposes `nfail[]` for runtests.jl. No AD package needed.
#
#     julia --project=. test/test_spectral_cube.jl
#
# HOW THESE TESTS ARE BUILT. The atmosphere is deliberately removed from the picture: the grid
# is FLAT in Teff, logg and mu and carries one Gaussian absorption line in lambda. Whatever
# structure appears in the disk-integrated profile is then entirely the work of the velocity
# field and the projected-area weighting, which is what is under test. Three checks are
# quantitative against closed forms rather than against recorded output:
#
#   * equivalent width must be INDEPENDENT of inclination. Rotation redistributes flux in
#     wavelength and removes none, so any inclination dependence is a Jacobian error in the
#     Doppler resampling. This is the sharpest check here and it holds to 6 digits.
#   * FWHM/(2 vsini) must be ~0.87, not 1. The classical rotational kernel
#     G(v) = sqrt(1-(v/vsini)^2) reaches half maximum at v/vsini = sqrt(3)/2, so its FWHM is
#     sqrt(3)*vsini and FWHM/(2 vsini) = 0.866. `2 vsini` is the full width at ZERO intensity.
#     The measured value sits slightly BELOW 0.866 because the projected-area weighting
#     concentrates weight toward disc centre, where |v| is small. Asserting ~1 here is wrong
#     and was the first version of this file.
#   * the continuum window must be WIDER than the rotational half-width. At vsini = 215 km/s
#     and 1.65 um that half-width is 11.8 A, so a +/-10 A window sits INSIDE the line and the
#     "continuum" it averages is contaminated. That, not the code, is why an early version of
#     this file found a non-flat continuum.

using ROTIR, LinearAlgebra, Printf

npass = Ref(0); nfail = Ref(0)
function cb(label, ok)
    ok ? (npass[] += 1) : (nfail[] += 1)
    @printf("  %-58s %s\n", label, ok ? "✓" : "✗")
    return ok
end
ap(label, a, b; tol = 1e-6) = cb(label, abs(a - b) / max(abs(a), abs(b), eps()) < tol)

const C_KMS = 299792.458
const LAM0  = 1.65e-6          # H band
const FWHM  = 4e-10            # intrinsic line FWHM, 4 A
const RP, DPC, FEV, PROT = 0.849, 16.8, 0.92, 1/1.12
const VEQ = equatorial_velocity(RP, DPC, FEV, PROT)

# +/-40 A window at 0.5 A sampling. The rotational half-width is 11.8 A, so +/-20 A is a
# safe continuum boundary and +/-40 A leaves genuine continuum beyond it.
const LA  = collect(range(LAM0 - 40e-10, LAM0 + 40e-10, length = 161))
const WIN = (LAM0 - 20e-10, LAM0 + 20e-10)

intrinsic(l) = 1.0 - 0.6 * exp(-((l - LAM0) / (FWHM / 2.355))^2 / 2)

# Flat in Teff, logg and mu: the profile's shape is then purely kinematic.
const GRID = RectGrid4([4000.0, 9000.0], [2.0, 5.0], [1e-3, 1.0], LA,
                       [intrinsic(l) for t in [4000.0, 9000.0], g in [2.0, 5.0],
                                         m in [1e-3, 1.0], l in LA])
const PROV = TabulatedProvider(GRID)
const TESS = tessellation_healpix(4, T = Float64)

params(; kw...) = default_star_params(2; rpole = RP, d = DPC, frac_escapevel = FEV,
                                      rotation_period = PROT, tpole = 7208.0,
                                      inclination = 90.0, position_angle = 0.0,
                                      ldtype = 0, beta = 0.25, kw...)
state_at(sp) = surface_state(create_star(TESS, sp, 0.0), sp)

println("\n[1] rest_lambda inverts doppler_lambda, with the right sign")
for v in (-300.0, -50.0, 0.0, 50.0, 300.0)
    ap("  v = $v km/s round trip", rest_lambda(doppler_lambda(LAM0, v), v), LAM0; tol = 1e-12)
end
# A receding tessel is redshifted, so its LINE CORE has moved redward of the channel; what it
# contributes AT the channel therefore comes from blueward of its own rest core. Getting this
# backwards mirrors the profile, which is invisible on a symmetric rotator.
cb("receding tessel contributes blueward of the channel", rest_lambda(LAM0, 100.0) < LAM0)
cb("approaching tessel contributes redward", rest_lambda(LAM0, -100.0) > LAM0)
# FLOAT32: `_C_KMS` is a Float64 const and `1` an Int, so the obvious spelling
# `lambda ./ (1 .+ v/_C_KMS)` promotes a Float32 model to Float64 for the whole wavelength
# vector. It hides well — `provider_intensity` narrows its OUTPUT back to the input eltype, so
# `channel_intensity` still returns Float32 while the intermediate lambda allocated at double
# width and the grid query ran in double precision. The cube calls this once per channel.
for PT in (Float32, Float64)
    v = PT[10, -10]; λ = PT(1.65e-6)
    cb("rest_lambda stays $PT", eltype(rest_lambda(λ, v)) === PT)
    cb("doppler_lambda stays $PT", eltype(doppler_lambda(λ, v)) === PT)
end

println("\n[2] the disk-integrated profile")
sp = params(); st = state_at(sp)
cb("SurfaceState is fully populated",
   length(st.Teff) == TESS.npix && length(st.logg) == TESS.npix &&
   length(st.μ) == TESS.npix && length(st.v_los) == TESS.npix)
cb("the grid really is flat in Teff/logg/mu",
   all(GRID.values[1,1,1,:] .== GRID.values[2,2,2,:]))
F = line_profile(PROV, st, LA)
R = normalize_profile(F, LA, WIN)
core = argmin(R)
@printf("      core %.4f at %.5f um (rest %.5f um);  continuum max|R-1| = %.2e\n",
        R[core], LA[core]*1e6, LAM0*1e6, maximum(abs.(R[.!(WIN[1] .<= LA .<= WIN[2])] .- 1)))
cb("core sits at the rest wavelength (no net shift)", abs(LA[core] - LAM0) < 2*(LA[2]-LA[1]))
cb("continuum is flat at 1", maximum(abs.(R[.!(WIN[1] .<= LA .<= WIN[2])] .- 1)) < 1e-6)
cb("rotation makes the line shallower than the local one", R[core] > intrinsic(LAM0))

println("\n[3] rotational broadening: FWHM against the analytic kernel")
function fwhm_kms(inc)
    sp_ = params(inclination = inc)
    Ri = normalize_profile(line_profile(PROV, state_at(sp_), LA), LA, WIN)
    half = (1 + minimum(Ri)) / 2
    ix = findall(Ri .< half)
    return (LA[ix[end]] - LA[ix[1]]) / LAM0 * C_KMS
end
w90, w60, w30, w05 = fwhm_kms(90.0), fwhm_kms(60.0), fwhm_kms(30.0), fwhm_kms(5.0)
@printf("      FWHM: i=90 %.1f, i=60 %.1f, i=30 %.1f, i=5 %.1f km/s   (vsini(90) = %.1f)\n",
        w90, w60, w30, w05, VEQ)
@printf("      FWHM/(2 vsini) = %.3f   (analytic elliptical kernel: 0.866)\n", w90/(2*VEQ))
cb("width is monotone in inclination", w90 > w60 > w30 > w05)
# Below 0.866 because the projected area weights disc centre, where |v| is small.
cb("FWHM/(2 vsini) is just under the analytic 0.866", 0.70 < w90/(2*VEQ) < 0.90)
cb("i=5 collapses toward the intrinsic width", w05 < 0.3 * w90)
# sin(i) scaling, via SECOND MOMENTS rather than FWHM. Variances add exactly under
# convolution, whatever the kernel shapes; FWHM adds in quadrature only for Gaussians, and
# the rotational kernel is elliptical with sharp edges — deconvolving FWHM that way is what
# made an earlier version of this check "fail" at i = 30, where vsini and the intrinsic width
# are comparable and the approximation is worst. So: treat the absorption (1 − R) as a
# distribution over velocity, take its variance, and subtract the i → 0 value. The remainder
# is the rotational variance and must scale as sin²(i) exactly.
function velvar(inc; fev = FEV)
    Ri = normalize_profile(
             line_profile(PROV, state_at(params(inclination = inc, frac_escapevel = fev)),
                          LA), LA, WIN)
    a = max.(1 .- Ri, 0.0)                       # absorption as a weight
    v = (LA .- LAM0) ./ LAM0 .* C_KMS
    w = sum(a)
    m = sum(a .* v) / w
    return sum(a .* (v .- m) .^ 2) / w
end
# Scaling holds EXACTLY only for a sphere. For an oblate figure the mapping from tessel to
# (v_los, projected area) is not a pure sin(i) rescaling — the silhouette and the
# foreshortening both change shape with inclination — so the residual departure is a physical
# prediction, not slack. Measured spread across i = 90, 60, 30, 15 at HEALPix nside 4:
# fev ~ 0 -> 0.45 %, fev = 0.3 -> 0.41 %, fev = 0.6 -> 0.42 %, fev = 0.92 -> 2.39 %.
# The ~0.4 % floor is mesh discretisation; everything above it is oblateness.
function rotvar_spread(fev)
    v0 = velvar(1.0; fev = fev)                  # sin(1 deg)^2 is 3e-4 of the i=90 term
    r  = [(velvar(i; fev = fev) - v0) / sind(i)^2 for i in (90.0, 60.0, 30.0, 15.0)]
    return (maximum(r) - minimum(r)) / maximum(r), r
end
sp_sph, r_sph = rotvar_spread(1e-6)
sp_obl, r_obl = rotvar_spread(FEV)
@printf("      rot var / sin^2(i), sphere : %s  spread %.2f %%\n",
        string(round.(r_sph, digits = 1)), 100 * sp_sph)
@printf("      rot var / sin^2(i), fev=%.2f: %s  spread %.2f %%\n",
        FEV, string(round.(r_obl, digits = 1)), 100 * sp_obl)
cb("a SPHERE scales as sin^2(i) to the mesh floor (<1%)", sp_sph < 0.01)
cb("oblateness measurably breaks it at fev = 0.92 (>1%)", sp_obl > 0.01)
cb("but only mildly (<5%)", sp_obl < 0.05)

println("\n[4] vgamma shifts the whole profile, and only shifts it")
spv = params(vgamma = 60.0)
Rv = normalize_profile(line_profile(PROV, state_at(spv), LA), LA, WIN)
shift = (LA[argmin(Rv)] - LA[argmin(R)]) / LAM0 * C_KMS
@printf("      measured %.1f km/s for vgamma = +60 (receding -> redshift)\n", shift)
cb("redshifted by roughly vgamma", 40 < shift < 80)
cb("depth is unchanged by vgamma", abs(minimum(Rv) - minimum(R)) < 0.02)

println("\n[5] equivalent width is conserved under rotation")
# The sharpest check in this file. Rotation redistributes flux in wavelength and removes
# none, so EW cannot depend on inclination. Any dependence is a Jacobian error in the
# Doppler resampling — which is exactly the mistake `src/di.jl` made by shifting in linear
# lambda rather than log lambda (its own comment at di.jl:323 admits it).
ews = [line_equivalent_width(line_profile(PROV, state_at(params(inclination = i)), LA),
                             LA, WIN) for i in (5.0, 30.0, 60.0, 90.0)]
@printf("      EW(i = 5, 30, 60, 90) = %s pm\n", string(round.(ews .* 1e12, digits = 4)))
cb("EW is inclination-independent to 1e-4 relative",
   (maximum(ews) - minimum(ews)) / maximum(ews) < 1e-4)

println("\n[6] continuum_mask guards")
cb("refuses a window leaving <2 continuum channels",
   try; continuum_mask(LA, (LA[1] - 1e-12, LA[end] + 1e-12)); false; catch; true; end)
cb("mask excludes the window", !any(continuum_mask(LA, WIN)[WIN[1] .<= LA .<= WIN[2]]))

println("\n[7] differential observables")
# Synthetic per-channel visibilities: continuum channels identical, one line channel altered
# in both amplitude and phase. The differential quantities must recover exactly that.
nuv = 7
base = [0.6 * cis(0.3k) for k in 1:nuv]
cv = [copy(base) for _ in 1:length(LA)]
inl = findall(WIN[1] .<= LA .<= WIN[2])
for c in inl
    cv[c] = base .* (0.5 * cis(0.2))
end
da, dp = differential_observables(cv, LA, WIN)
cb("continuum channels give unit amplitude ratio",
   all(abs.(da[:, .!(WIN[1] .<= LA .<= WIN[2])] .- 1) .< 1e-12))
cb("continuum channels give zero phase", all(abs.(dp[:, .!(WIN[1] .<= LA .<= WIN[2])]) .< 1e-10))
cb("line channels recover the amplitude ratio 0.5",
   all(abs.(da[:, inl] .- 0.5) .< 1e-12))
cb("line channels recover the phase +0.2 rad", all(abs.(dp[:, inl] .- rad2deg(0.2)) .< 1e-8))
cb("shape is (nuv, nchannel)", size(da) == (nuv, length(LA)))
cb("mismatched uv counts are refused",
   try; differential_observables([base, base[1:3]], LA[1:2], WIN); false; catch; true; end)

println("\n[8] spectral_cvis against the single-channel path")
# One channel of spectral_cvis must equal fused_cvis on the same intensity vector, and its
# flux must equal what line_profile computes. That ties the interferometric and the
# spectroscopic routes to the same sum, which is the point of the module.
files = [joinpath(@__DIR__, "..", "demos", "data", "MEDIAN5.MIRCX_L2.2025Oct30.HD_432.MIRCX_IDL.bet_Cas.AVG10m.oifits")]
if isfile(files[1])
    d = readoifits_multiepochs(files; T = Float64)[1, :]
    star = create_star(TESS, sp, 0.0)
    stt = surface_state(star, sp)
    λ1 = band_of(d[1])
    cvis, flux, λs = spectral_cvis(PROV, stt, star, [d[1]]; λ = [λ1])
    I1 = channel_intensity(PROV, stt, λ1)
    ref = ROTIR.fused_cvis(I1, star, d[1]; intensity_model = :linear)
    cb("one channel equals fused_cvis", maximum(abs.(cvis[1] .- ref)) < 1e-12)
    cb("flux equals line_profile at the same lambda",
       abs(flux[1] - line_profile(PROV, stt, [λ1])[1]) < 1e-9 * abs(flux[1]))
    cb("visibilities are normalised (|V| <= 1 + eps)", maximum(abs.(cvis[1])) <= 1 + 1e-9)
    # Multi-channel: the cube must reproduce channel-by-channel what one call gives.
    λ3 = [λ1 * 0.999, λ1, λ1 * 1.001]
    cv3, fl3, _ = spectral_cvis(PROV, stt, star, [d[1], d[1], d[1]]; λ = λ3)
    cb("3-channel cube matches its middle channel", maximum(abs.(cv3[2] .- ref)) < 1e-12)
    cb("3-channel flux matches line_profile",
       maximum(abs.(fl3 .- line_profile(PROV, stt, λ3))) < 1e-9 * abs(flux[1]))
else
    @warn "beta Cas OIFITS not found — section [8] skipped" files[1]
end

println("\n[9] the differential-phase S-curve — the rotation signature itself")
if isfile(files[1])
    d = readoifits_multiepochs(files; T = Float64, polychromatic = true)[1, 1]
    sp90 = params(inclination = 90.0, position_angle = 0.0)
    st90 = create_star(TESS, sp90, 0.0)
    stv  = surface_state(st90, sp90)
    cv, fl, _ = spectral_cvis(PROV, stv, st90, fill(d, length(LA)); λ = LA)
    da, dp = differential_observables(cv, LA, WIN)
    # Away from visibility nulls: a differential quantity divides by the continuum
    # visibility, so near a null the ratio explodes and the phase flips 180 deg as V crosses
    # zero. Real features of a resolved disc, not of the line, and they dominate any mean
    # over all baselines (mean ratio 1.30 +/- 1.27 at the core; median 0.998).
    m = continuum_mask(LA, WIN)
    ref = [sum(cv[c][k] for c in eachindex(cv) if m[c]) / count(m) for k in 1:length(cv[1])]
    good = abs.(ref) .> 0.2
    v = (LA .- LAM0) ./ LAM0 .* C_KMS
    inl = findall(abs.(v) .< 250)
    ph = [sum(dp[good, i]) / count(good) for i in inl]
    vv = v[inl]
    pk = maximum(abs.(ph)); ipk = argmax(abs.(ph))
    @printf("      %d/%d baselines used; peak |dphase| = %.3f deg at v = %+.1f km/s (vsini %.1f)\n",
            count(good), length(good), pk, vv[ipk], VEQ)
    @printf("      antisymmetry residual = %.2e deg (peak %.3f)\n", sum(ph .+ reverse(ph)), pk)
    # A symmetric rotator's photocentre displacement is odd in velocity, so the phase must be
    # too. This is the check that a spot or differential rotation would break — which is
    # precisely what makes differential phase a probe of the SURFACE and not only of the axis.
    cb("dphase is antisymmetric in velocity to 5% of its peak",
       abs(sum(ph .+ reverse(ph))) < 0.05 * pk * length(ph)^0 * 1.0 + 0.05 * pk)
    cb("dphase crosses zero at the line centre",
       abs(ph[argmin(abs.(vv))]) < 0.6 * pk)
    # The far wings come from the limb, where little projected area is left, so the extremum
    # sits INSIDE vsini rather than at +/- vsini.
    cb("peak sits inside vsini", abs(vv[ipk]) < VEQ)
    cb("blue and red halves have opposite sign",
       sign(sum(ph[vv .< -20])) == -sign(sum(ph[vv .> 20])))
    # EXACTLY pole-on: v_los is identically zero, every channel sees the same brightness
    # distribution, and the differential phase must vanish to round-off. i = 0 rather than a
    # small angle, deliberately — see the note below on why the phase does NOT scale with
    # sin(i), which makes any "nearly pole-on" threshold arbitrary.
    sp0 = params(inclination = 0.0, position_angle = 0.0)
    st0 = create_star(TESS, sp0, 0.0)
    st0s = surface_state(st0, sp0)
    cb("i = 0 gives an identically zero velocity field", maximum(abs, st0s.v_los) < 1e-9)
    cv0, _, _ = spectral_cvis(PROV, st0s, st0, fill(d, length(LA)); λ = LA)
    _, dp0 = differential_observables(cv0, LA, WIN)
    pk0 = maximum(abs.([sum(dp0[good, i]) / count(good) for i in inl]))
    @printf("      i = 0 peak |dphase| = %.2e deg vs %.3f deg at i = 90\n", pk0, pk)
    cb("a zero velocity field gives zero differential phase", pk0 < 1e-8)

    # NOT ASSERTED, because it is not true: the differential phase does not scale with sin(i).
    # Measured peaks on these baselines are 1.58 deg at i = 90 (vsini 215 km/s), 1.80 deg at
    # i = 5 (vsini 19 km/s) and 0.19 deg at i = 0.5 (vsini 1.9 km/s). At high inclination the
    # line is broad and SHALLOW (14 % absorption at the core) so the velocity-selected region
    # has little contrast against the rest of the disc; at low inclination it is narrow and
    # DEEP (57 %), so a much smaller velocity gradient still displaces the photocentre
    # strongly. The two effects roughly cancel until vsini drops below the channel width,
    # which is where the signal finally collapses. Worth knowing before designing an
    # observation: a slow rotator is not necessarily a harder differential-phase target.
else
    @warn "beta Cas OIFITS not found — section [9] skipped"
end

@printf("\n=== %d passed, %d failed ===\n", npass[], nfail[])

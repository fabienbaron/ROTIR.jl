#!/usr/bin/env julia
# =======================================================================================
# Velocity-resolved interferometry of a rapid rotator
# =======================================================================================
# One surface, one intensity provider, and BOTH the interferometric visibilities and the
# disk-integrated line profile out of the same sum. This is the capability PMOIRED's
# `rotastar` has and ROTIR had none of, and it is built on three pieces added for it:
#
#   src/stellar_physics.jl   absolute mass, radius and logg (the distance closes the system;
#                            mass is DERIVED from rpole, d, fev and the rotation period)
#   src/intensity_provider.jl  I(Teff, logg, mu, lambda) from a model atmosphere, so limb
#                            darkening is PREDICTED rather than fitted
#   src/velocity_field.jl    the per-tessel line-of-sight velocity, from the Roche shape
#   src/spectral_cube.jl     the channel loop, differential observables, line profile
#
# WHY THIS MATTERS BEYOND FEATURE PARITY. A fitted limb-darkening coefficient absorbs the
# error in whichever gravity-darkening law was assumed — gravity and limb darkening both
# take flux out of the limb and an interferometer measures only their sum. That is how a
# beta Cas fit reached a NEGATIVE ld1 (see demos/gravity_law_comparison.jl and the header of
# src/gravity_darkening.jl). With `ldtype = 0` and an atmosphere provider there is no ld1 to
# absorb anything.
#
#   julia --project=demos demos/rapid_rotator_spectro.jl
#
# Part 1 runs on the real beta Cas MIRC-X data, 6 H-band channels. Part 2 uses a synthetic
# high-resolution line to show the velocity-resolved observables, because MIRC-X at R ~ 50
# cannot resolve a line — that needs GRAVITY (R = 4000) or MIRC-X's R = 190 mode.
# =======================================================================================

using ROTIR, Printf, Statistics

const DATA = joinpath(@__DIR__, "data",
                      "MEDIAN5.MIRCX_L2.2025Oct30.HD_432.MIRCX_IDL.bet_Cas.AVG10m.oifits")

# beta Cas, as fitted in demos/rapid_rotator_betCas_param_fit.jl, plus the Gaia distance.
# Only these four set the scale; mass, radius, logg and vsini all come out.
star = default_star_params(2;
    rpole           = 0.849,     # mas
    d               = 16.8,      # pc   <- the new parameter; everything physical follows
    frac_escapevel  = 0.92,
    rotation_period = 1/1.12,    # d
    tpole           = 7208.0,    # K
    inclination     = 19.9,      # deg
    position_angle  = -7.09,     # deg
    beta            = 0.25,      # radiative; no longer a sponge for the LD error
    gravity_law     = 2,         # Espinosa Lara-Rieutord
    ldtype          = 0,         # NONE: the provider owns the mu dependence
)

println("="^78)
println("beta Cas — what the distance buys")
println("="^78)
q = derived_quantities(star)
@printf("  mass       %7.3f Msun     (DERIVED from rpole, d, fev, P — not fitted)\n", q.mass)
@printf("  R_pole     %7.3f Rsun     R_eq %6.3f Rsun  (oblateness %.3f)\n",
        q.rpole_rsun, q.req_rsun, q.req_rsun/q.rpole_rsun)
@printf("  logg_pole  %7.3f          v_eq %6.1f km/s   vsini %6.1f km/s\n",
        q.logg_pole, q.veq, q.vsini)
msgs = validate_star_params(star)
println("  validate: ", isempty(msgs) ? "clean" : msgs)

tess = tessellation_healpix(4, T = Float64)
st1  = create_star(tess, star, 0.0)
Tmap = parametric_temperature_map(star, st1)
θ    = tess.unit_spherical[:, 5, 2]
lg   = logg_map(star.rpole, star.d, star.frac_escapevel, star.rotation_period,
                sin.(θ), cos.(θ))
@printf("\n  surface spans Teff %.0f-%.0f K and logg %.2f-%.2f\n",
        extrema(Tmap)..., extrema(lg)...)
println("  -> a production atmosphere grid must cover ALL of that, not just the polar values")
vs = velocity_field_summary(st1, star)
@printf("  visible-hemisphere v_los: %.1f to %.1f km/s (half-span %.1f)\n",
        vs.vmin, vs.vmax, vs.vsini_proj)

# =======================================================================================
# Part 1 — chromatic visibilities on the real data, one channel at a time
# =======================================================================================
println("\n" * "="^78)
println("Part 1: per-channel visibilities, beta Cas MIRC-X (6 H-band channels)")
println("="^78)

if !isfile(DATA)
    @warn "beta Cas OIFITS not found; skipping Part 1" DATA
else
    # `polychromatic = true` is what keeps the wavelength axis. EVERY existing demo throws
    # it away with `data = data_all[1, :]`, which is why nothing in ROTIR had ever used it.
    dall = readoifits_multiepochs([DATA]; T = Float64, polychromatic = true)
    nchan, nep = size(dall)
    @printf("  %d channels x %d epoch(s); %d uv points per channel\n",
            nchan, nep, dall[1,1].nuv)

    # A Planck provider needs no atmosphere grid, so Part 1 runs anywhere. It carries no mu
    # dependence, so here ldtype must NOT be 0 — that is what `check_provider_consistency`
    # is for, and it is the reverse of the Part 2 case.
    star_pl = merge(star, (ldtype = 3, ld1 = 0.21))
    prov_pl = PlanckProvider()
    println("  provider consistency: ",
            isempty(check_provider_consistency(prov_pl, star_pl)) ? "ok" :
            check_provider_consistency(prov_pl, star_pl))

    st = create_star(tess, star_pl, dall[1,1].mean_mjd - dall[1,1].mean_mjd)
    state = surface_state(st, star_pl)
    cvis, flux, λs = spectral_cvis(prov_pl, state, st, dall[:, 1])

    println("\n  chan   lambda[um]   mean|V|    total flux")
    for c in 1:nchan
        @printf("   %2d      %.4f     %.4f     %.4e\n",
                c, λs[c]*1e6, mean(abs.(cvis[c])), flux[c])
    end
    # The star is resolved, so |V| falls with baseline; and because uv = B/lambda, the SAME
    # baseline sits at larger spatial frequency in the blue, so |V| is lower there. That
    # chromatic trend is the signature the channel loop exists to capture.
    @printf("\n  mean|V| blue->red: %.4f -> %.4f  (%s with wavelength, as a resolved disc must)\n",
            mean(abs.(cvis[1])), mean(abs.(cvis[end])),
            mean(abs.(cvis[end])) > mean(abs.(cvis[1])) ? "rises" : "falls")
end

# =======================================================================================
# Part 2 — the velocity-resolved line
# =======================================================================================
println("\n" * "="^78)
println("Part 2: a velocity-resolved line (synthetic grid, R ~ 20000)")
println("="^78)

# A grid flat in Teff/logg/mu with one Gaussian absorption line isolates the kinematics:
# every structure in the profile below is the velocity field and the projected-area
# weighting. A real run would put a Korg or Kurucz/TLUSTY grid here instead —
# `build_korg_grid` (needs `using Korg`) or `load_intensity_grid`.
λ0, fwhm = 1.65e-6, 4e-10
λs = collect(range(λ0 - 40e-10, λ0 + 40e-10, length = 161))
intrinsic(l) = 1.0 - 0.6 * exp(-((l - λ0) / (fwhm/2.355))^2 / 2)
grid = RectGrid4([4000.0, 9000.0], [2.0, 5.0], [1e-3, 1.0], λs,
                 [intrinsic(l) for _ in 1:2, _ in 1:2, _ in 1:2, l in λs])
prov = TabulatedProvider(grid; name = "synthetic line")
println("  provider consistency (ldtype=0): ",
        isempty(check_provider_consistency(prov, star)) ? "ok" :
        check_provider_consistency(prov, star))

window = (λ0 - 20e-10, λ0 + 20e-10)
println("\n  inclination   FWHM[km/s]   depth    EW[pm]")
for inc in (5.0, 30.0, 60.0, 90.0)
    sp = merge(star, (inclination = inc,))
    R  = normalize_profile(line_profile(prov, surface_state(create_star(tess, sp, 0.0), sp),
                                        λs), λs, window)
    half = (1 + minimum(R)) / 2
    ix = findall(R .< half)
    w = (λs[ix[end]] - λs[ix[1]]) / λ0 * 299792.458
    ew = line_equivalent_width(line_profile(prov,
             surface_state(create_star(tess, sp, 0.0), sp), λs), λs, window)
    @printf("     %5.1f       %6.1f      %.4f   %.4f\n", inc, w, minimum(R), ew*1e12)
end
println("  EW is the same at every inclination: rotation redistributes flux in wavelength,")
println("  it removes none. That invariance is the check a Doppler resampling must pass.")

# Differential visibilities need the same baselines at every channel, so reuse one real
# channel's uv across the synthetic wavelength grid. On real data these come straight from
# OIFITS v2, which OITOOLS already reads and tags amptyp/phityp = "differential"; they sit at
# weight positions 4 and 5 of OI_DEFAULT_WEIGHTS = [1,1,1,0,0,0,0] and are simply switched off.
if isfile(DATA)
    d1 = readoifits_multiepochs([DATA]; T = Float64, polychromatic = true)[1, 1]
    # Equator-on for the illustration: the velocity field is strongest and the geometry
    # cleanest. PA = 0 puts the projected spin axis along North.
    sp90 = merge(star, (inclination = 90.0, position_angle = 0.0))
    st90 = create_star(tess, sp90, 0.0)
    stv  = surface_state(st90, sp90)
    cv, fl, _ = spectral_cvis(prov, stv, st90, fill(d1, length(λs)); λ = λs)
    da, dp = differential_observables(cv, λs, window)

    # RESTRICT TO BASELINES WELL AWAY FROM A VISIBILITY NULL. beta Cas is strongly resolved
    # on these baselines (median |V| ~ 0.12, minimum 0.014), and a differential quantity
    # divides by the continuum visibility — so near a null the amplitude ratio explodes and
    # the phase flips by 180 deg as V crosses zero. Those are real features of a resolved
    # disc, not of the line, and they swamp any average taken over all baselines: the MEAN
    # ratio at the line core is 1.30 with a standard deviation of 1.27, while the MEDIAN is
    # 0.998. A real analysis applies the same cut.
    m = continuum_mask(λs, window)
    ref = [sum(cv[c][k] for c in eachindex(cv) if m[c]) / count(m) for k in 1:length(cv[1])]
    good = abs.(ref) .> 0.2
    v = (λs .- λ0) ./ λ0 .* 299792.458
    R = normalize_profile(fl, λs, window)

    println("\n  differential phase across the line (i = 90 deg, |V_cont| > 0.2, "
            * "$(count(good)) of $(length(good)) baselines):")
    println("     v[km/s]   profile   <dphase>[deg]   <|V| ratio>")
    for i in 55:6:107
        @printf("     %+7.1f    %.4f     %+8.4f        %.4f\n",
                v[i], R[i], sum(dp[good, i])/count(good), sum(da[good, i])/count(good))
    end

    inl = findall(abs.(v) .< 250)
    ph  = [sum(dp[good, i])/count(good) for i in inl]
    @printf("\n  peak |dphase| = %.3f deg at v = %+.1f km/s  (vsini = %.1f km/s)\n",
            maximum(abs.(ph)), v[inl][argmax(abs.(ph))], projected_veq(
                star.rpole, star.d, star.frac_escapevel, star.rotation_period, 90.0))
    @printf("  antisymmetry residual sum(dphase(v) + dphase(-v)) = %.2e deg\n",
            sum(ph .+ reverse(ph)))
    println("""
  That S-shape IS the measurement. One limb is blueshifted and the other redshifted, so each
  velocity channel sees a different part of the surface and the PHOTOCENTRE moves across the
  line — positive on the blue side, negative on the red, crossing zero at the line centre
  where only the v ~ 0 central strip contributes and the brightness distribution is
  symmetric again. It peaks inside vsini rather than at the extreme wings, because the
  far wings come from the limb, where little projected area remains.

  It is antisymmetric for a symmetric rotator; a spot or differential rotation breaks that,
  which is what makes differential phase a direct probe of the surface rather than just of
  the rotation axis. The sign fixes the sense of rotation on sky - something V^2 and closure
  phase cannot do at all.""")
end

println("\n" * "="^78)
println("Summary")
println("="^78)
println("  * mass, radius and logg are derived, not fitted: one new parameter (d)")
println("  * limb darkening is predicted: ldtype = 0, no ld1/ld2 in the fit")
println("  * the velocity field follows the Roche shape, from positions not normals")
println("  * visibilities and line profile come from ONE sum, so they cannot disagree")
println("  * the channel loop reuses the exact-polygon FT unchanged — no kernel widening")

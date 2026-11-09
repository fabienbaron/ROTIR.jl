# Spica (alpha Vir) — tessellated binary model at the joint OI + RV best fit
#
#   julia -t 8 --project=demos demos/spica_tessellated_model.jl
#
# Roche surfaces with von Zeipel gravity darkening and mutual irradiation, compared epoch by
# epoch with 2007_2012_2015.Spica.oifits; then a two-map image reconstruction and per-epoch
# temperature maps. The orbit, radii and mass ratio come from
# demos/spica_tessellated_joint_fit.jl (results/spica_tessellated_joint_best.txt).

using ROTIR
using OITOOLS: split_by_time
using DelimitedFiles, PythonPlot, PythonCall
import Statistics: mean

# The merged OIFITS file contains epochs from 2007, 2012, and 2015 campaigns
oifitsfile = joinpath(@__DIR__, "data", "2007_2012_2015.Spica.oifits")
data_all = readoifits(oifitsfile)[1,1]

# Split into epochs wherever consecutive observations are more than half a day apart.
# epoch_mean_mjd is the mean time of each epoch's points (all observables).
data, epoch_mean_mjd = split_by_time(data_all; gap = 0.5)
nepochs = length(data)
tepochs = epoch_mean_mjd .+ 2400000.5
for i in 1:nepochs
    println("Epoch $i: MJD=$(round(epoch_mean_mjd[i], digits=4)), nV2=$(data[i].nv2), nT3=$(data[i].nt3phi)")
end
bands = band_of.(data)    # mean wavelength of each epoch (m), for the Planck intensity

# set up parameters as constants for both stars
include(joinpath(@__DIR__, "spica_params.jl"))
# Orbit and radii from a joint fit of this tessellated model (below) to this OIFITS file and the
# radial velocities in demos/data/spica_rv.csv (Keplerian RVs; K1, K2 and gamma solved by linear
# least squares, q = K1/K2). P is the anomalistic period (periastron to periastron), fixed at
# 4.01454 d (Shobbrook et al. 1972, as in Aufdenberg et al. 2015); apsidal motion is free and
# comes out at U = 136 yr (literature 128-139 yr). K1 = 123.6, K2 = 193.5, gamma = -0.1 km/s,
# parallax 12.07 mas, M1 = 11.0, M2 = 7.05 Msun.
P_FIT = 4.01454
I_FIT, OMEGA_FIT, OMEGAP_FIT = 116.55023101969826, 310.20744924267166, 73.38916471090769
A_FIT, E_FIT, T0_FIT, DOMEGA_FIT = 1.5657838589894628, 0.12617824588937568, 2.454191182745739e6, 0.007262223743318985
Q_FIT = 0.6384678295749032   # M2/M1 = K1/K2
# Polar radii (mas) from the same fit. Limb darkening is the H-band power-law exponent
# predicted by Kurucz ATLAS9 intensities (demos/data/spica_H_kurucz.fits) at each star's
# Teff and logg; the data alone do not constrain it. Rotation rates follow the observed
# v sin i = 161 / 70 km/s at the fitted radii (spica_params.jl derives F1/F2 from 7.40/3.74 Rsun).
RPOLE1_FIT, RPOLE2_FIT = 0.40609325547403324, 0.20087736894609348
LD1_KURUCZ, LD2_KURUCZ = 0.077, 0.068    # fixed, not fitted
F1_FIT, F2_FIT = F1_SYNC * RPOLE1 / RPOLE1_FIT, F2_SYNC * RPOLE2 / RPOLE2_FIT
orbit_fit = (inclination = 180.0 - I_FIT, position_angle = OMEGA_FIT - 180.0,
             i = I_FIT, Ω = OMEGA_FIT, ω = OMEGAP_FIT, a = A_FIT, e = E_FIT, T0 = T0_FIT,
             P = P_FIT, dω = DOMEGA_FIT)
p1 = merge(SPICA_P1, orbit_fit, (rpole = RPOLE1_FIT, q = Q_FIT, ld1 = LD1_KURUCZ, rotation_period = P_FIT / F1_FIT))
p2 = merge(SPICA_P2, orbit_fit, (rpole = RPOLE2_FIT, q = 1 / Q_FIT, ld1 = LD2_KURUCZ, rotation_period = P_FIT / F2_FIT))
s1p = starparameters(RPOLE1_FIT, TPOLE1, 0.0, LDTYPE, LD1_KURUCZ, 0.0, BETA, 0.0,
                     180.0 - I_FIT, OMEGA_FIT - 180.0, 0.0, P_FIT / F1_FIT)
s2p = starparameters(RPOLE2_FIT, TPOLE2, 0.0, LDTYPE, LD2_KURUCZ, 0.0, BETA, 0.0,
                     180.0 - I_FIT, OMEGA_FIT - 180.0, 0.0, P_FIT / F2_FIT)
bp = binaryparameters(s1p, s2p, D_PC, I_FIT, OMEGA_FIT, OMEGAP_FIT, P_FIT, A_FIT, E_FIT,
                      T0_FIT, Q_FIT, [1.0, 1.0], 0.0, DOMEGA_FIT)
t1, t2 = tessellation_healpix(4), tessellation_healpix(3)

# for alternating params:
#p1 = roche_params(RPOLE1, TPOLE1, Q_BIN,     beta, a_mas)
#p2 = roche_params(RPOLE2, TPOLE2, 1 / Q_BIN, beta, a_mas)   # q is M_companion/M_self
#s1p = starparameters(RPOLE1, TPOLE1, 0.0, 3, 0.30, 0.0, beta, 0.0, INC, OMEGA, 0.0, PORB)
#s2p = starparameters(RPOLE2, TPOLE2, 0.0, 3, 0.30, 0.0, beta, 0.0, INC, OMEGA, 0.0, PORB)
#bp = binaryparameters(s1p, s2p, DIST, INC, OMEGA, 90.0, PORB, a_mas, 0.0, 0.0, QBIN, [-1.0, -1.0], 0.0, 0.0)

function build(p1, p2, t1::tessellation, t2::tessellation, ts; vcons=true, npoints=9)
    o1t = vcons ? roche_omega_table(t1, p1, bp; secondary=false, npoints=npoints) : nothing
    o2t = vcons ? roche_omega_table(t2, p2, bp; secondary=true,  npoints=npoints) : nothing
    n = length(ts)
    # Concretely typed (as create_star_multiepochs does): a Vector{Any} here makes the
    # threaded loop in setup_oi! type-unstable.
    s1 = Array{stellar_geometry}(undef, n); s2 = Array{stellar_geometry}(undef, n)
    for i in 1:n
        D = compute_separation(bp, ts[i])
        s1[i], s2[i] = create_binary_geometry(t1, p1, t2, p2, bp, ts[i];
                                              omega1 = o1t === nothing ? nothing : o1t(D),
                                              omega2 = o2t === nothing ? nothing : o2t(D))
    end
    return s1, s2
end
s1, s2 = build(p1, p2, t1, t2, tepochs)

# The secondary at each uv point's own observation time: it moves 0.24-0.52 mas during a
# 2015 night, which one position per epoch would smear out of the closure phases.
phases = [binary_phase_shift(data[i], bp) for i in 1:nepochs]

# Epochs where the stars overlap on the sky (flagged with a warning).
overlap = check_binary_overlap(s1, s2, bp, tepochs)
# Mutual occultation is off for now; to model it at the overlapping epochs, use
# occ = [o ? :exact : false for o in overlap]
occ = fill(false, nepochs)

# von Zeipel maps at each epoch's JD
tm1 = [temperature_map_vonZeipel_roche_single(p1, s1[i], tepochs[i]) for i in 1:nepochs]
tm2 = [temperature_map_vonZeipel_roche_single(p2, s2[i], tepochs[i]; secondary=true)
       for i in 1:nepochs]
# mutual irradiation
kern = [reflection_kernels(s1[i], s2[i]) for i in 1:nepochs]

function heat(i, A; t1=tm1, t2=tm2, k=kern)
    A == 0 && return t1[i], t2[i]
    handle_reflection(s1[i], t1[i], s2[i], t2[i];
                      albedo1=A, albedo2=A, kernels=k[i])
end

function chi2_split(A; imodel=:planck)
    cv2 = ca = cp = 0.0
    for i in 1:nepochs
        h1, h2 = heat(i, A)
        v2, t3a, t3p = binary_observables(h1, s1[i], h2, s2[i], data[i], phases[i];
                                          intensity_model=imodel, band=bands[i],
                                          occultation=occ[i])
        cv2 += sum(abs2, (v2 .- data[i].v2) ./ data[i].v2_err)
        ca  += sum(abs2, (t3a .- data[i].t3amp) ./ data[i].t3amp_err)
        cp  += sum(abs2, mod360(t3p .- data[i].t3phi) ./ data[i].t3phi_err)
        print("Chi2 = ", cv2," ", ca," ", cp, "\n")
        ndata_epoch = data[i].nv2 + data[i].nt3amp + data[i].nt3phi
    end
    return cv2, ca, cp
end

cv2, ca, cp = chi2_split(ALBEDO_DEFAULT)




################### COMPUTE CHI2 PER EPOCH ###################
total_chi2 = 0.0
model_obs = Vector{NamedTuple}(undef, nepochs)

for i in 1:nepochs
    # Use precise mean MJD from the data (not the rounded filter value)
    tepoch = epoch_mean_mjd[i] + 2400000.5
    # the geometries, irradiated maps and phase ramps built above
    star1, star2 = s1[i], s2[i]
    m1, m2 = heat(i, ALBEDO_DEFAULT)
    offset_x, offset_y = orbit_to_rotir_offset(bp, tepoch)    # for the separation printed below
    phase = phases[i]

    v2_model, t3amp_model, t3phi_model = binary_observables(m1, star1, m2, star2, data[i], phase;
                                                            intensity_model=:planck, band=bands[i],
                                                            occultation=occ[i])
    model_obs[i] = (v2=v2_model, t3amp=t3amp_model, t3phi=t3phi_model, visamp=Float64[], visphi=Float64[])

    chi2_epoch = binary_chi2_f(m1, star1, m2, star2, data[i], phase;
                               intensity_model=:planck, band=bands[i], occultation=occ[i])
    global total_chi2 += chi2_epoch

    sep = sqrt(offset_x^2 + offset_y^2)
    ndata_epoch = data[i].nv2 + data[i].nt3amp + data[i].nt3phi
    println("Epoch $i (MJD $(round(epoch_mean_mjd[i], digits=4))): sep=$(round(sep,digits=3)) mas, chi2r=$(round(chi2_epoch/ndata_epoch, digits=2))")
end
ndata_total = sum(d.nv2 + d.nt3amp + d.nt3phi for d in data)
println("\nTotal chi2 = $(round(total_chi2, digits=1)), reduced = $(round(total_chi2/ndata_total, digits=2))")


############ IMAGE RECONSTRUCTION ############
# One geometry pair and one phase ramp PER EPOCH: s1, s2 from build() and phases from the
# orbit (both above), so length(phases) == length(data) == nepochs.
setup_oi!(data, s1)
setup_oi!(data, s2)
# starting maps: the von Zeipel temperatures (the same tessellation at every epoch)
x = vcat(tm1[1], tm2[1])
# Quadratic gradient penalty ("sobel2") for smooth temperature maps, on all pixels.
# At 1e3: chi2r 3.65, primary 20-30 kK, secondary 18-24 kK. Lower weights fit noise
# (at 1e-5 most visible pixels go to 0 K and a few to ~5e6 K); 1e4 and above barely move
# from the von Zeipel start (chi2r 3.95-4.05). "sobel" instead keeps sharp edges.
reg1 = [["sobel2", 1e3, sobel_gradient_healpix(4), 1:length(tm1[1])]]
reg2 = [["sobel2", 1e3, sobel_gradient_healpix(3), 1:length(tm2[1])]]
tmap = binary_reconstruct_oi(x, data, s1, s2, phases; maxiter=500, regularizers1=reg1, regularizers2=reg2)
r1, r2 = copy.(split_binary_map(tmap, s1))
plot2d_binary(r1, r2, s1[1], s2[1], bp, tepochs[1], colormap="RdBu", rotation_axis=true,
        rotation_arrow=true, figtitle="Reconstruction, JD $(tepochs[1]) (Epoch 1)")


############ PLOT PER EPOCH WITH ALL EFFECTS ############
function timevar(t1, p1, t2, p2, bp, tepoch; irradiate = true)    
    star1, star2 = create_binary_geometry(t1, p1, t2, p2, bp, tepoch)
    # von Zeipel on the Roche surface: T = tpole*(g/g_pole)^beta, computed from the actual
    # distorted geometry, so beta = 0 gives a uniform map on a still-distorted star.
    m1 = temperature_map_vonZeipel_roche_single(p1, star1, tepoch)
    m2 = temperature_map_vonZeipel_roche_single(p2, star2, tepoch; secondary = true)
    if irradiate
        m1, m2 = handle_reflection(star1, m1, star2, m2;
                                   albedo1 = ALBEDO_DEFAULT, albedo2 = ALBEDO_DEFAULT)
    end
    return star1, star2, m1, m2, bp
end
for i in 1:nepochs
    tepoch = epoch_mean_mjd[i] + 2400000.5
    star1, star2, m1, m2 = timevar(t1, p1, t2, p2, bp, tepoch)   # bp is unchanged: not reassigned
    plot2d_binary(m1, m2, star1, star2, bp, tepoch, colormap="RdBu", rotation_axis=true,
        rotation_arrow=true, figtitle="JD $tepoch (Epoch $i)")
    mkpath(joinpath(@__DIR__, "results", "spica_tess"))
    savefig(joinpath(@__DIR__, "results", "spica_tess", "spica_prox_epoch$i.png"))
end


############# RASTERIZATION ##############
function binary_image(jd, which, λ, nx; pixsize = PIXSIZE)
    a, m1, b, m2 = models_at(jd)[which]
    # No occultation for SPHERES: `occultation_weights` reads `center_offsets`, which
    # `create_star` leaves at ZERO, so the pair comes out concentric and a disc is punched
    # out of the middle of the image.
    ow1, ow2 = which === :spheres ? (nothing, nothing) :
               occultation_weights(a, b; method = :exact)[1:2]
    i1 = a.index_quads_visible; i2 = b.index_quads_visible
    w1 = intensity(m1, :planck, λ)[i1] .* a.vis_weights[i1] .* a.ldmap[i1]
    w2 = intensity(m2, :planck, λ)[i2] .* b.vis_weights[i2] .* b.ldmap[i2]
    ow1 === nothing || (w1 = w1 .* ow1[i1])
    ow2 === nothing || (w2 = w2 .* ow2[i2])
    ox, oy = orbit_to_rotir_offset(bp, jd)
    c0 = nx ÷ 2 + 1
    return (rasterize_polygon_image(a.proj_west[i1, :], a.proj_north[i1, :], w1, pixsize, nx;
                                    cx = c0, cy = c0) .+
            rasterize_polygon_image(b.proj_west[i2, :], b.proj_north[i2, :], w2, pixsize, nx;
                                    cx = c0 + ox / pixsize, cy = c0 + oy / pixsize)) |>
           img -> img ./ sum(img)
end

# Spica (alpha Vir) — tessellated binary model and image reconstruction at the joint fit
#
#   julia -t 8 --project=demos demos/spica_full_imaging.jl
#
# Roche surfaces with von Zeipel gravity darkening and mutual irradiation, compared epoch by
# epoch with 2007_2012_2015.Spica.oifits; then a two-map image reconstruction at each epoch
# and per-epoch temperature maps. The orbit, radii and mass ratio come from
# demos/spica_refine_orbit_fit.jl (results/spica_orbitfit_joint_rv_best_levenberg.txt).
# Figures go to demos/results/spica_tess/.

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
# The orbit, radii, mass ratio and limb darkening of the tessellated joint OI + RV fit
# (demos/spica_refine_orbit_fit.jl, which describes the model, the fixed quantities and
# the result), read from demos/results/spica_orbitfit_joint_rv_best_levenberg.txt.
fit, p1, p2, bp = spica_tessellated_system()
# The meshes of demos/spica_refine_orbit_fit.jl, so the parametric chi2 below is the fit's.
t1, t2 = tessellation_healpix(4), tessellation_healpix(3)
OUT = joinpath(@__DIR__, "results", "spica_tess")
mkpath(OUT)

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
tm2 = [temperature_map_vonZeipel_roche_single(p2, s2[i], tepochs[i]; secondary=true) for i in 1:nepochs]
# mutual irradiation
kern = [reflection_kernels(s1[i], s2[i]) for i in 1:nepochs]

function heat(i, A; t1=tm1, t2=tm2, k=kern)
    A == 0 && return t1[i], t2[i]
    return handle_reflection(s1[i], t1[i], s2[i], t2[i]; albedo1=A, albedo2=A, kernels=k[i])
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
    # The model the reconstruction starts from: this epoch's Roche geometries from build()
    # and the irradiated von Zeipel maps.
    tepoch = tepochs[i]
    star1, star2 = s1[i], s2[i]
    m1, m2 = heat(i, ALBEDO_DEFAULT)

    offset_x, offset_y = orbit_to_rotir_offset(bp, tepoch)
    phase = phases[i]                          # per-observation times, as above

    v2_model, t3amp_model, t3phi_model = binary_observables(m1, star1, m2, star2, data[i], phase; intensity_model=:planck, band=bands[i], occultation=occ[i])
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
# One reconstruction per epoch: the Roche shapes and the irradiation change along the
# eccentric orbit, so each epoch gets its own pair of maps, fitted to that epoch's data only
# with that epoch's geometries (s1[i], s2[i]) and phase ramp (phases[i]).
# binary_reconstruct_oi takes vectors of epochs, so a single epoch is passed as 1-element
# vectors: data[i:i], [s1[i]], [s2[i]], [phases[i]].
#
# The reconstructed map is a linear surface brightness (binary_chi2_fg has no Planck model),
# started from the irradiated von Zeipel temperatures; in the H band I is close to
# proportional to T at 20-30 kK, so the values read approximately as temperatures.
#
# Gaussian prior ("prior") toward the irradiated von Zeipel map h1, h2 (the known
# temperatures), sigma = SIGMA_T per tessel, weight 1. The prior is in chi2 units, so the same
# weight applies at every epoch and epochs with more data move further from it. In a scan
# with HEALPix 3 for both stars (sigma = 1.5, 2.5 and 4 kK), adding "sobel2" at weight 10
# changed neither chi2r nor the maps at any epoch.
ndata_epoch = [d.nv2 + d.nt3amp + d.nt3phi for d in data]
SIGMA_T = 500.0
recon1 = Vector{Vector{Float64}}(undef, nepochs)
recon2 = Vector{Vector{Float64}}(undef, nepochs)
chi2_recon = zeros(nepochs)
for i in 1:nepochs
    setup_oi!(data[i:i], [s1[i]])
    setup_oi!(data[i:i], [s2[i]])
    h1, h2 = heat(i, ALBEDO_DEFAULT)
    x = vcat(h1, h2)
    reg1 = [["prior", 1.0, (x0 = h1, σ = SIGMA_T), 1:length(h1)]]
    reg2 = [["prior", 1.0, (x0 = h2, σ = SIGMA_T), 1:length(h2)]]
    xr = binary_reconstruct_oi(x, data[i:i], [s1[i]], [s2[i]], [phases[i]]; maxiter = 500,
                               verbose = false, regularizers1 = reg1, regularizers2 = reg2)
    r1, r2 = copy.(split_binary_map(xr, s1[i]))
    recon1[i], recon2[i] = r1, r2
    chi2_recon[i] = binary_chi2_f(r1, s1[i], r2, s2[i], data[i], phases[i])
    println("Epoch $i: reconstruction chi2r = $(round(chi2_recon[i] / ndata_epoch[i], digits = 2))")
end
println("\nReconstruction: total chi2 = $(round(sum(chi2_recon), digits = 1)), " *
        "reduced = $(round(sum(chi2_recon) / sum(ndata_epoch), digits = 2))")

# Plot every epoch on the same colour scale (min/max of the visible tessels over all epochs
# and both stars) and the same axis extent (the widest separation, plus a 0.3 mas margin).
visible_values = reduce(vcat, vcat(recon1[i][s1[i].index_quads_visible],
                                   recon2[i][s2[i].index_quads_visible]) for i in 1:nepochs)
vmin_all, vmax_all = extrema(visible_values)
star_extent(s) = maximum(hypot.(s.vertices_xyz[:, :, 1], s.vertices_xyz[:, :, 2]))
axis_max_all = maximum(1:nepochs) do i
    ow, on = orbit_to_rotir_offset(bp, tepochs[i])
    max(star_extent(s1[i]), star_extent(s2[i]) + max(abs(ow), abs(on)))
end + 0.3
for i in 1:nepochs
    plot2d_binary(recon1[i], recon2[i], s1[i], s2[i], bp, tepochs[i], colormap = "RdBu",
                  rotation_axis = true, rotation_arrow = true,
                  vmin = vmin_all, vmax = vmax_all, axis_max = axis_max_all,
                  figtitle = "Reconstruction, MJD $(round(epoch_mean_mjd[i], digits = 3)) (Epoch $i)")
    savefig(joinpath(OUT, "spica_reconstruction$i.png"))
end



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
    savefig(joinpath(OUT, "spica_prox_epoch$i.png"))
end

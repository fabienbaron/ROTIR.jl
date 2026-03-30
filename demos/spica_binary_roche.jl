using ROTIR
using DelimitedFiles, PyPlot, PyCall
import Statistics: mean

# =============================================================================
# 1. LOAD SPICA OIFITS DATA
# =============================================================================
# The merged OIFITS file contains epochs from 2007, 2012, and 2015 campaigns
oifitsfile = "./data/2007_2012_2015.Spica.oifits"
data_all = readoifits(oifitsfile)[1,1]

# Automatically identify epochs from MJD gaps in the data
all_v2_mjds = sort(data_all.v2_mjd)
gap_threshold = 0.5  # days — observations separated by > half a day are different epochs
jumps = findall(diff(all_v2_mjds) .> gap_threshold)
epoch_starts = all_v2_mjds[[1; jumps .+ 1]]
epoch_ends   = all_v2_mjds[[jumps; length(all_v2_mjds)]]
nepochs = length(epoch_starts)

# Split into per-epoch data
data = Vector{typeof(data_all)}(undef, nepochs)
epoch_mean_mjd = zeros(nepochs)
for i in 1:nepochs
    idx = set_data_filter(data_all; mjd_range=[epoch_starts[i] - 0.01, epoch_ends[i] + 0.01])
    data[i] = filter_data(data_all, idx)
    epoch_mean_mjd[i] = mean(data[i].v2_mjd)
    println("Epoch $i: MJD=$(round(epoch_mean_mjd[i], digits=4)), nV2=$(data[i].nv2), nT3=$(data[i].nt3phi)")
end

# =============================================================================
# 2. SPICA ROCHE PARAMETERS
# =============================================================================
# Spica is a close B-type binary (P=4.01d, e=0.12, q=0.62)
# Both stars are tidally locked and fill a significant fraction of their Roche lobes.
# Aufdenberg+2015: angular diameters ~0.93 mas (primary), ~0.57 mas (secondary)
# Tkachenko+2016: Tpole ~25300 K (primary), ~20585 K (secondary)
# beta = 0.25 for radiative envelopes (standard von Zeipel for B stars)

# Shared orbital parameters
P_orb    = 4.0145     # days
a_orb    = 1.54       # mas (semi-major axis of relative orbit)
e_orb    = 0.123
T0_orb   = 2454189.40 # JD (time of periastron)
q_binary = 0.6188     # M2/M1
i_orb    = 116.0      # degrees (>90 = retrograde)
Omega    = 309.938     # degrees
omega    = 255.0       # degrees (argument of periapsis, relative orbit)

# Primary Roche parameters
# q = M2/M1 for the primary potential
roche_params_1 = (
    surface_type             = 3,           # Roche lobe
    rpole                    = 0.93/2,      # mas (polar radius)
    tpole                    = 25300.0,     # K
    ldtype                   = 3,           # Hestroffer power law
    ld1                      = 0.15,
    ld2                      = 0.0,
    inclination              = 180.0 - i_orb,  # equivalent prograde viewing angle
    position_angle           = Omega - 180.0,  # rotation axis PA on sky
    rotation_period          = P_orb,       # tidally locked
    beta                     = 0.25,        # von Zeipel (radiative envelope)
    d                        = 77.0,        # pc
    q                        = q_binary,    # M2/M1
    fillout_factor_primary   = -1,          # use rpole to define potential
    fillout_factor_secondary = -1,
    i = i_orb, Ω = Omega, ω = omega,
    P = P_orb, a = a_orb, e = e_orb, T0 = T0_orb,
    dP = 0.0, dω = 0.0
)

# Secondary Roche parameters
# For compute_potential_secondary: q = M_companion/M_self = M1/M2 = 1/q_binary
roche_params_2 = (
    surface_type             = 3,
    rpole                    = 0.57/2,      # mas
    tpole                    = 20585.0,     # K
    ldtype                   = 3,
    ld1                      = 0.15,
    ld2                      = 0.0,
    inclination              = 180.0 - i_orb,
    position_angle           = Omega - 180.0,
    rotation_period          = P_orb,       # tidally locked
    beta                     = 0.25,
    d                        = 77.0,
    q                        = 1.0/q_binary, # M1/M2 (inverted for secondary potential)
    fillout_factor_primary   = -1,
    fillout_factor_secondary = -1,
    i = i_orb, Ω = Omega, ω = omega,
    P = P_orb, a = a_orb, e = e_orb, T0 = T0_orb,
    dP = 0.0, dω = 0.0
)

# Binary orbital parameters (for orbit computation and plotting)
star1p = starparameters(0.93/2, 25300.0, 0.0, 3, 0.15, 0.0, 0.25, 0.0,
                        180.0-i_orb, Omega-180.0, 0.0, P_orb)
star2p = starparameters(0.57/2, 20585.0, 0.0, 3, 0.15, 0.0, 0.25, 0.0,
                        180.0-i_orb, Omega-180.0, 0.0, P_orb)
bparams = binaryparameters(star1p, star2p,
    77.0, i_orb, Omega, omega, P_orb, a_orb, e_orb, T0_orb, q_binary,
    [1.0, 1.0], 0.0, 0.0)

# =============================================================================
# 3. TESSELLATION & GEOMETRY
# =============================================================================
n1 = 4  # HEALPix level for primary (3072 pixels)
n2 = 3  # Lower resolution for smaller secondary (768 pixels)
tessels1 = tessellation_healpix(n1)
tessels2 = tessellation_healpix(n2)

# Epoch times: use dummy t=0 for all epochs since the Roche shape depends on
# the instantaneous separation D(t), which is computed internally from the orbital elements.
# For spheres D doesn't matter, but for Roche it does (eccentric orbit → varying D).
# However, create_star uses `t` for both Roche radii AND rotation angle.
# Since we compute the binary offset separately, set t=0 for geometry.
tepochs_dummy = zeros(nepochs)

println("\nCreating primary Roche geometry...")
stars1 = create_star_multiepochs(tessels1, roche_params_1, tepochs_dummy)
println("Creating secondary Roche geometry...")
stars2 = create_star_multiepochs(tessels2, roche_params_2, tepochs_dummy; secondary=true)

# Set up polygon FTs for each star against the data
setup_oi!(data, stars1)
setup_oi!(data, stars2)

# Temperature maps (von Zeipel gravity darkening)
tmap1 = parametric_temperature_map(roche_params_1, stars1[1])
tmap2 = parametric_temperature_map(roche_params_2, stars2[1]; secondary=true)

println("Primary Tmap range: $(minimum(tmap1)) - $(maximum(tmap1)) K")
println("Secondary Tmap range: $(minimum(tmap2)) - $(maximum(tmap2)) K")

# =============================================================================
# 4. FORWARD MODEL: COMPUTE CHI2 PER EPOCH
# =============================================================================
println("\n=== Binary Roche Forward Model ===")
total_chi2 = 0.0
model_obs = Vector{NamedTuple}(undef, nepochs)

for i in 1:nepochs
    # Use precise mean MJD from the data (not the rounded filter value)
    tepoch_jd = epoch_mean_mjd[i] + 2400000.5
    offset_x, offset_y = orbit_to_rotir_offset(bparams, tepoch_jd)
    phase = binary_phase_shift(data[i].uv, offset_x, offset_y)

    v2_model, t3amp_model, t3phi_model = binary_observables(tmap1, stars1[i], tmap2, stars2[i], data[i], phase)
    model_obs[i] = (v2=v2_model, t3amp=t3amp_model, t3phi=t3phi_model, visamp=Float64[], visphi=Float64[])

    chi2_epoch = binary_chi2_f(tmap1, stars1[i], tmap2, stars2[i], data[i], phase)
    global total_chi2 += chi2_epoch

    sep = sqrt(offset_x^2 + offset_y^2)
    ndata_epoch = data[i].nv2 + data[i].nt3amp + data[i].nt3phi
    println("Epoch $i (MJD $(round(epoch_mean_mjd[i], digits=4))): sep=$(round(sep,digits=3)) mas, chi2r=$(round(chi2_epoch/ndata_epoch, digits=2))")
end
ndata_total = sum(d.nv2 + d.nt3amp + d.nt3phi for d in data)
println("\nTotal chi2 = $(round(total_chi2, digits=1)), reduced = $(round(total_chi2/ndata_total, digits=2))")

# =============================================================================
# 5. PLOT: BINARY IMAGE (epoch 6 = widest separation, 1.56 mas)
# =============================================================================
i_epoch = 6
tepoch_jd = epoch_mean_mjd[i_epoch] + 2400000.5
inc_star = 180.0 - i_orb  # stellar inclination = 64°
pa_star  = Omega - 180.0   # position angle of spin axis = 129.938°
plot2d_binary(tmap1, tmap2, stars1[i_epoch], stars2[i_epoch], bparams, tepoch_jd,
    rotation_axis=true, graticules=true,
    inclination1=inc_star, position_angle1=pa_star,
    inclination2=inc_star, position_angle2=pa_star,
    star_params1=roche_params_1, star_params2=roche_params_2,
    figtitle="Spica Binary (Roche) — Epoch $i_epoch")

# Debug: plot each component separately
plot2d(tmap1, stars1[i_epoch], intensity=true, graticules=true, rotation_axis=true, rotation_arrow=true,
    inclination=inc_star, position_angle=pa_star, star_params=roche_params_1, figtitle="Primary (Roche)")
plot2d(tmap2, stars2[i_epoch], intensity=true, graticules=true, rotation_axis=true, rotation_arrow=true,
    inclination=inc_star, position_angle=pa_star, star_params=roche_params_2, figtitle="Secondary (Roche)")

# =============================================================================
# 6. PLOT: ORBITAL DIAGRAM
# =============================================================================
fig1, ax1 = subplots(1, 1, figsize=(8, 8))
ax1.set_aspect("equal", adjustable="box")
patches_mpl = pyimport("matplotlib.patches")

tepochs_orbit = bparams.T0 .+ collect(range(0.0, stop=bparams.P, length=500))
orbit_x = zeros(500)
orbit_y = zeros(500)
for (j, t) in enumerate(tepochs_orbit)
    x1, y1, z1, x2, y2, z2 = binary_orbit_abs(bparams, t)
    orbit_x[j] = y2 - y1  # East offset
    orbit_y[j] = x2 - x1  # North offset
end
ax1.plot(orbit_x, orbit_y, "b-", linewidth=1, alpha=0.5)

# Draw primary at origin (filled circle with correct angular radius)
rpri = bparams.star1.rpole  # angular radius in mas
rsec = bparams.star2.rpole
c_pri = patches_mpl.Circle((0, 0), rpri, facecolor="gold", edgecolor="black", linewidth=0.5, zorder=5)
ax1.add_patch(c_pri)

# Mark observed epochs with secondary disk at correct scale
for i in 1:nepochs
    t_jd = epoch_mean_mjd[i] + 2400000.5
    x1, y1, z1, x2, y2, z2 = binary_orbit_abs(bparams, t_jd)
    east = y2 - y1; north = x2 - x1
    c_sec = patches_mpl.Circle((east, north), rsec, facecolor="lightskyblue",
        edgecolor="black", linewidth=0.5, zorder=4, alpha=0.8)
    ax1.add_patch(c_sec)
    ax1.annotate("$i", (east, north), textcoords="offset points", xytext=(5,5), fontsize=8)
end

ax1.invert_xaxis()
ax1.set_xlabel("East offset (mas)")
ax1.set_ylabel("North offset (mas)")
ax1.set_title("Spica Binary Orbit (secondary relative to primary)")
ax1.grid(true, alpha=0.3)

# =============================================================================
# 7. PLOT: RADIAL VELOCITIES
# =============================================================================
data_rv1 = readdlm("./data/all_rv_1_ORIG.txt")
data_rv2 = readdlm("./data/all_rv_2_ORIG.txt")

plot_rv(bparams, K1=123.9, K2=198.8, γ=0.0,
    rv_data1=data_rv1, rv_data2=data_rv2, figtitle="Spica Radial Velocities")

# =============================================================================
# 8. PLOT: V2 AND T3PHI MODEL VS DATA (per epoch)
# =============================================================================
for i in 1:nepochs
    plot_residuals(data[i], model_obs[i], figsize=(12, 8))
    suptitle("Epoch $i — MJD $(round(epoch_mean_mjd[i], digits=4))", fontsize=14)
end

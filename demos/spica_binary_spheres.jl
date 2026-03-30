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
data = Vector{typeof(data_all)}(undef, nepochs);
epoch_mean_mjd = zeros(nepochs)
for i in 1:nepochs
    idx = set_data_filter(data_all; mjd_range=[epoch_starts[i] - 0.01, epoch_ends[i] + 0.01])
    data[i] = filter_data(data_all, idx)
    epoch_mean_mjd[i] = mean(data[i].v2_mjd)
    println("Epoch $i: MJD=$(round(epoch_mean_mjd[i], digits=4)), nV2=$(data[i].nv2), nT3=$(data[i].nt3phi)")
end

# =============================================================================
# 2. SPICA PARAMETERS
# =============================================================================
# Orbital elements from Aufdenberg+2015
P_orb    = 4.0145     # days
a_orb    = 1.54       # mas (semi-major axis of relative orbit)
e_orb    = 0.123
T0_orb   = 2454189.40 # JD (time of periastron)
q_binary = 0.6188     # M2/M1
i_orb    = 116.0      # degrees (>90 = retrograde)
Omega    = 309.938     # degrees (longitude of ascending node)
omega    = 255.0       # degrees (argument of periapsis, relative orbit)

# Stellar orientation (spin aligned with orbital angular momentum, tidally locked)
inc_star = 180.0 - i_orb   # prograde equivalent viewing angle = 64°
pa_star  = Omega - 180.0    # position angle of spin axis on sky = 129.938°

# Star parameters as NamedTuples (for create_star / parametric_temperature_map)
# Angular diameters from Aufdenberg+2015: primary ~0.93 mas, secondary ~0.57 mas
# Inclination/PA don't affect sphere shape but are used for graticules and decorations
star1_params = (
    surface_type    = 0,           # Sphere
    radius          = 0.93/2,      # mas (angular radius)
    tpole           = 25300.0,     # K (Tkachenko+2016)
    ldtype          = 3,           # Hestroffer power law
    ld1             = 0.15,
    ld2             = 0.0,
    inclination     = inc_star,
    position_angle  = pa_star,
    rotation_period = P_orb        # tidally locked
)

star2_params = (
    surface_type    = 0,
    radius          = 0.57/2,      # mas
    tpole           = 20585.0,     # K
    ldtype          = 3,
    ld1             = 0.15,
    ld2             = 0.0,
    inclination     = inc_star,
    position_angle  = pa_star,
    rotation_period = P_orb
)

# Binary orbital parameters
#                    rpole       tpole   ω    ld  ld1   ld2  β_vZ  B_rot inc       PA        rot_off  P_rot
star1p = starparameters(0.93/2, 25300.0, 0.0, 3, 0.15, 0.0, 0.205, 0.0, inc_star, pa_star, 0.0, P_orb)
star2p = starparameters(0.57/2, 20585.0, 0.0, 3, 0.15, 0.0, 0.205, 0.0, inc_star, pa_star, 0.0, P_orb)

# omega = argument of periapsis of the *relative orbit* (= secondary's, astrometric convention)
bparams = binaryparameters(star1p, star2p,
    77.0,          # d: distance (pc)
    i_orb,         # i: orbital inclination (degrees; >90 = retrograde)
    Omega,         # Omega: longitude of ascending node (degrees)
    omega,         # omega: argument of periapsis (degrees)
    P_orb,         # P: orbital period (days)
    a_orb,         # a: semi-major axis (mas)
    e_orb,         # e: eccentricity
    T0_orb,        # T0: time of periastron (JD)
    q_binary,      # q: mass ratio M2/M1
    [1.0, 1.0],    # fillout factor (unused for spheres)
    0.0,           # dP: period change
    0.0            # domega: apsidal motion
)

# =============================================================================
# 3. TESSELLATION & GEOMETRY
# =============================================================================
n1 = 3  # HEALPix level for primary (768 pixels)
n2 = 2  # Lower resolution for smaller secondary (192 pixels)
tessels1 = tessellation_healpix(n1)
tessels2 = tessellation_healpix(n2)

# For spheres with no rotation, geometry is the same at every epoch
# Use dummy tepochs (all zeros)
tepochs_dummy = zeros(nepochs)
stars1 = create_star_multiepochs(tessels1, star1_params, tepochs_dummy)
stars2 = create_star_multiepochs(tessels2, star2_params, tepochs_dummy)

# Set up polygon FTs for each star against the data
setup_oi!(data, stars1)
setup_oi!(data, stars2)

# Parametric temperature maps (uniform for spheres)
tmap1 = parametric_temperature_map(star1_params, stars1[1])
tmap2 = parametric_temperature_map(star2_params, stars2[1])

# =============================================================================
# 4. FORWARD MODEL: COMPUTE CHI2 PER EPOCH
# =============================================================================
println("\n=== Binary Forward Model ===")
total_chi2 = 0.0
offsets = zeros(nepochs, 2)
model_obs = Vector{NamedTuple}(undef, nepochs)

for i in 1:nepochs
    # Use precise mean MJD from the data (not the rounded filter value)
    tepoch_jd = epoch_mean_mjd[i] + 2400000.5

    # Get secondary offset in ROTIR coordinates (West, North) in mas
    offset_x, offset_y = orbit_to_rotir_offset(bparams, tepoch_jd)
    offsets[i, :] = [offset_x, offset_y]

    # Compute phase shift for this epoch's UV coverage
    phase = binary_phase_shift(data[i].uv, offset_x, offset_y)

    # Compute model observables
    v2_model, t3amp_model, t3phi_model = binary_observables(tmap1, stars1[i], tmap2, stars2[i], data[i], phase)
    model_obs[i] = (v2=v2_model, t3amp=t3amp_model, t3phi=t3phi_model, visamp=Float64[], visphi=Float64[])

    # Compute chi2
    chi2_epoch = binary_chi2_f(tmap1, stars1[i], tmap2, stars2[i], data[i], phase)
    global total_chi2 += chi2_epoch

    sep = sqrt(offset_x^2 + offset_y^2)
    println("Epoch $i (MJD $(round(epoch_mean_mjd[i], digits=4))): sep=$(round(sep,digits=3)) mas, chi2r=$(round(chi2_epoch/(data[i].nv2+data[i].nt3amp+data[i].nt3phi), digits=2))")
end
ndata_total = sum(d.nv2 + d.nt3amp + d.nt3phi for d in data)
println("\nTotal chi2 = $(round(total_chi2, digits=1)), reduced = $(round(total_chi2/ndata_total, digits=2))")

# =============================================================================
# 5. PLOT: BINARY IMAGE (epoch 6 = widest separation)
# =============================================================================
i_epoch = 6
tepoch_jd = epoch_mean_mjd[i_epoch] + 2400000.5
plot2d_binary(tmap1, tmap2, stars1[i_epoch], stars2[i_epoch], bparams, tepoch_jd,
    intensity=true, graticules=true, compass=true,
    inclination1=inc_star, position_angle1=pa_star,
    inclination2=inc_star, position_angle2=pa_star,
    star_params1=star1_params, star_params2=star2_params,
    figtitle="Spica Binary (spheres) — Epoch $i_epoch")

# Debug: plot each component separately
plot2d(tmap1, stars1[i_epoch], intensity=true, graticules=true, rotation_axis=true, rotation_arrow=true, compass=true,
    inclination=inc_star, position_angle=pa_star, star_params=star1_params,
    figtitle="Primary (sphere)")
plot2d(tmap2, stars2[i_epoch], intensity=true, graticules=true, rotation_axis=true, rotation_arrow=true, compass=true,
    inclination=inc_star, position_angle=pa_star, star_params=star2_params,
    figtitle="Secondary (sphere)")

# =============================================================================
# 6. PLOT: ORBITAL DIAGRAM
# =============================================================================
fig1, ax1 = subplots(1, 1, figsize=(8, 8))
ax1.set_aspect("equal", adjustable="box")
patches_mpl = pyimport("matplotlib.patches")

# Full orbit curve
tepochs_orbit = bparams.T0 .+ collect(range(0.0, stop=bparams.P, length=500))
orbit_x = zeros(500)
orbit_y = zeros(500)
for (j, t) in enumerate(tepochs_orbit)
    x1, y1, z1, x2, y2, z2 = binary_orbit_abs(bparams, t)
    orbit_x[j] = y2 - y1  # East offset (for plotting)
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

ax1.invert_xaxis()  # East to the left (astronomical convention)
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

# Single-epoch imaging of Polaris (α UMi)
#
# Polaris is a Cepheid supergiant — its angular diameter (~3.3 mas) is resolved
# by optical interferometers, but it rotates too slowly for Doppler or
# multi-epoch rotational imaging.  This script demonstrates how to reconstruct
# a surface map from a single interferometric snapshot on a spherical star.

using ROTIR

# --- 1. Load data (single OIFITS file → single epoch) -----------------------
oifitsfiles = ["./data/polaris.oifits"]
data_all = readoifits_multiepochs(oifitsfiles, warn=false, verbose=false, T=Float32)
data = data_all[1, :]          # first wavelength bin, all epochs
nepochs = length(data)         # should be 1
tepochs = [0.0f0]             # single epoch: t = 0

# --- 2. Stellar geometry (sphere) -------------------------------------------
n = 4                          # HEALPix level → 3072 pixels
tessels = tessellation_healpix(n)

star_params = (
    surface_type    = 0,       # sphere
    radius          = 1.6,     # mas (approximate angular radius)
    tpole           = 6000.0,  # K (effective temperature)
    ldtype          = 3,       # Hestroffer power-law limb darkening
    ld1             = 0.24,    # LD coefficient
    ld2             = 0.0,
    inclination     = 90.0,    # degrees — pole-on for a non-rotating model
    position_angle  = 0.0,     # degrees — arbitrary for a sphere
    rotation_period = 1.0      # days — irrelevant for single epoch
)

stars = create_star_multiepochs(tessels, star_params, tepochs)
setup_oi!(data, stars)

# --- 3. Starting map ---------------------------------------------------------
tmap_start = parametric_temperature_map(star_params, stars[1])

# --- 4. Regularization -------------------------------------------------------
# TV2 (quadratic total variation) encourages smooth maps while preserving
# large-scale structure — good for resolving spots on a slow rotator.
regularizers = [["tv2", 1e-4, tv_neighbors_healpix(n), 1:length(tmap_start)]]

# --- 5. Reconstruct ----------------------------------------------------------
tmap = image_reconstruct_oi(tmap_start, data, stars,
    maxiter=500, regularizers=regularizers, verbose=true)

# --- 6. Evaluate fit ---------------------------------------------------------
chi2 = image_reconstruct_oi_chi2(tmap, data, stars, verbose=true)

# --- 7. Plot -----------------------------------------------------------------
plot2d(tmap, stars[1], intensity=true,
    graticules=true, compass=true,
    figtitle="Polaris – single-epoch reconstruction")

plot_mollweide(tmap, stars[1])

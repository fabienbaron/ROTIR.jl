# demo_rasterize_nfft.jl
#
# Demonstrates the two new forward-modeling methods added to ROTIR:
#   1. Rasterization: exact polygon-to-pixel rendering via Sutherland-Hodgman clipping
#   2. NFFT: complex visibilities on a regular Fourier grid via Gauss-Legendre + NFFT
#
# Both are alternatives to the dense polyft matrix for producing images.
# The polyft matrix approach is O(Npix * Nk * 4) trig operations and builds a
# large dense matrix; rasterization and NFFT avoid the dense matrix entirely.

using ROTIR
using FFTW
using PyPlot

# --- Create a rapid rotator star with gravity darkening --------------------
n = 4  # HEALPix refinement level
tessels = tessellation_healpix(n)

star_params = (
    surface_type    = 2,        # Rapid Rotator
    rpole           = 0.5,      # polar radius (mas)
    tpole           = 8000.0,   # polar temperature (K)
    frac_escapevel  = 0.9,      # omega = fraction of critical rotation
    ldtype          = 3,        # power-law limb darkening (Hestroffer)
    ld1             = 0.22886,
    ld2             = 0.0,
    beta            = 0.25,     # von Zeipel gravity-darkening exponent
    inclination     = 60.0,     # degrees
    position_angle  = 30.0,     # degrees
    rotation_period = 1.0,      # days (arbitrary for single epoch)
)

tepochs = [0.0]
stars = create_star_multiepochs(tessels, star_params, tepochs)
star = stars[1]

# Temperature map with von Zeipel gravity darkening
tmap = parametric_temperature_map(star_params, star)

# --- Prepare the weighted pixel intensities --------------------------------
# Only visible pixels contribute; the intensity is temperature * visibility weight
indx = star.index_quads_visible
pw = star.proj_west[indx, :]    # (nvis, 4) in mas
pn = star.proj_north[indx, :]   # (nvis, 4) in mas
x_weighted = tmap[indx] .* star.vis_weights[indx]

# --- Image parameters ------------------------------------------------------
nx = 128          # pixels per side
pixsize = 0.015f0 # mas/pixel (chosen to comfortably fit the star)

# --- Method 1: Direct rasterization ----------------------------------------
println("Rasterizing...")
@time img_raster = rasterize_polygon_image(pw, pn, x_weighted, pixsize, nx)

# --- Method 2: NFFT --------------------------------------------------------
println("Computing NFFT...")
@time img_nfft = polyft_nfft_image(pw, pn, x_weighted, pixsize, nx; ngauss=6)

# --- Plot the results -------------------------------------------------------
half_fov = nx * pixsize / 2  # mas

fig, axes = subplots(1, 2, figsize=(12, 5))

axes[1].imshow(img_raster, origin="lower",
               extent=[-half_fov, half_fov, -half_fov, half_fov],
               cmap="hot", interpolation="nearest")
axes[1].set_title("Rasterization")
axes[1].set_xlabel("West (mas)")
axes[1].set_ylabel("North (mas)")

axes[2].imshow(img_nfft, origin="lower",
               extent=[-half_fov, half_fov, -half_fov, half_fov],
               cmap="hot", interpolation="nearest")
axes[2].set_title("NFFT (ngauss=6)")
axes[2].set_xlabel("West (mas)")
axes[2].set_ylabel("North (mas)")

suptitle("Rapid Rotator (omega=0.9, inc=60 deg) -- Rasterize vs NFFT")
tight_layout()
savefig("demo_rasterize_nfft.png", dpi=150)
println("Saved demo_rasterize_nfft.png")

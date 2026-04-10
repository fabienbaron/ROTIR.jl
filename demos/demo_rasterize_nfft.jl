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

# --- Method 3: Step-by-step Fourier pipeline via rfftfreq/fftfreq -----------
# Demonstrates the full pipeline explicitly:
#   1. Define the frequency grid using rfftfreq + fftfreq
#   2. Compute complex visibilities via NFFT on that grid
#   3. Inverse real-FFT (irfft) back to a real-space image

# Step 1: Build the frequency grids corresponding to the image grid.
# For an image of nx pixels with spacing pixsize (mas/pixel), the spatial
# sampling rate is 1/pixsize (samples/mas), giving frequencies in cycles/mas.
u_freq = rfftfreq(nx, 1 / pixsize)    # non-negative freqs, length nx÷2+1
v_freq = fftfreq(nx, 1 / pixsize)     # full freqs,         length nx

# Step 2: Compute complex visibilities via adjoint NFFT.
# polyft_nfft_forward returns F of shape (nx÷2+1, nx) in standard rfft layout:
#   dim 1 (rows) matches rfftfreq — non-negative u frequencies
#   dim 2 (cols) matches fftfreq  — standard FFT order
println("Computing visibilities on rfftfreq/fftfreq grid via NFFT...")
@time F = polyft_nfft_forward(pw, pn, x_weighted, pixsize, nx; ngauss=6)

# Step 3: Inverse real-FFT to recover the image from the visibilities.
img_irfft = fftshift(irfft(F, nx))

# For display, fftshift dim 2 so zero-frequency is centered
F_display = fftshift(F, 2)
v_freq_shifted = fftshift(v_freq)

# --- Plot the results -------------------------------------------------------
half_fov = nx * pixsize / 2  # mas

# Frequency axis extents for the visibility plot (fftshifted v on x-axis)
u_extent = [Float64(v_freq_shifted[1]), Float64(v_freq_shifted[end]),
            Float64(u_freq[1]),         Float64(u_freq[end])]

fig, axes = subplots(2, 2, figsize=(12, 10))

# Top-left: rasterization
axes[1,1].imshow(img_raster, origin="lower",
                 extent=[-half_fov, half_fov, -half_fov, half_fov],
                 cmap="hot", interpolation="nearest")
axes[1,1].set_title("Rasterization")
axes[1,1].set_xlabel("West (mas)")
axes[1,1].set_ylabel("North (mas)")

# Top-right: NFFT convenience image (polyft_nfft_image)
axes[1,2].imshow(img_nfft, origin="lower",
                 extent=[-half_fov, half_fov, -half_fov, half_fov],
                 cmap="hot", interpolation="nearest")
axes[1,2].set_title("NFFT image (ngauss=6)")
axes[1,2].set_xlabel("West (mas)")
axes[1,2].set_ylabel("North (mas)")

# Bottom-left: visibility amplitudes |F| on the rfftfreq x fftfreq grid
axes[2,1].imshow(log10.(abs.(F_display) .+ 1f-12), origin="lower",
                 extent=u_extent, aspect="auto",
                 cmap="inferno", interpolation="nearest")
axes[2,1].set_title("log₁₀|V| on rfft grid")
axes[2,1].set_xlabel("v  (cycles/mas)")
axes[2,1].set_ylabel("u  (cycles/mas)")

# Bottom-right: image recovered via explicit irfft
axes[2,2].imshow(img_irfft, origin="lower",
                 extent=[-half_fov, half_fov, -half_fov, half_fov],
                 cmap="hot", interpolation="nearest")
axes[2,2].set_title("irfft(F) image")
axes[2,2].set_xlabel("West (mas)")
axes[2,2].set_ylabel("North (mas)")

suptitle("Rapid Rotator (ω=0.9, inc=60°) — Rasterize / NFFT / Fourier pipeline")
tight_layout()
savefig("demo_rasterize_nfft.png", dpi=150)
println("Saved demo_rasterize_nfft.png")

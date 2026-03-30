using ROTIR
using PyPlot

# Simple limb-darkened sphere: inclination 35°, PA 20°
star_params = (
    surface_type    = 0,       # Sphere
    radius          = 1.0,     # mas
    tpole           = 10000.0, # K
    ldtype          = 3,       # Hestroffer power law
    ld1             = 0.3,
    ld2             = 0.0,
    inclination     = 35.0,    # degrees from line of sight
    position_angle  = 20.0,    # degrees, N through E
    rotation_period = 1.0      # day (arbitrary)
)

n = 4  # HEALPix level (3072 pixels)
tessels = tessellation_healpix(n)
star = create_star(tessels, star_params, 0.0)
tmap = parametric_temperature_map(star_params, star)

plot2d(tmap, star, intensity=true,
    graticules=true, rotation_axis=true, rotation_arrow=true, compass=true,
    inclination=35.0, position_angle=20.0,
    star_params=star_params,
    figtitle="Sanity check: inc=35° PA=20°")

# Interactive 3-D views of ROTIR surfaces.
#
#     julia --project=. demos/view3d.jl
#
# GLMakie, not CairoMakie: this demo is about the mouse. Left-drag orbits the camera about the
# star, scroll zooms, right-drag pans. Every figure opens looking down the line of sight —
# East left, North up, exactly the view `plot2d_makie` draws — so the first thing you see is
# the sky projection, and turning it away from that is what shows you the far hemisphere.
using ROTIR, GLMakie

# ── 1. a rapid rotator ───────────────────────────────────────────────────────────────────
rr = (surface_type = 2, rpole = 1.5, tpole = 8000.0, ldtype = 3, ld1 = 0.2, ld2 = 0.0,
      inclination = 55.0, position_angle = 25.0, rotation_period = 1.0, beta = 0.25,
      frac_escapevel = 0.85, B_rot = 0.0)
star = create_star(tessellation_healpix(4), rr, 0.0)
tmap = temperature_map_vonZeipel_rapid_rotator(rr, star)

fig1, _ = plot3d_makie(tmap, star; figtitle = "rapid rotator — drag to rotate")
display(GLMakie.Screen(), fig1)

# ── 2. a Roche binary, with a phase slider ───────────────────────────────────────────────
# Spica: P = 4.0145 d, e = 0.123, q = 0.62, both components tidally locked.
P, a, e, T0, q, inc, Ω, ω = 4.0145, 1.54, 0.123, 2454189.40, 0.6188, 116.0, 309.938, 255.0
base = (surface_type = 3, ldtype = 3, ld1 = 0.15, ld2 = 0.0, inclination = 180.0 - inc,
        position_angle = Ω - 180.0, rotation_period = P, beta = 0.25, d = 77.0,
        fillout_factor_primary = -1, fillout_factor_secondary = -1,
        i = inc, Ω = Ω, ω = ω, P = P, a = a, e = e, T0 = T0, dP = 0.0, dω = 0.0)
p1 = merge(base, (rpole = 0.93/2, tpole = 25300.0, q = q))
p2 = merge(base, (rpole = 0.57/2, tpole = 20585.0, q = 1/q))
sp1 = starparameters(0.93/2, 25300.0, 0.0, 3, 0.15, 0.0, 0.25, 0.0, 180.0-inc, Ω-180.0, 0.0, P)
sp2 = starparameters(0.57/2, 20585.0, 0.0, 3, 0.15, 0.0, 0.25, 0.0, 180.0-inc, Ω-180.0, 0.0, P)
bparams = binaryparameters(sp1, sp2, 77.0, inc, Ω, ω, P, a, e, T0, q, [1.0, 1.0], 0.0, 0.0)

tess1, tess2 = tessellation_healpix(4), tessellation_healpix(3)
fig2, _ = plot3d_binary_makie(
    parametric_temperature_map(p1, create_star(tess1, p1, T0)),
    parametric_temperature_map(p2, create_star(tess2, p2, T0; secondary = true);
                               secondary = true),
    create_star(tess1, p1, T0), create_star(tess2, p2, T0; secondary = true),
    bparams, T0; figtitle = "Spica — periastron")
display(GLMakie.Screen(), fig2)

# ── 3. the same binary driven by a slider ────────────────────────────────────────────────
#
# The meshes are created ONCE and only their vertex Observables are reassigned as the phase
# moves. That is not merely an optimisation: inserting a new plot into a scene that is already
# on screen reallocates GPU buffers, which is the failure mode the GUI has to avoid entirely.
# Doing it the right way here keeps this demo honest about the pattern the GUI will use.
fig3 = Figure(size = (900, 760))
sl   = Makie.Slider(fig3[2, 1], range = range(0, 1, length = 201), startvalue = 0.0)
phase = sl.value

meshes = map(phase) do φ
    t  = T0 + φ * P
    s1 = create_star(tess1, p1, t)
    s2 = create_star(tess2, p2, t; secondary = true)
    x1, y1, z1, x2, y2, z2 = binary_orbit_abs(bparams, t)
    m1, c1 = star_mesh(s1; values = parametric_temperature_map(p1, s1))
    m2, c2 = star_mesh(s2; offset = (-(y2 - y1), x2 - x1, -(z2 - z1)),
                       values = parametric_temperature_map(p2, s2; secondary = true))
    (m1, c1, m2, c2)
end

_, ls = ROTIR.scene3d(fig3, fig3[1, 1]; title = "Spica — drag the slider through one period",
                      radius = a)
crange = (20000.0, 25600.0)
cmap   = Makie.cgrad(Makie.cgrad(:gist_heat)[range(0.15, 0.95, length = 256)])
mesh!(ls, map(m -> m[1], meshes); color = map(m -> m[2], meshes),
      colormap = cmap, colorrange = crange, shading = NoShading)
mesh!(ls, map(m -> m[3], meshes); color = map(m -> m[4], meshes),
      colormap = cmap, colorrange = crange, shading = NoShading)
lines!(ls, ROTIR.relative_orbit_track(bparams); color = (:grey35, 0.8), linewidth = 1.2)
Colorbar(fig3[1, 2]; colormap = cmap, colorrange = crange, label = "T (K)")
display(GLMakie.Screen(), fig3)

println("Three windows open. Left-drag rotates, scroll zooms, right-drag pans.")
println("Press Enter to close.")
readline()

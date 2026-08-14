# =====================================================================
# Plotting regression tests
# =====================================================================
# Every plotting entry point ROTIR exports, called with every decoration and option,
# asserting on the STRUCTURE of what lands on the axes.
#
# Why structural assertions and not image comparison: reference PNGs are hostage to the
# matplotlib version, the font stack and the freetype build, so they produce false
# failures on upgrade and get disabled within a release or two. What actually breaks in
# practice is the Julia/Python boundary: a colour array arriving as `Vector{Py}` instead
# of an (n,4) float array, a
# `Matrix` silently read column-major, a `zip` becoming a numpy *structured* array, a
# Julia `Symbol` reaching matplotlib as `:black` instead of `"black"`, a removed API
# (`matplotlib.cm.get_cmap`). Every one of those either raises or produces a collection
# with the wrong shape. Both are caught here; neither needs a reference image.
#
# So the assertions are of two kinds:
#   1. it runs at all                    -> catches removed APIs and type-conversion errors
#   2. the artists have the right shape  -> catches silent conversion corruption
#      (one PolyCollection per star, one path per visible tessel, facecolors (n,4))
#
# Runs by default: Agg needs no display and the whole file is a few seconds.
# Set ROTIR_TEST_FIGURES=1 to also write the PNGs to test/figures/plotting/ for eyeballing.
#
# NOT covered, deliberately:
#   * src/lciplot.jl and src/specutils.jl — neither is `include`d by src/ROTIR.jl, so
#     their functions do not exist at runtime. Wire them in and they belong here too.
#   * plot_v2_residuals / plot_t3*_residuals — re-exported from OITOOLS, upstream's to test.
#   * draw_graticules on surface_type 1 (ellipsoid) and 2 (rapid rotator): only the
#     surface_type 0 exact branch, surface_type 3 (Roche) and the mesh-derived path for a
#     star with no star_params are exercised below.

using Test
using ROTIR
using LinearAlgebra
using PythonCall
using PythonPlot

const mpl = pyimport("matplotlib")
mpl.use("Agg")                     # force, whatever MPLBACKEND resolved to
const plt = PythonPlot.pyplot

const SAVE_FIGS = get(ENV, "ROTIR_TEST_FIGURES", "0") == "1"
const FIGDIR    = joinpath(@__DIR__, "figures", "plotting")
if SAVE_FIGS
    # Wipe the directory first. A renamed or deleted test case would otherwise leave its
    # PNG behind, indistinguishable from a current one — so the gallery would show a
    # figure produced by code that no longer exists, and comparing it against a live
    # figure invents a disagreement that is not in the code.
    isdir(FIGDIR) && for f in readdir(FIGDIR)
        endswith(f, ".png") && rm(joinpath(FIGDIR, f))
    end
    mkpath(FIGDIR)
end

# ---------------------------------------------------------------------
# Introspection helpers
# ---------------------------------------------------------------------
ncoll(ax)  = pylen(ax.collections)
nlines(ax) = pylen(ax.lines)
ntexts(ax) = pylen(ax.texts)
nfigs()    = pylen(plt.get_fignums())
npaths(c)  = pylen(c.get_paths())
fcshape(c) = pyconvert(Vector{Int}, c.get_facecolor().shape)
pytype(o)  = pyconvert(String, pybuiltins.type(o).__name__)

"""
First axes of the current figure — the star. ROTIR builds its colorbars with
`make_axes_locatable(...).append_axes(...)`, so `plt.gca()` afterwards is the COLORBAR
axes, whose QuadMesh/LineCollection would be mistaken for the surface.
"""
main_ax() = plt.gcf().get_axes()[0]

"""
The first *filled* PolyCollection on `ax` — the tessellated surface. Skips colorbar
artists (QuadMesh/LineCollection) and stroked decorations such as graticules and rotation
arrows, which are also PolyCollections but carry an empty facecolor array.
"""
function poly_collection(ax)
    for k in 0:(ncoll(ax) - 1)
        c = ax.collections[k]
        pytype(c) == "PolyCollection" && fcshape(c)[1] > 0 && return c
    end
    return nothing
end

"""
Assert `ax` carries a tessellated surface of `n` polygons whose facecolors are a proper
RGBA array — either one row per tessel, or a single row when the whole surface is drawn
in one uniform colour (the wireframe case).

This is the most valuable check in the file: it separates "the call returned without
raising" from "the colours actually arrived across the language boundary".
"""
function check_surface(ax, n; per_tessel = true)
    c = poly_collection(ax)
    @test c !== nothing
    c === nothing && return
    @test npaths(c) == n
    fs = fcshape(c)
    @test length(fs) == 2 && fs[2] == 4          # RGBA, not a flattened/structured array
    @test per_tessel ? fs[1] == n : fs[1] in (1, n)
    check_visible(ax, c)
end

"""
The surface must lie INSIDE the viewport. Structure alone does not prove visibility: a
collection with a perfect (n,4) colour array drawn outside the axis limits renders a
blank figure and every structural assertion still passes. Overlap, not containment —
`plot2d(...; xlim=)` deliberately crops.
"""
function check_visible(ax, c)
    dl = c.get_datalim(ax.transData)
    x0, x1 = pyconvert(Float64, dl.x0), pyconvert(Float64, dl.x1)
    y0, y1 = pyconvert(Float64, dl.y0), pyconvert(Float64, dl.y1)
    xl = extrema(pyconvert(Vector{Float64}, ax.get_xlim()))
    yl = extrema(pyconvert(Vector{Float64}, ax.get_ylim()))
    @test x1 >= xl[1] && x0 <= xl[2]
    @test y1 >= yl[1] && y0 <= yl[2]
end


"""
A bare axes set up the way ROTIR's own plots are: equal aspect, symmetric limits, and the
SAME tick spacing on both axes via `set_tick_spacing`. Without this matplotlib picks ticks
per-axis (0.5 in x, 0.25 in y), which makes a circular graticule look non-isotropic even
though the aspect is equal.
"""
function deco_axes(; amax = 1.5, star = nothing)
    fig, ax = plt.subplots()
    ax.set_aspect("equal")
    ax.set_xlim(amax, -amax)      # E to the left, as in plot2d
    ax.set_ylim(-amax, amax)
    ROTIR.set_tick_spacing(ax, amax)
    # Pass `star` to outline the disk. A decoration drawn on a blank axes cannot be
    # judged: whether the spin arrow sits at the right pole, or the axis pierces the
    # limb where it should, is only answerable against the star's silhouette.
    star === nothing || draw_limb!(ax, star; color = "lightgrey")
    return fig, ax
end

"""
Output resolution for the saved figures, in dots per inch.

`bbox_inches="tight"` crops to the axes, so the equal-aspect decoration panels come out
small: at 90 dpi `deco_arrow_retrograde.png` was only 385x377, which is why zooming into
it to check an arrowhead showed nothing but interpolation. Override with
`ROTIR_TEST_FIGURE_DPI` — note the pixel count, and roughly the file size, grows with the
SQUARE of this: 300 dpi is ~11x the data of 90 (about 30 MB for the set), 600 dpi ~44x
(about 63 MB).
"""
const FIG_DPI = parse(Float64, get(ENV, "ROTIR_TEST_FIGURE_DPI", "300"))

"Close every figure, optionally saving the current one first."
function finish(name)
    if SAVE_FIGS && nfigs() > 0
        plt.savefig(joinpath(FIGDIR, "$(name).png"), dpi = FIG_DPI, bbox_inches = "tight")
    end
    plt.close("all")
end

# ---------------------------------------------------------------------
# Fixtures — small meshes, at the low end of the useful nside range
# ---------------------------------------------------------------------
sph(; r = 1.0, pa = 30.0, inc = 60.0, ldtype = 1) =
    (surface_type = 0, radius = r, tpole = 4000.0, ldtype = ldtype,
     ld1 = 0.4, ld2 = 0.1, inclination = inc, position_angle = pa,
     rotation_period = 10.0)

"""
A map with no symmetry: an m=3 sectoral pattern plus an off-centre cool spot. A flipped
axis, transposed colour array or scrambled face order shows up by eye in the saved PNGs,
and the asymmetry guarantees adjacent tessels differ so the 'colours vary' check is sharp.
"""
function testmap(star)
    θ = star.vertices_spherical[:, 5, 2]
    φ = star.vertices_spherical[:, 5, 3]
    return 4000.0 .+ 700.0 .* sin.(3φ) .* sin.(θ).^2 .-
           500.0 .* exp.(-((θ .- 0.7).^2 .+ (φ .- 1.0).^2) ./ 0.05)
end

const TES  = tessellation_healpix(3)
const PAR  = sph()
const STAR = create_star(TES, PAR, 0.0)
const TMAP = testmap(STAR)
const NVIS = STAR.nquads_visible

# lat/long mesh — the other tessellation_type branch
const NTHETA, NPHI = 20, 40
const STARLL = create_star(tessellation_latlong(NTHETA, NPHI), PAR, 0.0)
const TMAPLL = testmap(STARLL)

# Spica-like binary, reused for plot2d_binary / plot_rv / animation
const S1 = starparameters(0.45, 25300.0, 0.0, 1, 0.0, 0.0, 0.25, 0.0, 64.0, 130.0, 0.0, 4.0145)
const S2 = starparameters(0.23, 24000.0, 0.0, 1, 0.0, 0.0, 0.25, 0.0, 64.0, 130.0, 0.0, 4.0145)
const BP = binaryparameters(S1, S2, 77.0, 116.0, 309.9, 255.0, 4.0145, 1.54, 0.123,
                            2454189.4, 0.6188, [1.0, 1.0], 0.0, 0.0)
# Component temperatures: 25300 K and 24000 K.
#
# Two constraints pull against each other here, and this pair satisfies both.
#   * They must not be EQUAL. sph()'s 4000 K default gave both components the same
#     temperature, the maps spanned ~7 K end to end, the shared colour scale collapsed and
#     gist_heat rendered both stars white on a white background.
#   * They must not be FAR APART. Spica's true 25300/20585 pair spans 4715 K, which on a
#     shared scale compresses each star's own structure to 2.5% (primary) and 21.7%
#     (secondary) of the colour range — the surfaces read as flat.
# At 25300/24000 the shared range is 1509 K, and irradiation structure occupies 13.9% of
# it on the primary and 45.9% on the secondary, so both surfaces are legible while the two
# components remain clearly distinct. This is a FIGURE-LEGIBILITY choice for the plotting
# tests, not Spica's real parameters — nothing here asserts physics.
sp1() = merge(sph(r = 0.45, pa = 130.0), (tpole = 25300.0,))
sp2() = merge(sph(r = 0.23, pa = 130.0), (tpole = 24000.0,))
rp(r, tpole) = merge(sph(r = r, pa = 130.0),
              (tpole = tpole,
               i = 116.0, Omega = 309.9, omega = 255.0, P = 4.0145, a = 1.54, e = 0.123,
               T0 = 2454189.4, q = 0.6188, dP = 0.0, dOmega = 0.0))
const RP1, RP2 = rp(0.45, 25300.0), rp(0.23, 24000.0)

# A single, strongly distorted Roche star — the docs' `surface_roche.png` geometry, at the
# minimum useful nside. rpole/a = 0.355 against an Eggleton lobe radius of 0.38a for q = 1,
# so it is close to filling and the tidal teardrop is pronounced: this is the case the
# graticule fallback could not draw.
const ROCHE_PAR = (surface_type = 3, rpole = 0.355, tpole = 4800.0, ldtype = 3,
                   ld1 = 0.46, ld2 = 0.0, inclination = 60.0, position_angle = 0.0,
                   rotation_period = 5.0, beta = 0.08, d = 77.0, q = 1.0,
                   fillout_factor_primary = -1, i = 0.0, Ω = 0.0, ω = 0.0, P = 5.0,
                   a = 1.0, e = 0.0, T0 = 0.0, dP = 0.0, dω = 0.0)
const ROCHE_STAR = create_star(tessellation_healpix(3), ROCHE_PAR, 0.0)
const TEP = 2454189.4
const BS1, BS2 = create_binary_geometry(tessellation_healpix(3), RP1,
                                        tessellation_healpix(3), RP2, BP, TEP)
# Use the IRRADIATED maps. A surface_type=0 sphere with no gravity darkening is exactly
# isothermal (spread 0.0 K), so the intrinsic maps render each component as one flat
# colour and show nothing about the surface. Mutual irradiation is what the binary path
# exists for, and it puts real one-sided structure on the maps: +1051.8 K across the
# secondary, +120.8 K across the primary.
const BT1, BT2 = handle_reflection(BS1, parametric_temperature_map(sp1(), BS1),
                                   BS2, parametric_temperature_map(sp2(), BS2;
                                                                   secondary = true);
                                   albedo1 = 0.6, albedo2 = 0.6)

@testset "plotting" begin

    # -----------------------------------------------------------------
    @testset "decorations in isolation" begin
        ROTIR.set_oiplot_defaults()   # unexported; OITOOLS exports a same-named function

        # Compass: two arrows drawn as annotations, so they live in ax.texts, not ax.lines.
        fig, ax = deco_axes()          # limits come from deco_axes; amax must match below
        draw_compass(ax, 1.5)
        @test ntexts(ax) >= 2                            # "N" and "E"
        finish("deco_compass")

        # The two ways of getting the spin axis MUST agree. They did not: the mesh
        # fallback took the first/last pixels BY INDEX, which is the pole only under
        # HEALPix RING ordering. ROTIR is NESTED, where pixels 1:4 sit near the EQUATOR
        # (colatitude ~1.3-1.5 rad), so the axis came out 76.5 deg off with the north
        # arrow pointing away from the observer. Pin them together.
        let (na, sa) = ROTIR._spin_axis(STAR, PAR, 60.0, 30.0),      # analytic
            (nm, sm) = ROTIR._spin_axis(STAR, PAR, NaN, NaN)         # mesh-derived
            va = normalize(na .- sa)
            vm = normalize(nm .- sm)
            @test acosd(clamp(dot(va, vm), -1, 1)) < 1.0
            @test sign(va[3]) == sign(vm[3])          # same hemisphere, not anti-parallel
        end

        # Same story for the graticule: `star_params` carries inclination/position_angle,
        # so passing it alone must orient the grid identically to passing them explicitly.
        # It did not — star_params was read for the SHAPE but not the ORIENTATION, so the
        # graticule came out pole-up on an inclined star.
        function graticule_signature(kw)
            fig, ax = deco_axes()
            draw_graticules(ax, STAR; kw...)
            paths = ax.collections[0].get_paths()   # a PolyCollection, not a LineCollection
            n = pylen(paths)
            s = sum(sum(abs, pyconvert(Array, paths[k].vertices)) for k in 0:(n - 1))
            plt.close("all")
            return (n, s)
        end
        # All THREE routes to the orientation must agree: explicit angles, angles read from
        # star_params, and — with neither available — the rotation recovered from the mesh
        # itself. Defaulting the last one to 0/0 would draw a pole-up graticule on a star
        # inclined at 60 deg, i.e. a different star from the one plotted.
        let explicit = graticule_signature((star_params = PAR, inclination = 60.0,
                                            position_angle = 30.0)),
            from_par = graticule_signature((star_params = PAR,)),
            from_mesh = graticule_signature((;))
            # star_params carries no orientation authority — with no explicit angles both
            # of these take the Procrustes fit, so the ORIENTATION is shared
            # exactly and the curves must lie on top of each other.
            #
            # Not bitwise equal, though, and the residual is instructive. `from_par` knows
            # the surface is a sphere and clips at z > 0, which for a sphere is the exact
            # terminator. `from_mesh` has no closed form, so it interpolates the mesh — and
            # clips on the mesh's own normals, which on an nside=3 sphere are the quad
            # diagonal cross products, radial only to discretisation. The two terminators
            # differ by a fraction of a tessel. That is the intended trade: the mesh normal
            # is what `draw_limb!` thresholds on, so the graticule ends where the drawn
            # silhouette does rather than where an ideal sphere's would.
            @test isapprox(from_par[2], from_mesh[2]; rtol = 5e-3)
            @test from_par[1] == from_mesh[1]
            # ...and must reproduce the explicit angles to within the Float32 mesh fit.
            @test from_mesh[1] == explicit[1]                       # same arc count
            @test isapprox(from_mesh[2], explicit[2]; rtol = 5e-3)  # same geometry
        end

        # Rotation axis: mesh-derived branch and the analytic inclination/PA branch
        for (nm, kw) in (("mesh", (;)),
                         ("analytic", (inclination = 60.0, position_angle = 30.0,
                                       star_params = PAR)))
            fig, ax = deco_axes(star = STAR)
            n0 = nlines(ax)                              # deco_axes already drew the limb
            draw_rotation_axis(ax, STAR; kw...)
            @test nlines(ax) - n0 == 3                   # south tip / interior / north tip
            @test ntexts(ax) >= 1                        # arrowhead annotation
            finish("deco_axis_$nm")
        end

        # Rotation arrow: both poles, both branches, plus the retrograde sweep flip
        for pole in ("N", "S"), (nm, kw) in (("mesh", (;)),
                                             ("analytic", (inclination = 60.0,
                                                           position_angle = 30.0,
                                                           star_params = PAR)))
            fig, ax = deco_axes(star = STAR)
            n0 = nlines(ax)
            draw_rotation_arrow(ax, STAR; pole = pole, kw...)
            @test nlines(ax) - n0 == 2                   # front (solid) + behind (dashed)
            finish("deco_arrow_$(pole)_$nm")
        end
        # Retrograde rotation reverses the sweep direction (e1/e2 swap)
        fig, ax = deco_axes(star = STAR)
        nr0 = nlines(ax)
        retro = merge(PAR, (rotation_period = -10.0,))
        draw_rotation_arrow(ax, STAR; star_params = retro)
        @test nlines(ax) - nr0 == 2
        finish("deco_arrow_retrograde")

        # Graticules: exact (surface_type known), mesh fallback, analytic, and dense
        for (nm, kw) in (("fallback", (;)),
                         ("exact", (star_params = PAR,)),
                         ("analytic", (inclination = 60.0, position_angle = 30.0,
                                       star_params = PAR)),
                         ("dense", (nlat = 9, nlon = 16, star_params = PAR)),
                         ("styled", (color = "cyan", linewidth = 1.5, alpha = 0.9,
                                     star_params = PAR)))
            fig, ax = deco_axes()
            draw_graticules(ax, STAR; kw...)
            @test ncoll(ax) >= 1                         # PolyCollection of visible arcs
            finish("deco_grat_$nm")
        end

        # --- Roche graticules: the surface read off the mesh -------------------------
        #
        # `draw_graticules` has closed forms for surface types 0, 1 and 2. Everything else
        # interpolates the mesh's own r(θ, φ). A biaxial ellipsoid fitted to mesh-averaged
        # polar and equatorial radii will NOT do: it is axisymmetric about the spin axis and
        # therefore structurally incapable of showing the one feature a Roche figure exists
        # to show.
        let fld = ROTIR._mesh_surface_field(ROCHE_STAR)
            r_point = ROTIR._interp_surface(fld, π/2, 0.0)[1]      # toward the companion
            r_back  = ROTIR._interp_surface(fld, π/2, π)[1]        # away from it
            r_side  = ROTIR._interp_surface(fld, π/2, π/2)[1]      # along the orbit
            r_pole  = ROTIR._interp_surface(fld, 1e-6, 0.0)[1]

            # The Roche ordering. Any axisymmetric fit gives r_point == r_back == r_side,
            # so this single chain is what separates the mesh path from an axisymmetric one.
            @test r_point > r_back > r_side > r_pole
            @test r_point / r_side > 1.2               # a pronounced teardrop, not a nudge
            @test isapprox(r_pole, ROCHE_PAR.rpole; rtol = 2e-2)

            # Interpolation is a normalised weighted mean, so it can never leave the range
            # of the samples — no overshoot on the steep gradient near L1.
            allr = ROCHE_STAR.vertices_spherical[:, 5, 1]
            for (θq, φq) in ((0.3, 2.0), (1.1, 4.5), (2.4, 0.9), (π/2, 0.05))
                @test minimum(allr) <= ROTIR._interp_surface(fld, θq, φq)[1] <= maximum(allr)
            end

            # Querying exactly at a sample must return that sample, not a smoothed value.
            # (Float64 out, Float32 mesh in, so compare by value not by ===.)
            k = 7
            let (rk, nk) = ROTIR._interp_surface(fld, ROCHE_STAR.vertices_spherical[k, 5, 2],
                                                      ROCHE_STAR.vertices_spherical[k, 5, 3])
                @test rk == Float64(allr[k])
                @test nk == Float64(ROCHE_STAR.normals[k, 3])
            end

            # Normalised weights reproduce a constant exactly: a sphere stays a sphere,
            # so routing type 0 through the mesh path would not have degraded it.
            let sfld = ROTIR._mesh_surface_field(STAR)
                @test all(isapprox(ROTIR._interp_surface(sfld, θq, φq)[1], PAR.radius;
                                   rtol = 1e-10)
                          for (θq, φq) in ((0.4, 1.0), (1.9, 3.3), (2.8, 5.7)))
            end
        end

        # Visibility comes from the interpolated NORMAL, not from z > 0. Those differ on a
        # distorted body, and the normal is what `draw_limb!` clips on, so the curves stop
        # exactly where the drawn silhouette does. Passing the star's own angles must not
        # change that — it reconstructs the same rotation, which the code checks for.
        let sig(kw) = begin
                fig, ax = deco_axes()
                draw_graticules(ax, ROCHE_STAR; kw...)
                paths = ax.collections[0].get_paths()
                n = pylen(paths)
                s = sum(sum(abs, pyconvert(Array, paths[j].vertices)) for j in 0:(n - 1))
                plt.close("all")
                (n, s)
            end
            a = sig((star_params = ROCHE_PAR,))
            b = sig((star_params = ROCHE_PAR, inclination = ROCHE_PAR.inclination,
                     position_angle = ROCHE_PAR.position_angle))
            @test a[1] == b[1]                          # same arcs, same clip
            @test isapprox(a[2], b[2]; rtol = 1e-4)     # Procrustes vs rot_vertex only
        end

        fig, ax = deco_axes()
        draw_graticules(ax, ROCHE_STAR; star_params = ROCHE_PAR)
        @test ncoll(ax) >= 1
        finish("deco_grat_roche")

        # Tick spacing helper across four decades of axis size (unexported)
        for amax in (0.05, 0.5, 1.5, 15.0, 150.0)
            fig, ax = deco_axes()
            ROTIR.set_tick_spacing(ax, amax)
            @test pylen(ax.get_xticks()) >= 2
            plt.close("all")
        end

        # add_tessel_collection!: scalar colour and per-tessel colours
        fig, ax = deco_axes()
        add_tessel_collection!(ax, STAR, "white"; edgecolors = "black", linewidths = 0.5)
        check_surface(ax, NVIS; per_tessel = false)
        finish("deco_tessels_scalar")

        fig, ax = deco_axes()
        add_tessel_collection!(ax, STAR, [[0.2, 0.4, 0.6, 1.0] for _ in 1:NVIS];
                               plotmesh = true)
        check_surface(ax, NVIS)
        finish("deco_tessels_vector")

        # Single-quad helper (unexported)
        fig, ax = deco_axes()
        ROTIR.plot2Dquad(STAR, 1)
        plt.close("all")
    end

    # -----------------------------------------------------------------
    @testset "plot2d" begin
        cases = Dict(
            "default"       => (;),
            "intensity"     => (intensity = true,),
            "mesh"          => (plotmesh = true,),
            "compass_off"   => (compass = false,),
            "graticules"    => (graticules = true, star_params = PAR),
            "rot_axis"      => (rotation_axis = true, star_params = PAR),
            "rot_arrow"     => (rotation_arrow = true, star_params = PAR),
            "all_deco"      => (graticules = true, rotation_axis = true,
                                rotation_arrow = true, compass = true, star_params = PAR),
            "flipx"         => (flipx = true,),
            "limits"        => (xlim = [-2.0, 2.0], ylim = [-2.0, 2.0]),
            "pad"           => (pad = 1.5,),
            "background"    => (background = "black",),
            "colormap"      => (colormap = "viridis",),
            "contours"      => (contours = [3600.0, 4000.0, 4400.0],),
            "contour_style" => (contours = [4000.0], contour_color = "cyan",
                                contour_labels = false, contour_fontsize = 8),
            "planck"        => (intensity_model = :planck, band = 1.65e-6),
            "title"         => (figtitle = "regression",),
            "analytic_geom" => (inclination = 60.0, position_angle = 30.0,
                                graticules = true, rotation_axis = true,
                                rotation_arrow = true, star_params = PAR),
        )
        for (nm, kw) in cases
            plot2d(TMAP, STAR; kw...)
            check_surface(main_ax(), NVIS)
            finish("plot2d_$nm")
        end

        # The colours must actually VARY across tessels. A conversion bug that collapses
        # or broadcasts one colour would still give the right array shape.
        plot2d(TMAP, STAR)
        fc = pyconvert(Array, poly_collection(main_ax()).get_facecolor())
        plt.close("all")
        @test size(fc) == (NVIS, 4)
        @test length(unique(eachrow(fc))) > NVIS ÷ 10

        # The intensity model must actually bite — guards against the toggle silently
        # doing nothing. Note it only applies when `intensity=true`: with `intensity=false`
        # the surface is coloured by TEMPERATURE, which is band-independent by definition,
        # so linear and Planck must agree there and differ once brightness is plotted.
        facecolors(; kw...) = begin
            plot2d(TMAP, STAR; kw...)
            a = pyconvert(Array, poly_collection(main_ax()).get_facecolor())
            plt.close("all"); a
        end
        t_lin = facecolors(intensity = false, intensity_model = :linear)
        t_pla = facecolors(intensity = false, intensity_model = :planck, band = 1.65e-6)
        @test t_lin == t_pla                       # temperature map: band-independent

        i_lin = facecolors(intensity = true, intensity_model = :linear)
        i_H   = facecolors(intensity = true, intensity_model = :planck, band = 1.65e-6)
        i_V   = facecolors(intensity = true, intensity_model = :planck, band = 0.55e-6)
        @test size(i_lin) == size(i_H) == size(i_V)
        @test i_lin != i_H                         # linear proxy != Planck
        @test i_H != i_V                           # and Planck depends on the band
        @test t_lin != i_lin                       # temperature != brightness
    end

    # -----------------------------------------------------------------
    @testset "plot3d" begin
        plot3d(TMAP, STAR)
        # Poly3DCollection broke hardest in the migration: unlike the 2-D PolyCollection
        # it needs an explicit pylist nesting, and a Julia Matrix of facecolors is read
        # column-major. Assert one polygon per pixel and one RGBA row per polygon.
        ax = main_ax()
        @test ncoll(ax) >= 1
        c = ax.collections[0]
        @test npaths(c) == STAR.npix
        fs = fcshape(c)
        @test length(fs) == 2 && fs[2] == 4
        @test fs[1] == STAR.npix
        finish("plot3d")
    end

    # -----------------------------------------------------------------
    @testset "plot2d_wireframe" begin
        for (nm, kw) in (("default", (;)),
                         ("no_compass", (compass = false,)),
                         ("rot_axis", (rotation_axis = true,)),
                         ("opaque", (hidden = false,)),
                         ("styled", (hidden_color = "pink", front_color = "navy",
                                     linewidth = 1.0)))
            plot2d_wireframe(STAR; kw...)
            ax = main_ax()
            # With hidden=true there are TWO collections: far side (zorder 1) then near
            # side (zorder 2). Every tessel is accounted for exactly once between them.
            want_hidden = get(kw, :hidden, true)   # NamedTuple get, not Dict
            @test ncoll(ax) == (want_hidden ? 2 : 1)
            total = sum(npaths(ax.collections[k]) for k in 0:(ncoll(ax) - 1))
            @test total == (want_hidden ? STAR.npix : NVIS)
            finish("wireframe_$nm")
        end
    end

    # -----------------------------------------------------------------
    @testset "plot2d_allepochs" begin
        stars = create_star_multiepochs(TES, PAR, [0.0, 2.5, 5.0])
        tm    = testmap(stars[1])
        for (nm, kw) in (("default", (;)),
                         ("mesh", (plotmesh = true,)),
                         ("epochs", (tepochs = [0.0, 2.5, 5.0],)),
                         ("no_compass", (compass = false,)),
                         ("colormap", (colormap = "viridis",)))
            plot2d_allepochs(tm, stars; kw...)
            axes = plt.gcf().get_axes()
            @test pylen(axes) >= 3                       # one panel per epoch
            for k in 0:2
                c = poly_collection(axes[k])
                @test c !== nothing
                if c !== nothing
                    @test npaths(c) == stars[k + 1].nquads_visible
                    @test fcshape(c) == [stars[k + 1].nquads_visible, 4]
                end
            end
            finish("allepochs_$nm")
        end
    end

    # -----------------------------------------------------------------
    @testset "plot2d_binary" begin
        SP1, SP2 = sp1(), sp2()
        cases = Dict(
            "default"     => (;),
            "intensity"   => (intensity = true,),
            "mesh"        => (plotmesh = true,),
            "graticules"  => (graticules = true, star_params1 = SP1, star_params2 = SP2),
            "rot_axis"    => (rotation_axis = true, star_params1 = SP1, star_params2 = SP2),
            "rot_arrow"   => (rotation_arrow = true, star_params1 = SP1, star_params2 = SP2),
            # NB: no inclination/position_angle here. A binary component is oriented by
            # the shared binary_frame (the ORBIT), not by its own inclination/PA, so
            # passing single-star angles would describe a star 102.6 deg away from the
            # one drawn. The decorations must take the orientation from the mesh.
            "mesh_orient" => (graticules = true, rotation_axis = true,
                              star_params1 = SP1, star_params2 = SP2),
            # intensity=true is REQUIRED: plot2d_binary only consults intensity_model
            # inside the `if intensity` branch, so without it this case was byte-identical
            # to "default" and tested nothing.
            "planck"      => (intensity = true, intensity_model = :planck, band = 1.65e-6),
            "fixed_scale" => (vmin = 15000.0, vmax = 30000.0),
            "no_colorbar" => (colorbar_on = false,),
            "axis_max"    => (axis_max = 2.0,),
            "background"  => (background = "black",),
            "colormap"    => (colormap = "viridis",),
            "title"       => (figtitle = "binary",),
        )
        for (nm, kw) in cases
            plot2d_binary(BT1, BT2, BS1, BS2, BP, TEP; kw...)
            ax = main_ax()
            # Two components, each its own PolyCollection, z-ordered by depth. Decorations
            # (graticules, rotation arrows) also land as PolyCollections but with an EMPTY
            # facecolor array — they are stroked, not filled — so select on having colours.
            filled = [c for c in (ax.collections[k] for k in 0:(ncoll(ax) - 1))
                      if pytype(c) == "PolyCollection" && fcshape(c)[1] > 0]
            @test length(filled) >= 2
            for c in filled
                fs = fcshape(c)
                @test length(fs) == 2 && fs[2] == 4
                @test fs[1] == npaths(c)
                check_visible(ax, c)
            end
            # The two components are 25300 K and 20585 K, so on a shared colour scale they
            # MUST come out visibly different. This is what catches a collapsed
            # normalisation — the failure mode where both stars render the same flat
            # colour and the figure looks empty against a white background.
            m1 = vec(sum(pyconvert(Array, filled[1].get_facecolor())[:, 1:3], dims = 2))
            m2 = vec(sum(pyconvert(Array, filled[2].get_facecolor())[:, 1:3], dims = 2))
            @test abs(sum(m1)/length(m1) - sum(m2)/length(m2)) > 0.05
            finish("binary_$nm")
        end

        # Both components keep their limb outline WHATEVER else is switched on. This is a
        # regression guard: `draw_graticules` draws its own limb at `zorder - 1`, which is
        # right for the single-star default of 5 but not for the per-component fractional
        # zorders here (zord + 0.5 - 1 lands UNDER that star's own surface). Skipping the
        # explicit limb whenever graticules were requested therefore removed the outline
        # from every binary drawn with `graticules = true`, silently — the artist was still
        # in the axes, just painted behind the star. Counting lines catches the omission;
        # only the zorder check catches the burial.
        for (nm, kw) in (("plain", (;)),
                         ("graticules", (graticules = true, star_params1 = sp1(),
                                         star_params2 = sp2())))
            plot2d_binary(BT1, BT2, BS1, BS2, BP, TEP; kw...)
            ax = main_ax()
            surf_z = sort([pyconvert(Float64, ax.collections[k].get_zorder())
                           for k in 0:(ncoll(ax) - 1)
                           if pytype(ax.collections[k]) == "PolyCollection" &&
                              fcshape(ax.collections[k])[1] > 0])
            limb_z = sort([pyconvert(Float64, ax.lines[k].get_zorder())
                           for k in 0:(nlines(ax) - 1)])
            @test length(limb_z) >= 2                    # one hull per component
            # Each limb must sit above its OWN surface: pair them up in depth order.
            @test limb_z[1] > surf_z[1] && limb_z[2] > surf_z[2]
            # ...and the far limb must stay below the near surface, or it paints the far
            # star's outline across the near one.
            @test limb_z[1] < surf_z[2]
            plt.close("all")
        end

        # Drawing into a caller-supplied axes — the path binary_movie(:compare) uses
        fig, axs = plt.subplots(1, 2, figsize = (12, 6))
        for k in 0:1
            plot2d_binary(BT1, BT2, BS1, BS2, BP, TEP;
                          ax = axs[k], colorbar_on = false, vmin = 15000.0, vmax = 30000.0)
        end
        @test pylen(fig.get_axes()) == 2                 # no stray axes created
        finish("binary_supplied_ax")
    end

    # -----------------------------------------------------------------
    @testset "plot_rv" begin
        K1, K2, γ = 120.0, 190.0, 1.0
        t = BP.T0 .+ collect(range(0, 1, length = 12)) .* BP.P
        rv1, rv2 = binary_RV(BP, t; K1 = K1, K2 = K2, γ = γ)

        fig, ax = plot_rv(BP; K1 = K1, K2 = K2, γ = γ)
        @test nlines(ax) >= 2                            # two model curves
        finish("rv_model_only")

        # 2-column data (scatter) and 3-column data (errorbar) take different branches
        fig, ax = plot_rv(BP; K1 = K1, K2 = K2, γ = γ,
                          rv_data1 = hcat(t, rv1), rv_data2 = hcat(t, rv2))
        @test ncoll(ax) >= 2                             # two scatters
        finish("rv_data_2col")

        σ = fill(3.0, length(t))
        fig, ax = plot_rv(BP; K1 = K1, K2 = K2, γ = γ,
                          rv_data1 = hcat(t, rv1, σ), rv_data2 = hcat(t, rv2, σ),
                          figtitle = "with errors")
        @test nlines(ax) >= 2
        finish("rv_data_3col")
    end

    # -----------------------------------------------------------------
    @testset "plot_mollweide" begin
        # healpix branch, every keyword. `incl` is the path where a Julia Symbol
        # (`c=:black`) can reach matplotlib, which PythonCall will not convert.
        cases = Dict(
            "default"    => (;),
            "range"      => (vmin = 3500.0, vmax = 4500.0),
            "colormap"   => (colormap = "viridis",),
            "title"      => (figtitle = "moll",),
            "incl"       => (incl = 45.0,),
            "nomask"     => (mask_unobserved = false,),
            "badcolor"   => (bad_color = "blue",
                             visible_pixels = collect(1:div(STAR.npix, 2))),
            "gridcolors" => (lon_color = "black", lat_color = "white"),
            "visible"    => (visible_pixels = collect(1:div(STAR.npix, 2)),),
            "combined"   => (incl = 30.0, vmin = 3500.0, vmax = 4500.0,
                             visible_pixels = collect(1:div(STAR.npix, 2)),
                             mask_unobserved = true, bad_color = "blue",
                             colormap = "viridis"),
        )
        for (nm, kw) in cases
            plot_mollweide(TMAP, STAR; kw...)
            @test nfigs() >= 1
            finish("mollweide_$nm")
        end

        # lat/long branch — dispatches through _latlong_dims, which must recover the grid
        # dimensions from the mesh because `stellar_geometry` does not store them.
        @test ROTIR._latlong_dims(STARLL) == (NTHETA, NPHI)
        for (nm, kw) in (("default", (;)),
                         ("explicit", (ntheta = NTHETA, nphi = NPHI)),
                         ("range", (vmin = 3500.0, vmax = 4500.0)),
                         ("incl", (incl = 45.0,)))
            plot_mollweide(TMAPLL, STARLL; kw...)
            @test nfigs() >= 1
            finish("mollweide_ll_$nm")
        end
    end

    # -----------------------------------------------------------------
    @testset "animation" begin
        t1, t2 = tessellation_healpix(2), tessellation_healpix(2)

        # Frame maps with and without the reflection step
        for refl in (false, true)
            m = binary_frame_maps(t1, RP1, t2, RP2, BP, TEP; reflection = refl)
            @test length(m) >= 2
        end

        mktempdir() do dir
            for panels in (:single, :compare)
                binary_movie(BP, t1, RP1, t2, RP2;
                             nframes = 2, panels = panels, outdir = dir,
                             prefix = "t_$(panels)", encode = false, verbose = false)
                pngs = filter(f -> startswith(f, "t_$(panels)") && endswith(f, ".png"),
                              readdir(dir))
                @test length(pngs) == 2
            end
            # ffmpeg absent must warn and return nothing, not throw
            if Sys.which("ffmpeg") === nothing
                @test frames_to_movie(dir, "t_single"; verbose = false) === nothing
            end
        end
        plt.close("all")
    end

    # -----------------------------------------------------------------
    @testset "no figure leaks" begin
        plt.close("all")
        @test nfigs() == 0
    end
end

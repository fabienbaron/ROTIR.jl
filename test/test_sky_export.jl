# Sky cubes for OITOOLS: the cube, fed to OITOOLS, must reproduce ROTIR's own visibilities.
#
# A LOPSIDED surface on purpose (a bright spot and a dark one, at non-symmetric places): a
# symmetric star has real visibilities, and a transposed or mirrored image would pass every
# check below. Closure phases are what catch a wrong layout; the complex visibilities are what
# catch a one-pixel centring error.
using Test, ROTIR, OITOOLS, Statistics, LinearAlgebra, FITSIO

const ZETD = joinpath(pkgdir(ROTIR), "demos", "zetAnd")
files = sort(filter(f -> endswith(f, ".oifits"), readdir(ZETD)))[1:2]
data = [readoifits(joinpath(ZETD, f); verbose = false)[1] for f in files]
tep = [minimum(d.v2_mjd) for d in data]
tess = tessellation_healpix(3; T = Float64)
rmin = 2.495 * 1.0054 / 2
p = (surface_type = 1, radius_x = rmin * 1.063, radius_y = rmin, radius_z = rmin, tpole = 4372.0,
     ldtype = 3, ld1 = 0.269, ld2 = 0.0, beta = 0.25, inclination = 69.3, position_angle = 124.0,
     rotation_period = 17.769426, t0 = 49992.281 + 17.769426 / 2)
lat = 90 .- rad2deg.(tess.unit_spherical[:, 5, 2]); lon = rad2deg.(tess.unit_spherical[:, 5, 3])
x = fill(4400.0, tess.npix)
x[(abs.(lat .- 20) .< 15) .& (abs.(mod.(lon .- 30, 360) .- 180) .> 150)] .= 3700.0   # dark spot
x[(abs.(lat .+ 10) .< 12) .& (abs.(mod.(lon .- 200, 360) .- 180) .> 160)] .= 4800.0  # bright spot

g = sky_grid(p, tess, tep; data = data)
@testset "sky_grid" begin
    @test iseven(g.nx) && g.nx * g.pixsize >= 6 * g.diameter * 0.99
    @test occursin("Nyquist", g.reason)
end

@testset "band-limited cube reproduces ROTIR's visibilities in OITOOLS" begin
    for model in (:linear, :planck)
        λs = model === :linear ? [1.6e-6] : [1.55e-6, 1.65e-6]
        cube, meta = render_sky_cube(p, x, tess; tepochs = tep, wavelengths = λs,
                                     pixsize = g.pixsize, nx = g.nx, intensity_model = model)
        @test size(cube) == (g.nx, g.nx, length(λs), length(tep))
        @test all(abs(sum(cube[:, :, i, j]) - 1) < 1e-4 for i in axes(cube, 3), j in axes(cube, 4))
        chk = check_sky_cube(cube, meta, x, tess, data)
        @test chk.worst.rel_dV < 5e-5             # ~3e-5 at the default 6x field
        @test chk.worst.rms_dT3phi_deg < 0.06      # closure phases: the layout is right
        @test chk.worst.chi2dist < 5e-5            # the image adds no χ² worth the name
    end
end

@testset "the layout is the one OITOOLS needs, not one of its seven neighbours" begin
    cube, meta = render_sky_cube(p, x, tess; tepochs = tep[1:1], wavelengths = [1.6e-6],
                                 pixsize = g.pixsize, nx = g.nx)
    good = check_sky_cube(cube, meta, x, tess, data[1:1]).worst
    for bad in (permutedims(cube[:, :, 1, 1]), reverse(cube[:, :, 1, 1]; dims = 1),
                circshift(cube[:, :, 1, 1], (1, 0)))
        c2 = reshape(bad, size(bad)..., 1, 1)
        @test check_sky_cube(c2, meta, x, tess, data[1:1]).worst.rel_dV > 100 * good.rel_dV
    end
end

@testset "area rendering is positive, flux-conserving, and measurably less exact" begin
    cA, mA = render_sky_cube(p, x, tess; tepochs = tep[1:1], wavelengths = [1.6e-6],
                             pixsize = g.pixsize, nx = g.nx, rendering = :area)
    @test all(>=(0), cA)
    cB, mB = render_sky_cube(p, x, tess; tepochs = tep[1:1], wavelengths = [1.6e-6],
                             pixsize = g.pixsize, nx = g.nx)
    @test any(<(0), cB)                            # ringing at the limb: part of being exact
    eA = check_sky_cube(cA, mA, x, tess, data[1:1]).worst.rel_dV
    eB = check_sky_cube(cB, mB, x, tess, data[1:1]).worst.rel_dV
    @test eA > 10eB
end

@testset "relative normalisation keeps the spectrum" begin
    cube, meta = render_sky_cube(p, x, tess; tepochs = tep[1:1], wavelengths = [1.5e-6, 2.2e-6],
                                 pixsize = g.pixsize, nx = g.nx, intensity_model = :planck,
                                 normalize = :relative)
    f = [sum(cube[:, :, i, 1]) for i in 1:2]
    @test maximum(f) ≈ 1 rtol = 1e-5
    # A 4400 K photosphere is past its peak at 1.5 µm: the flux per unit λ falls toward 2.2 µm.
    @test f[2] < f[1]
end

@testset "FITS round trip and header" begin
    cube, meta = render_sky_cube(p, x, tess; tepochs = tep, wavelengths = [1.55e-6, 1.65e-6],
                                 pixsize = g.pixsize, nx = g.nx, mjd = tep)
    path = tempname() * ".fits"
    write_sky_cube(path, cube, meta)
    r = read_sky_cube(path)
    @test r.cube == cube
    @test r.pixsize ≈ g.pixsize && r.wavelengths ≈ meta.wavelengths && r.mjd ≈ tep
    # OITOOLS reads the pixel scale back in mas from the WCS alone.
    @test fits_pixsize(path) ≈ g.pixsize rtol = 1e-9
    rm(path)
end

@testset "binaries: both components in one image reproduce ROTIR's binary model" begin
    p2 = merge(p, (radius_x = 0.45, radius_y = 0.45, radius_z = 0.45, inclination = 60.0))
    x2 = fill(6000.0, tess.npix)
    for (sep, label) in ((1.8, "close"), (12.0, "wide"))
        offs = [(sep * cosd(40 + 15k), sep * sind(40 + 15k)) for k in eachindex(tep)]   # (west, north)
        comps = [(; params = p, x, tess), (; params = p2, x = x2, tess, secondary = true, offsets = offs)]
        gb = sky_grid(comps, tep; data = data)
        @test gb.nx * gb.pixsize >= 2gb.extent           # the whole system is in the field
        for model in (:linear, :planck)
            cube, meta = render_sky_cube(comps; tepochs = tep, wavelengths = [1.6e-6],
                                         pixsize = gb.pixsize, nx = gb.nx, intensity_model = model)
            chk = check_sky_cube(cube, meta, comps, data)
            # Two limbs ring instead of one: measured 5e-5 to 1.2e-4 at the default margin
            # (Planck, with a 6000 K companion, the worse), falling to ~3e-5 at twice the
            # margin. Either way the image adds ~1e-4 to a χ² per point — nothing.
            @test chk.worst.rel_dV < 2e-4
            @test chk.worst.rms_dT3phi_deg < 0.06
            @test chk.worst.chi2dist < 2e-4
        end
        # And the file says where the companion is.
        cube, meta = render_sky_cube(comps; tepochs = tep, wavelengths = [1.6e-6],
                                     pixsize = gb.pixsize, nx = gb.nx, mjd = tep,
                                     bandwidths = [3.5e-8])
        path = tempname() * ".fits"
        write_sky_cube(path, cube, meta)
        FITSIO.FITS(path) do f
            ep = f["EPOCHS"]
            @test read(ep, "OFF2_W") ≈ first.(offs) && read(ep, "OFF2_N") ≈ last.(offs)
            @test read(f["CHANNELS"], "EFF_BAND") ≈ [3.5e-8]
            @test FITSIO.read_header(f[1])["NCOMP"] == 2
            @test any(startswith("c2."), read(f["PARAMS"], "NAME"))
        end
        rm(path)
    end
end

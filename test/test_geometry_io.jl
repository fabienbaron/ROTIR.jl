# The geometry and map FITS files: an exact round trip, and the one thing FITS cannot carry.
#
# `save_star_geometry` exists because a parameter-only file reproduces the geometry OF THE CODE
# THAT READS IT — the Roche equipotential is a root solve, the rapid rotator's radii come from
# a cubic, the visibility clip is a sigmoid with a κ. So the assertions here are about the mesh
# coming back bit-for-bit, not about it being close.

@testset "the geometry as a file" begin
    tess = tessellation_healpix(3)

    @testset "a single star round-trips exactly" begin
        p = default_star_params(2)                     # rapid rotator: an oriented mesh
        star = create_star(tess, p, 0.7)
        f = joinpath(mktempdir(), "geom.fits")
        @test save_star_geometry(f, star, p; nside_exp = 3) == f
        g = load_star_geometry(f)

        # BIT FOR BIT, on every array the struct carries. A mesh that came back "close" would
        # defeat the point of writing it: the file is the record of what the geometry WAS.
        @test g.star.vertices_xyz == star.vertices_xyz
        @test g.star.vertices_spherical == star.vertices_spherical
        @test g.star.normals == star.normals
        @test g.star.proj_west == star.proj_west
        @test g.star.proj_north == star.proj_north
        @test g.star.ldmap == star.ldmap
        @test g.star.vis_weights == star.vis_weights
        @test g.star.sig_args == star.sig_args
        @test g.star.center_offsets == star.center_offsets
        @test g.star.index_quads_visible == star.index_quads_visible
        @test g.star.nquads_visible == star.nquads_visible
        @test g.star.surface_type == star.surface_type
        @test g.star.tessellation_type == star.tessellation_type
        @test g.star.npix == star.npix
        @test g.star.t == star.t

        # The element type travels too. Reading a Float32 mesh as Float64 would be a quiet
        # lie about the precision the geometry was computed at.
        @test eltype(g.star.vertices_xyz) === eltype(star.vertices_xyz)
        # `polyft` is deliberately absent — it belongs to a dataset's uv points, not to the
        # star — and comes back empty, which is what `create_star` leaves there as well.
        @test isempty(g.star.polyft)

        @test g.params == p
        @test g.nside_exp == 3 && g.tessellation === :healpix
        # The epoch comes back as the star CARRIED it, not as it was typed: `stellar_geometry`
        # stores `t::T`, so a Float32 mesh holds 0.7 as 0.69999999 and the file records that
        # rather than inventing back the Float64 the caller passed.
        @test g.t == Float64(star.t)
        @test g.t ≈ 0.7 atol = 1e-6
        @test g.star2 === nothing && g.params2 === nothing
        @test g.tepochs === nothing && g.mjd === nothing
    end

    @testset "a binary, and the Unicode orbital elements" begin
        # A Roche model is the case that matters, twice over: it is the only surface type
        # whose parameters include the orbit, and those are spelled Ω, ω and dω — Unicode by
        # necessity, since `compute_coeff` reads the fields by those exact names. A FITS
        # string column is ASCII, so writing them verbatim throws
        # "FITS file format accepts ASCII strings only" on every binary, which is what the
        # alias table in surface_map_io.jl is for.
        p  = default_star_params(3)
        p2 = merge(default_star_params(3), (rpole = 0.3, q = 1 / p.q))
        s1 = create_star(tess, p, 0.0)
        s2 = create_star(tess, p2, 0.0; secondary = true)
        f = joinpath(mktempdir(), "binary.fits")
        save_star_geometry(f, s1, p; nside_exp = 3, star2 = s2, params2 = p2,
                           place = :offset, offset = (3.0, 1.5, 0.25),
                           tepochs = [0.0, 1.5], mjd = [55806.0, 55807.5],
                           comment = "a test binary")
        g = load_star_geometry(f)

        @test haskey(g.params, :Ω) && haskey(g.params, :ω) && haskey(g.params, :dω)
        @test g.params == p
        @test g.params2 == p2
        @test g.params2.q ≈ 1 / g.params.q          # the secondary's own convention survives
        @test g.star.vertices_xyz == s1.vertices_xyz
        @test g.star2.vertices_xyz == s2.vertices_xyz
        @test g.star2.ldmap == s2.ldmap
        @test g.place === :offset
        @test g.offset == (3.0, 1.5, 0.25)
        @test g.tepochs == [0.0, 1.5]
        @test g.mjd == [55806.0, 55807.5]

        # The same aliasing on the MAP file, which had the identical latent problem and had
        # simply never been given a Roche model to write.
        fm = joinpath(mktempdir(), "roche_map.fits")
        x = Float64.(parametric_temperature_map(p, s1))
        save_surface_map(fm, x, p; nside_exp = 3)
        mp = load_surface_map(fm)
        @test mp.params == p
        @test mp.x ≈ x
    end

    @testset "the rebuild, measured rather than assumed" begin
        # What the file is FOR: comparing a rebuild against the mesh instead of trusting it.
        p = default_star_params(0)
        star = create_star(tess, p, 0.0)
        f = joinpath(mktempdir(), "sphere.fits")
        save_star_geometry(f, star, p; nside_exp = 3)
        g = load_star_geometry(f)
        r = create_star(tessellation_healpix(g.nside_exp), g.params, g.t)
        # A sphere is a multiplication, so it rebuilds exactly. A Roche surface is a root
        # solve and lands within Float32 round-off of the stored mesh — which is the number
        # the GUI reports on load, and the reason it reports one at all.
        @test maximum(abs.(r.vertices_xyz .- g.star.vertices_xyz)) == 0
    end

    @testset "what it refuses" begin
        p = default_star_params(0)
        star = create_star(tess, p, 0.0)
        d = mktempdir()
        @test_throws ArgumentError save_star_geometry(joinpath(d, "a.fits"), star, p;
                                                      nside_exp = 3,
                                                      tessellation = :nosuchgrid)
        @test_throws ArgumentError save_star_geometry(joinpath(d, "b.fits"), star, p;
                                                      nside_exp = 3, place = :somehow)
        # Half a binary is not a binary: a companion mesh with no parameters could be loaded
        # but never rebuilt, and the failure would come much later.
        @test_throws ArgumentError save_star_geometry(joinpath(d, "c.fits"), star, p;
                                                      nside_exp = 3, star2 = star)
    end
end

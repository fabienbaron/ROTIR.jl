@compile_workload begin
    # Exercise the core computational path with a small tessellation (Float32 default).
    # No plotting: every drawing function is a stub here, with methods only in
    # ROTIRPythonPlotExt or ROTIRMakieExt, and those extensions carry their own workloads.
    tessels = tessellation_healpix(1)
    star_params = (
        surface_type    = 0,
        radius          = 1.0f0,
        tpole           = 10000.0f0,
        ldtype          = 3,
        ld1             = 0.3f0,
        ld2             = 0.0f0,
        inclination     = 35.0f0,
        position_angle  = 20.0f0,
        rotation_period = 1.0f0
    )
    star = create_star(tessels, star_params, 0.0f0)
    tmap = parametric_temperature_map(star_params, star)

    # Tessellation utilities
    tessellation_healpix(2)
    tessellation_latlong(8, 16)

    # THE VISIBILITY PATH, which used to be absent from this workload entirely — the geometry
    # and the temperature map were exercised and then nothing that turns them into data. The
    # default forward kernel is now `:t3` (src/type3_nufft.jl), which is heavily specialised:
    # `Type3Plan{T,FP,W}` carries the FFTW plan type and the kernel width as type parameters,
    # and the stencils are `Val`-dispatched NTuples. That is a lot to compile on first call.
    #
    # SYNTHETIC uv COORDINATES rather than a file: this workload reads nothing from disk, and
    # what is being compiled does not depend on the numbers being real. The scale is chosen to
    # look like an interferometric dataset in radians-per-mas so the plan sizing takes its
    # normal branch rather than the degenerate one.
    let idx = star.index_quads_visible,
        pw = Matrix(star.proj_west[idx, :]), pn = Matrix(star.proj_north[idx, :]),
        xw = star.vis_weights[idx] .* star.ldmap[idx],
        kx = collect(range(-0.9f0, 0.9f0, length = 12)),
        ky = collect(range(0.7f0, -0.8f0, length = 12))

        ng, ns = quadrature_for_type3(pw, pn, kx, ky)
        en, ew = t3_gauss_rule(ng, ns, Float64)
        F = type3_cvis(pw, pn, xw, kx, ky)
        g = Vector{Float32}(undef, length(xw))
        adj = Complex{Float32}.(range(-1.0f0, 1.0f0, length = 12),
                                range(0.5f0, -0.5f0, length = 12))
        type3_cvis_adj!(g, pw, pn, adj, kx, ky)

        # The plan-level entry points too, in Float64, which is what the wrappers above build
        # internally — and the generic point-set route, which is a separate specialisation.
        p64 = plan_type3(1.0, Float64.(kx), Float64.(ky))
        out64 = Vector{ComplexF64}(undef, length(kx))
        type3_quads!(out64, p64, Float64.(pw), Float64.(pn), Float64.(xw), en, ew)
        g64 = Vector{Float64}(undef, length(xw))
        type3_quads_adj!(g64, p64, Float64.(pw), Float64.(pn), ComplexF64.(adj), en, ew)
        type3_points!(out64, p64, [0.3, -0.4], [0.2, 0.5], [1.0, 0.7])

        # And the exact kernel the other backends use, so a cross-check is not a cold start.
        nuv = length(kx)
        Fs = Vector{Complex{Float32}}(undef, nuv); pf = zeros(Float32, length(xw))
        k2 = precompute_k2_inv_im(kx, ky)
        compute_polyflux_and_cvis!(Fs, pf, kx, ky, k2, pw, pn, xw)
        gs = Vector{Float32}(undef, length(xw))
        compute_adjoint_cvis!(gs, adj, kx, ky, k2, pw, pn, pf)
    end
    # CLEARED, and this is not tidiness. A `Type3Plan` holds an FFTW plan — a C pointer — and
    # global state mutated inside `@compile_workload` is serialised into the precompiled
    # image, so leaving a plan here would hand every later session a dangling pointer.
    empty!(TYPE3_PLANS)
end

# The geometry builders in both element types, as explicit hints.
#
# The matplotlib hints that used to sit here — `plot2d`, `draw_compass(::Py, …)` and the rest —
# are gone. They named functions that are now STUBS in this module, so `precompile` had nothing
# to compile and returned false silently; and `Py` is no longer a name this package has, since
# PythonCall became a weak dependency.
# Primary: Float32 (default)
let T = Float32, NT = @NamedTuple{surface_type::Int, radius::Float32, tpole::Float32,
        ldtype::Int, ld1::Float32, ld2::Float32,
        inclination::Float32, position_angle::Float32, rotation_period::Float32}
    precompile(create_star, (tessellation{T}, NT, T))
    precompile(parametric_temperature_map, (NT, stellar_geometry{T}))
end
# The type-3 plan and its two directions, as explicit hints for the width the wrappers use.
# `W` is a type parameter, so each width is a separate specialisation; 11 is the default and
# the only one anything reaches unless a caller asks otherwise.
precompile(quadrature_for_type3, (Matrix{Float32}, Matrix{Float32}, Vector{Float32}, Vector{Float32}))
precompile(quadrature_for_type3, (Matrix{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}))
precompile(t3_gauss_rule, (Int, Int, Type{Float64}))
precompile(type3_cvis, (Matrix{Float32}, Matrix{Float32}, Vector{Float32}, Vector{Float32}, Vector{Float32}))
precompile(type3_cvis, (Matrix{Float64}, Matrix{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}))
precompile(type3_cvis_adj!, (Vector{Float32}, Matrix{Float32}, Matrix{Float32},
                             Vector{Complex{Float32}}, Vector{Float32}, Vector{Float32}))
precompile(type3_cvis_adj!, (Vector{Float64}, Matrix{Float64}, Matrix{Float64},
                             Vector{ComplexF64}, Vector{Float64}, Vector{Float64}))

# Secondary: Float64 opt-in
let T = Float64, NT = @NamedTuple{surface_type::Int, radius::Float64, tpole::Float64,
        ldtype::Int, ld1::Float64, ld2::Float64,
        inclination::Float64, position_angle::Float64, rotation_period::Float64}
    precompile(create_star, (tessellation{T}, NT, T))
    precompile(parametric_temperature_map, (NT, stellar_geometry{T}))
end

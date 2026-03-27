@compile_workload begin
    # Exercise the core computational path with a small tessellation (Float32 default).
    # Plotting is excluded — PyCall/matplotlib segfaults during precompilation.
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
end

# Explicit precompile hints for plotting functions.
# These can't run in @compile_workload (PyCall segfaults during precompilation)
# but the hints still cache inference/compilation results.
# Primary: Float32 (default)
let T = Float32, NT = @NamedTuple{surface_type::Int, radius::Float32, tpole::Float32,
        ldtype::Int, ld1::Float32, ld2::Float32,
        inclination::Float32, position_angle::Float32, rotation_period::Float32}
    precompile(plot2d, (Vector{T}, stellar_geometry{T}))
    precompile(plot2d_wireframe, (stellar_geometry{T},))
    precompile(draw_compass, (PyCall.PyObject, T))
    precompile(draw_rotation_axis, (PyCall.PyObject, stellar_geometry{T}))
    precompile(draw_rotation_arrow, (PyCall.PyObject, stellar_geometry{T}))
    precompile(draw_graticules, (PyCall.PyObject, stellar_geometry{T}))
    precompile(create_star, (tessellation{T}, NT, T))
    precompile(parametric_temperature_map, (NT, stellar_geometry{T}))
end
# Secondary: Float64 opt-in
let T = Float64, NT = @NamedTuple{surface_type::Int, radius::Float64, tpole::Float64,
        ldtype::Int, ld1::Float64, ld2::Float64,
        inclination::Float64, position_angle::Float64, rotation_period::Float64}
    precompile(create_star, (tessellation{T}, NT, T))
    precompile(parametric_temperature_map, (NT, stellar_geometry{T}))
end

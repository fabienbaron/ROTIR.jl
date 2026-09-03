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
# Secondary: Float64 opt-in
let T = Float64, NT = @NamedTuple{surface_type::Int, radius::Float64, tpole::Float64,
        ldtype::Int, ld1::Float64, ld2::Float64,
        inclination::Float64, position_angle::Float64, rotation_period::Float64}
    precompile(create_star, (tessellation{T}, NT, T))
    precompile(parametric_temperature_map, (NT, stellar_geometry{T}))
end

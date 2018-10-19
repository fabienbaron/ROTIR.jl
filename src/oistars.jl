mutable struct starparameters
  radius::Float64 # milliarcseconds (at pole)
  temperature::Float64 # Kelvin (at pole)
  frac_escapevel::Float64 # unitless; fractional rotational velocity

  ld::Array{Float64,1} # limb darkening,first coefficient is for LD law type, then LD coefficients
  beta::Float64 # exponent for von Zeipel law
  B_rot::Float64 # 2nd constant for rotational velocity
  #C_rot::Float64 # 3rd constant for rotational velocity -- not in use anymore! Probably only used for the Sun

  inclination::Float64 # degrees
  position_angle::Float64 # degrees
  selfrotangle::Float64 # degrees
  rotation_period::Float64 # days
end

mutable struct binaryparameters
  star1::starparameters
  star2::starparameters
  
  # parameters for single star or binary
  orbit_incl::Float64 # degrees
  long_ascending_node::Float64 # degrees (longitude of ascending node for binaries)
  arg_pariapsis::Float64 # degrees (argument of periapsis for binaries)
  binary_period::Float64 # days (period for binaries)

  true_anomaly::Float64 # degrees
  separation::Float64 # in milliarcseconds (aka semi-major axis)
  eccentricity::Float64 # unitless
  mass_ratio::Float64 # unitless; M2/M1
  async_ratio::Float64 # unitless; Ratio self-rotation/orbital revolution period
  fillout_factor::Float64 # unitless; value of the potential at Roche lobe divided by value of potential at the surface
end
mutable struct starparameters{T}
  rpole::T # milliarcseconds (at pole)
  tpole::T # Kelvin (at pole)
  frac_escapevel::T # unitless; fractional rotational velocity
  ldtype::Int64
  ld1::T # limb darkening,first coefficient is for LD law type, then LD coefficients
  ld2::T
  beta_vZ::T # exponent for von Zeipel law
  B_rot::T # 2nd constant for rotational velocity
  inclination::T # degrees
  position_angle::T # degrees
  rotation_offset::T # degrees # offsets the rotation by a fixed angle
  rotation_period::T # Rotation period days
  #T_long0::T # Reference time for longitude=0
end

mutable struct binaryparameters{T}
  star1::starparameters
  star2::starparameters
  d::T # distance (parsecs)
  # parameters for binary
  i::T # degrees
  Ω::T # degrees
  ω::T # degrees
  P::T # days
  a::T # in milliarcseconds (aka semi-major axis)
  e::T # unitless
  T0::T # JD; time of periastron
  q::T # unitless; M2/M1
  fillout_factor::Array{T,1} # unitless; value of the potential at Roche lobe divided by value of potential at the surface
  dP::T # days/day; linear change of the binary period
  dω::T # periapsis change (degrees/day)
end

Base.iterate(bparameters::binaryparameters, state = 1) = state <= 1 ? (bparameters, state+1) : nothing
Base.length(bparameters::binaryparameters) = 1
 


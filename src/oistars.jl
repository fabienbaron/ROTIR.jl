"""
    starparameters(rpole, tpole, frac_escapevel, ldtype, ld1, ld2, beta_vZ, B_rot,
                   inclination, position_angle, rotation_offset, rotation_period)

Construct a stellar parameters NamedTuple. Fields:
- `rpole`: polar radius (mas)
- `tpole`: polar temperature (K)
- `frac_escapevel`: fractional rotational velocity (0–1)
- `ldtype`: limb-darkening law (1=linear, 2=quadratic, 3=Hestroffer)
- `ld1`, `ld2`: limb-darkening coefficients
- `beta_vZ`: von Zeipel gravity-darkening exponent (0.25 radiative, 0.08 convective)
- `B_rot`: differential rotation coefficient
- `inclination`: spin-axis inclination (degrees)
- `position_angle`: spin-axis PA (degrees)
- `rotation_offset`: fixed rotation offset (degrees)
- `rotation_period`: rotation period (days)
"""
function starparameters(rpole, tpole, frac_escapevel, ldtype, ld1, ld2, beta_vZ, B_rot,
                        inclination, position_angle, rotation_offset, rotation_period)
    return (rpole=rpole, tpole=tpole, frac_escapevel=frac_escapevel,
            ldtype=ldtype, ld1=ld1, ld2=ld2, beta_vZ=beta_vZ, B_rot=B_rot,
            inclination=inclination, position_angle=position_angle,
            rotation_offset=rotation_offset, rotation_period=rotation_period)
end

"""
    binaryparameters(star1, star2, d, i, Ω, ω, P, a, e, T0, q, fillout_factor, dP, dω)

Construct a binary parameters NamedTuple. Fields:
- `star1`, `star2`: NamedTuple from `starparameters()` for each component
- `d`: distance (pc)
- `i`: orbital inclination (degrees)
- `Ω`: longitude of ascending node (degrees)
- `ω`: argument of periapsis (degrees)
- `P`: orbital period (days)
- `a`: semi-major axis (mas)
- `e`: eccentricity
- `T0`: time of periastron passage (JD)
- `q`: mass ratio M₂/M₁
- `fillout_factor`: Roche lobe fillout factors [primary, secondary]
- `dP`: period derivative (days/day)
- `dω`: periapsis precession (degrees/day)
"""
function binaryparameters(star1, star2, d, i, Ω, ω, P, a, e, T0, q, fillout_factor, dP, dω)
    return (star1=star1, star2=star2, d=d, i=i, Ω=Ω, ω=ω, P=P, a=a, e=e,
            T0=T0, q=q, fillout_factor=fillout_factor, dP=dP, dω=dω)
end

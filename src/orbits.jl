"""
    kepler_E(M, e) -> E

Solve Kepler's equation `M = E − e·sin E` for the eccentric anomaly of an elliptical orbit
(`0 ≤ e < 1`), scalar in / scalar out. The result is wrapped to `[0, 2π)`.

Converged to machine precision by a Newton iteration that is *provably* monotone:

* `M` is wrapped to `[0, 2π)` and, if it exceeds π, reflected to `2π − M` (the equation is
  odd about π), so the solve always runs on `M ∈ [0, π]`.
* There `f(E) = E − e·sin E − M` is increasing with `f'' = e·sin E ≥ 0`, i.e. convex, and
  `f(0) ≤ 0 ≤ f(π)`. Newton started from any point with `f ≥ 0` therefore decreases
  monotonically onto the root with no bracketing or damping.
* The start `E₀ = min(π, M/(1−e))` satisfies `f(E₀) ≥ 0` because `sin E ≤ E` on `[0, π]`,
  and is much tighter than π alone when `M` is small — typically 3–5 iterations.

# Differentiability

`kepler_E` carries an analytic `ChainRulesCore.rrule` (and `frule`) obtained by
differentiating Kepler's equation implicitly rather than by unrolling the solver:

```
M = E − e·sin E   ⇒   dM = (1 − e·cos E)·dE − sin E·de
⇒   ∂E/∂M = 1/(1 − e·cos E),   ∂E/∂e = sin E/(1 − e·cos E)
```

Both derivatives are *exact* at the converged root and cost two trig evaluations, whereas
reverse-mode AD through the iteration would tape every Newton step, allocate per step, and
return a derivative only as accurate as the convergence tolerance. `1 − e·cos E = r/a > 0`
for `e < 1`, so the denominator never vanishes on an ellipse.

See also [`compute_eccentric_anomaly`](@ref), [`compute_true_anomaly`](@ref).
"""
function kepler_E(M::Real, e::Real)
    T = float(promote_type(typeof(M), typeof(e)))
    ee = T(e)
    (ee < 0 || ee >= 1) &&
        throw(DomainError(e, "kepler_E: elliptical orbits only, need 0 ≤ e < 1"))
    Mw = mod2pi(T(M))
    ee == 0 && return Mw
    flip = Mw > T(π)
    Mr   = flip ? T(2π) - Mw : Mw            # now Mr ∈ [0, π]
    E    = min(T(π), Mr / (1 - ee))          # f(E₀) ≥ 0, so Newton descends monotonically
    for _ in 1:100
        δ = (E - ee*sin(E) - Mr) / (1 - ee*cos(E))
        E -= δ
        abs(δ) <= 4*eps(T)*max(one(T), abs(E)) && break
    end
    return flip ? T(2π) - E : E
end

function ChainRulesCore.rrule(::typeof(kepler_E), M::Real, e::Real)
    E = kepler_E(M, e)
    function kepler_E_pullback(Ēraw)
        Ē   = unthunk(Ēraw)
        den = 1 - e*cos(E)                   # = r/a, strictly positive for e < 1
        return (NoTangent(), Ē / den, Ē * sin(E) / den)
    end
    return E, kepler_E_pullback
end

function ChainRulesCore.frule((_, ΔM, Δe), ::typeof(kepler_E), M::Real, e::Real)
    E   = kepler_E(M, e)
    den = 1 - e*cos(E)
    return E, (ΔM + sin(E)*Δe) / den
end

"""
    kepler_E_vec(M::AbstractArray, e) -> E

Array form of [`kepler_E`](@ref), with an **array-level** `rrule`.

Broadcasting the scalar `kepler_E` works and is correct, but reverse-mode AD then invokes
the scalar pullback once per element and allocates around each one; on a 55-epoch orbit that
made the orbit stage cost 37× its primal, against ~1.6× for a hand-written adjoint over the
same data. Differentiating the whole array in one pullback removes that: the derivative is
still exactly the implicit-function result, just applied as three vector operations.

`e` is shared by every element, so its cotangent is the *sum* over the array.
"""
kepler_E_vec(M::AbstractArray, e) = kepler_E.(M, e)

function ChainRulesCore.rrule(::typeof(kepler_E_vec), M::AbstractArray, e::Real)
    E = kepler_E_vec(M, e)
    function kepler_E_vec_pullback(Ēraw)
        Ē   = unthunk(Ēraw)
        den = 1 .- e .* cos.(E)              # = r/a, strictly positive for e < 1
        return (NoTangent(), Ē ./ den, sum(Ē .* sin.(E) ./ den))
    end
    return E, kepler_E_vec_pullback
end

"""
    compute_E_NR(M, e) -> E

Eccentric anomaly from the mean anomaly, elementwise over `M` (scalar or vector). The
result takes the float type of `M` and `e`.

!!! note "The `T` keyword is accepted and ignored"
    The working precision is derived from `M` and `e`, not dictated by a keyword. `T`
    survives only so existing call sites keep working.


Dispatches to [`kepler_E`](@ref) for a scalar and [`kepler_E_vec`](@ref) for an array — the
latter so that reverse-mode AD gets one array pullback rather than one per element. Kept as
a separate name because it is the historical entry point.

!!! note "Converges to machine precision"
    Results are not bitwise comparable with codes that stop at a fixed `1e-6` rad tolerance
    on the whole vector; expect differences up to ~1e-6 rad of eccentric anomaly, far below
    any measurement error.
"""
# `T` is accepted and ignored: kepler_E derives its working type from M and e, which is
# strictly better than being told. Kept in the signature so existing callers do not break.
compute_E_NR(M::Real, e; T = nothing)          = kepler_E(M, e)
compute_E_NR(M::AbstractArray, e; T = nothing) = kepler_E_vec(M, e)

"""
    omega_at(bparameters, tepoch) -> ω(t)  [radians]

Argument of periapsis at `tepoch`, including apsidal motion:

    ω(t) = ω₀ + dω·(t − T0)

with `bparameters.dω` in **degrees per day** and `tepoch`, `T0` in JD. Returns `ω₀` when
`dω = 0`, so circular / non-precessing systems are unaffected.

Apsidal motion is not a small correction for a system observed over years: Spica's apsidal period is U ≈ 139 yr
(Robinette & Aufdenberg 2015; 105 yr in Wages & Aufdenberg, ~110 yr in Aufdenberg+2015),
i.e. dω/dt ≈ 2.6°/yr, so ω advances ~21° across the 2007–2015 CHARA campaigns. Ignoring
it displaces the predicted secondary position by up to 0.44 mas — ten to forty times the
MIRC astrometric precision.

Sign convention: positive `dω` advances the periapsis in the direction of orbital motion,
matching the standard definition of the apsidal period.

Accepts a scalar or a vector of epochs. With `dω = 0` it returns the scalar `ω₀`
regardless, so downstream broadcasting is unchanged for non-precessing systems.
"""
@inline function omega_at(bparameters, tepoch)
    ω0 = bparameters.ω * pi / 180.0
    dω = hasproperty(bparameters, :dω) ? bparameters.dω : 0.0
    dω == 0 && return ω0
    return ω0 .+ (dω * pi / 180.0) .* (tepoch .- bparameters.T0)
end

# Compute coefficients for orbital equations
function compute_coeff(Omega, inclination, omega)
    L1 = cos(Omega)*cos(omega) - sin(Omega)*sin(omega)*cos(inclination);
    M1 = sin(Omega)*cos(omega) + cos(Omega)*sin(omega)*cos(inclination);
    N1 = sin(omega)*sin(inclination);
    L2 = -cos(Omega)*sin(omega) - sin(Omega)*cos(omega)*cos(inclination);
    M2 = -sin(Omega)*sin(omega) + cos(Omega)*cos(omega)*cos(inclination);
    N2 = cos(omega)*sin(inclination);
    return L1, M1, N1, L2, M2, N2
end

# Solve relative binary orbit
function binary_orbit_rel(bparameters,tepoch::Float64)
    Ω = bparameters.Ω*pi/180.; # longitude of ascending node
    i = bparameters.i*pi/180.;
    ω = omega_at(bparameters, tepoch); # argument of periapsis, advanced by apsidal motion
    a = bparameters.a;
    e = bparameters.e;
    # Eccentric anomaly
    E = compute_eccentric_anomaly(bparameters, tepoch);
    # compute orbital coefficients
    L1, M1, N1, L2, M2, N2 = compute_coeff(Ω, i, ω);
    x, y, z = compute_xyz_rel(a, sqrt(1.0 - e^2), e, L1, M1, N1, L2, M2, N2, cos(E), sin(E));
    return 0.0, 0.0, 0.0, x, y, z
end

function binary_orbit_rel_alt(bparameters,tepoch::Float64)
    Ω = bparameters.Ω*pi/180.; # longitude of ascending node
    i = bparameters.i*pi/180.;
    ω = omega_at(bparameters, tepoch); # argument of periapsis, advanced by apsidal motion
    υ = compute_true_anomaly(bparameters, tepoch);
    x = cos(Ω)*cos(ω+υ) - sin(Ω)*sin(ω+υ)*cos(i);
    y = sin(Ω)*cos(ω+υ) - cos(Ω)*sin(ω+υ)*cos(i);
    z = sin(ω+υ)*sin(i);
    return 0.0, 0.0, 0.0, x, y, z
end

"""
    binary_orbit_abs(bparameters, tepoch) -> (x1, y1, z1, x2, y2, z2)

Compute absolute positions of both stars in the observer's sky frame.

Convention: `bparameters.ω` is the argument of periapsis of the **relative orbit**
(= secondary's orbit around the primary), the standard astrometric convention.
At periastron (υ=0), the secondary is in the ω-direction from the center of mass.
"""
function binary_orbit_abs(bparameters,tepoch::Float64)
    Ω = bparameters.Ω*pi/180.0; # longitude of ascending node
    i = bparameters.i*pi/180.0;
    ω = omega_at(bparameters, tepoch); # argument of periapsis, advanced by apsidal motion
    q = bparameters.q;
    a = bparameters.a;
    e = bparameters.e;
    υ = compute_true_anomaly(bparameters, tepoch);
    D = a*(1.0 - e^2)./(1.0 .+ e* cos.(υ));
    # distance of objects from the center of mass
    r1 = D / (1/q+1);
    r2 = D / (1+q);
    L1, M1, N1, L2, M2, N2 = compute_coeff(Ω, i, ω);
    # Secondary at true anomaly υ (its periapsis is at υ=0),
    # primary on the opposite side (υ+π)
    x2, y2, z2 = compute_xyz_abs(L1, M1, N1, L2, M2, N2, υ, r2);
    if ((υ .>= 0.0) & (υ .<= pi))
        x1, y1, z1 = compute_xyz_abs(L1, M1, N1, L2, M2, N2, υ .+pi, r1);
    else
        x1, y1, z1 = compute_xyz_abs(L1, M1, N1, L2, M2, N2, υ .-pi, r1);
    end
    return x1, y1, z1, x2, y2, z2
end

# =====================================================================
# Orbital Coordinate Convention
# =====================================================================
# compute_xyz_abs/rel return positions in the observer's frame:
#   x = North (Dec increasing)
#   y = East  (RA increasing)
#   z = -toward_observer (positive = receding)
#
# To convert to ROTIR's projected frame:
#   proj_west  = -y_orbit   (West = -East)
#   proj_north =  x_orbit   (North = North)
# See orbit_to_rotir_offset() in oichi2_binary.jl.
# =====================================================================

"""
    compute_xyz_abs(L1, M1, N1, L2, M2, N2, nu, r_star) -> (north, east, z)

Returns absolute position in the observer's sky frame:
- x (north): declination direction
- y (east): right ascension direction
- z: line-of-sight (positive = receding from observer)
"""
function compute_xyz_abs(L1, M1, N1, L2, M2, N2, nu, r_star)
    astro_north = (L1.*r_star.*cos.(nu) .+ L2.*r_star.*sin.(nu));
    astro_east  = (M1.*r_star.*cos.(nu) .+ M2.*r_star.*sin.(nu));
    astro_z     = (N1.*r_star.*cos.(nu) .+  N2.*r_star.*sin.(nu));
    x = astro_north;
    y = astro_east;
    z = -astro_z;
    return x, y, z
end

# Compute the x,y,z relative positions
function compute_xyz_rel(a, β, e, L1, M1, N1, L2, M2, N2, cos_E, sin_E)
    astro_north = a*(L1*cos_E + β*L2*sin_E - e*L1);
    astro_east  = a*(M1*cos_E + β*M2*sin_E - e*M1);
    astro_z     = a*(N1*cos_E + β*N2*sin_E - e*N1);
    x = astro_north;
    y = astro_east;
    z = -astro_z;
    return x, y, z
end

function compute_separation_alt(bparameters, tepoch) # uses true anomaly
  # dimentionless instantaneous separation of the centers of mass of the two stars
  # Multiply by a to find the real separation
  υ = compute_true_anomaly(bparameters, tepoch);
  e = bparameters.e;
  D = (1 - e^2)./(1 .+ e*cos.(υ)); 
  return D
end

function compute_separation(bparameters, tepoch) 
    # dimentionless instantaneous separation of the centers of mass of the two stars
    # Multiply by a to find the real separation
    E = compute_eccentric_anomaly(bparameters, tepoch);
    e = bparameters.e;
    D = 1 .- e*cos.(E)
    return D
end

# Calculate eccentric anomaly
"""
    compute_eccentric_anomaly(bparameters, tepoch; T=<derived>) -> E

Eccentric anomaly at `tepoch` (JD; scalar or vector).

`T` defaults to the float type promoted from `bparameters.P`, `.T0`, `.e` and `eltype(tepoch)`
rather than to `Float64`, so a Float32 parameter set stays Float32 through the orbit. Note
that `tepoch` is an ABSOLUTE JD and should stay `Float64` — see `binary_orbit_abs`, which
annotates it as such: `eps(Float32(2.45e6))` is 0.25 d.


`bparameters.dP` is the period derivative Ṗ in **days per day**, applied through the
standard quadratic ephemeris — the nth periastron falls at

    t_n = T0 + P·n + ½·(Ṗ·P)·n²

whose exact inversion is the `M` below (verified against solving the quadratic directly to
1e-10). Systems whose period drifts are common in interacting binaries: β Lyrae's is
*lengthening* at Ṗ ≈ +19 s/yr ≈ 6e-7 d/d from mass transfer, while orbits shrinking to
gravitational-wave or magnetic-braking losses have Ṗ < 0.

Both signs are handled: the guard is on `|Ṗ|`, since a one-sided `Ṗ > tol` test drops
*negative* Ṗ back onto a constant period — worth 0.6° of orbital phase after 1000 d and 65°
after 10⁴ d at β Lyr's rate. Note the sign convention: Ṗ > 0 is a LENGTHENING period.
"""
# `T` defaults to the float type the INPUTS already carry rather than to Float64. Hardcoding
# Float64 here silently widened every downstream quantity, so a Float32 parameter set came
# back Float64 from orbit_to_rotir_offset — which defeats carrying Float32 through the
# forward model at all. Float64 inputs are unaffected (the promotion returns Float64).
function compute_eccentric_anomaly(bparameters, tepoch;
                                   T = float(promote_type(typeof(bparameters.P),
                                                          typeof(bparameters.T0),
                                                          typeof(bparameters.e),
                                                          eltype(tepoch))))
    P = bparameters.P;
    dP = bparameters.dP;
    e = bparameters.e;
    T0 = bparameters.T0;
    # mean angular velocity (mean motion)
    n = T(2pi)/P;
    # Mean anomaly
    if abs(dP) > 1e-12
        # Quadratic ephemeris, either sign of Ṗ. The radicand is 1 + 2Ṗ(t−T0)/P; it can
        # only go negative if the extrapolated period would pass through zero, which is
        # unphysical and means dP or the epoch baseline is wrong.
        radicand = 1 .+ 2dP*(tepoch .- T0)/P
        if any(radicand .< 0)
            @warn "compute_eccentric_anomaly: quadratic ephemeris breaks down (period " *
                  "would reach zero) for dP = $dP over this epoch range; " *
                  "falling back to a constant period."
            M = n*(tepoch .- T0)
        else
            M = T(4pi)*(tepoch .- T0) .* (1 ./ (1 .+ sqrt.(radicand)))/P;
        end
    else
        M = n*(tepoch .- T0);
    end
    # Eccentric anomaly using Newton-Raphson method
    E = mod2pi.(compute_E_NR(M,e));
    return E
end


# Calculate true anomaly
function compute_true_anomaly(bparameters,tepoch)
    E = compute_eccentric_anomaly(bparameters,tepoch)
    e = bparameters.e;
    # calculate true anomaly, following Broucke, R. ; Cefola, P. 
    # A Note on the Relations between True and Eccentric Anomalies in the Two-Body Problem
    # Supposedly avoids numerical issues when the arguments are near ± π, as the two tangents
    # in the classic equation become infinite. 
     β = e/(1+sqrt(1-e^2))
     υ = E .+ 2*atan.((β*sin.(E)) ./ (1 .-  β*cos.(E)))
    return υ
end

function compute_masses(bparameters)
    G = 6.67408e-11; # m^3/kg/s^2
    M_sun = 1.98847e30; # kg
    rad_mas = 180/pi*3600*1000
    pc = 3.08567782e16
    q = bparameters.q;
    a = bparameters.a/rad_mas; # radians
    P = bparameters.P*86400.0; # seconds
    d = bparameters.d*pc; # km
    m1 = 4*pi^2*(a*d)^3/((G*P^2*(1+q)))/M_sun; # in solar mass
    m2 = m1*q
    return m1, m2
end

"""
    binary_RV(bparameters, tepoch; K1, K2, γ) -> (Vrad1, Vrad2)

Compute radial velocities for both stars of a binary system.
Positive = receding (redshift), following the standard spectroscopic convention.

Convention: `bparameters.ω` is the argument of periapsis of the **relative orbit**
(= secondary's). The primary's argument of periapsis is ω + π.
"""
function binary_RV(bparameters, tepoch::Union{Float64, Vector{Float64}}; K1::Float64, K2::Float64, γ::Float64)
    ω = omega_at(bparameters, tepoch); # argument of periapsis, advanced by apsidal motion
    e = bparameters.e;
    υ = compute_true_anomaly(bparameters, tepoch)
    # Primary uses ω+π (its own argument of periapsis)
    Vrad1 = γ .+ K1 .* (cos.(υ .+ ω .+ pi) .+ e .* cos.(ω .+ pi))
    # Secondary uses ω directly
    Vrad2 = γ .+ K2 .* (cos.(υ .+ ω) .+ e .* cos.(ω))
    return Vrad1, Vrad2
end

function binary_proj_plane(bparameters, tepochs) 
    # Project orbit into the x,y observer plane (= screen plane)
    # Primary is at (0,0)
    Ω = bparameters.Ω*pi/180.0; # longitude of ascending node
    i = bparameters.i*pi/180.0;
    ω = omega_at(bparameters, tepochs); # argument of periapsis, advanced by apsidal motion
    a = bparameters.a;
    e = bparameters.e;
    υ = compute_true_anomaly(bparameters, tepochs)
    r = a * (1.0 - e^2)./(1.0 .+ e*cos.(υ)); # Separation between the two stars
    O_w = atan.(sin.(υ .+ ω) .* cos(i) , cos.(υ .+ ω))
    θ = O_w .+ Ω;                            # Principal axis angle in cyclindrical coordinates
    ρ = abs.(r .* cos.(υ .+ ω) ./ cos.(O_w)) # Radius in the projected cylindrical coordinates
    x = ρ .* cos.(θ)
    y = ρ .* sin.(θ)
    return x, y, ρ, θ
end

# Solve Kepler's equation using Newton-Raphson
function compute_E_NR(M,e; T=Float64) # Union{Vector{Float},Float}
    # initial guess (from Smith (1979))
    E = M + e*sin.(M)./(1 .- sin.(M .+ e) + sin.(M));
    E_old = M*0
    # Thresholds
    orbit_thresh = T(1f-6);
    orbit_max_iter = 50;

    # Newton-Raphson method
    i = 0;
    if (e > 0.5)
        while((i < orbit_max_iter) & (norm(E - E_old) > orbit_thresh))
            i += 1;
            E_old = deepcopy(E);
            E = E_old + (M + e*sin.(E_old) - E_old)./(1 .- e*cos.(E_old));

        end
    else
        while((i < orbit_max_iter) & (norm(E - E_old) > orbit_thresh))
            i += 1;
            E_old = copy(E);
            E = M + e*sin.(E_old);
        end
    end
    return E
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
    ω = bparameters.ω*pi/180.; # argument of periapsis
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
    ω = bparameters.ω*pi/180.; # argument of periapsis
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
    ω = bparameters.ω*pi/180.0; # argument of periapsis (of relative orbit / secondary)
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
function compute_eccentric_anomaly(bparameters,tepoch; T=Float64) # Note: this can be called with tepochs = vector
    P = bparameters.P;
    dP = bparameters.dP;
    e = bparameters.e;
    T0 = bparameters.T0;
    # mean angular velocity (mean motion)
    n = T(2pi)/P;
    # Mean anomaly
    if (dP > 1e-12)
        # if orbit is decaying
        M = T(4pi)*(tepoch .- T0) .* (1 ./ (1 .+ sqrt.(1 .+ 2dP*(tepoch .- T0)/P)))/P;
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
    ω = bparameters.ω*pi/180.0; # argument of periapsis (of relative orbit / secondary)
    e = bparameters.e;
    υ = compute_true_anomaly(bparameters, tepoch)
    # Primary uses ω+π (its own argument of periapsis)
    Vrad1 = γ + K1 * (cos.(υ .+ (ω .+ pi)) .+ e * cos(ω .+ pi))
    # Secondary uses ω directly
    Vrad2 = γ + K2 * (cos.(υ .+ ω) .+ e * cos(ω))
    return Vrad1, Vrad2
end

function binary_proj_plane(bparameters, tepochs) 
    # Project orbit into the x,y observer plane (= screen plane)
    # Primary is at (0,0)
    Ω = bparameters.Ω*pi/180.0; # longitude of ascending node
    i = bparameters.i*pi/180.0;
    ω = bparameters.ω*pi/180.0; # argument of periapsis
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

include("orbits.jl")
import Base: Math as math

# The following function is for single visible roche lobes (= symbiotics or large stars with hidden companions) ONLY
#
function update_roche_radii(tessels::tessellation, roche_parameters, D; use_fillout_factor = false, secondary = false, T=Float32)
    # if wanting to call this for secondary=true, invert roche_parameters.q;
    secondary == false ? potential_function = compute_potential_primary : potential_function = compute_potential_secondary;
    secondary == false ? fillout_factor = roche_parameters.fillout_factor_primary : fillout_factor = roche_parameters.fillout_factor_secondary;
    async_ratio = roche_parameters.rotation_period/roche_parameters.P
    a = roche_parameters.a;
    q = roche_parameters.q;
    rpole = roche_parameters.rpole
    # Compute surface potentials and good init
    pot_surface, r_init = get_surface_potential(rpole/a, D, q, async_ratio, fillout_factor, use_fillout_factor = use_fillout_factor, secondary=secondary);
    # Compute R_L1 for the near-L1 interpolation patch (needed when fillout > 0.997)
    ff = use_fillout_factor ? fillout_factor : rpole_to_fillout(rpole/a, D, q, async_ratio, secondary=secondary)
    R_L1 = T(0)
    if ff > 0.99
        secondary == false ? rtry = radius_point_Pathania(q) : rtry = radius_point_Pathania(1/q)
        R_L1 = solve_R_L1(rtry, D, q, async_ratio, potential_function, secondary=secondary)
    end
    # Update the radii r(θ,ϕ) to match the surface potential
    npix = tessels.npix
    r = zeros(T, npix, 5);
    for ii = 1:npix
        for jj = 1:5
            r[ii,jj] = a*solve_radius(r_init, pot_surface, D, tessels.unit_spherical[ii,jj,2], tessels.unit_spherical[ii,jj,3], q, async_ratio, potential_function,
                                       verbose=false, fillout=ff, R_L1=R_L1, secondary=secondary);
        end
    end
    return r
end


# TODO: update with the new scheme
function update_roche_radii_binary(star1_geom::tessellation, star2_geom::tessellation, binary_parameters, D, use_fillout_factor = false) # updates both primary and secondary
## TEST star1_geom = star1; star2_geom =  star2; use_fillout_factor = false;
    fillout_factor1 = binary_parameters.fillout_factor[1];
    fillout_factor2 = binary_parameters.fillout_factor[1];
    async_ratio1 = binary_parameters.star1.rotation_period/binary_parameters.P
    async_ratio2 = binary_parameters.star2.rotation_period/binary_parameters.P
    a = binary_parameters.a
    q = binary_parameters.q;

    # Compute surface potentials and good init
    potS1, r_init1 = get_surface_potential(binary_parameters.star1.rpole/a, D, q, async_ratio1, fillout_factor1, use_fillout_factor = use_fillout_factor);
    potS2, r_init2 = get_surface_potential(binary_parameters.star2.rpole/a, D, 1/q, async_ratio2, fillout_factor2, use_fillout_factor = use_fillout_factor, secondary = true);

    # Update the radii r at each (θ,ϕ) to match the surface potential
    star1_roche_geom = update_roche_geom(star1_geom, r_init1, potS1, a, D, q, async_ratio1);
    star2_roche_geom = update_roche_geom(star2_geom, r_init2, potS2, a, D, 1/q, async_ratio2, secondary = true);
    return star1_roche_geom, star2_roche_geom
end

function get_surface_potential(rpole_a, D, q, async_ratio, fillout_factor; secondary = false, use_fillout_factor = false, T=Float32)
  secondary == false ? potential_function = compute_potential_primary : potential_function = compute_potential_secondary;
  if (use_fillout_factor == true)
    #
    # If Fillout factor defines the Roche Lobe
    #
    secondary == false ? rtry = radius_point_Pathania(q) : rtry=radius_point_Pathania(1/q)
    R_L1 = solve_R_L1(rtry, D, q, async_ratio, potential_function, secondary = secondary)
    pot_L1, ~ = potential_function(R_L1, D, Int(-2*(secondary == true)+1)*pi/2, 0.0, q, async_ratio) # Primary -> π/2, Secondary -> -π/2
    potS = (pot_L1 + q * q / 2(1 + q)) / fillout_factor - q * q / 2(1 + q) # eq 6, Leahy paper (Mochnacki 1984 definition)
    return potS, rtry
  else
    #  The radius at the North pole defines the potential
    potS, ~ = potential_function(rpole_a, D, T(0), T(0), q, async_ratio)
    return potS, rpole_a
end
end

function filllout_to_rpole(fillout, D, q, async_ratio; secondary = false)
    # Note we expect calls to this with q = M2/M1 for primary, and = M1/M2 for secondary
    secondary == false ? potential_function = compute_potential_primary : potential_function = compute_potential_secondary;
    # Finds which (dimensionless) rpole corresponds to the fillout
    # Multiply by a to find the size in mas
    secondary == false ? rtry = radius_point_Pathania(q) : rtry=radius_point_Pathania(1/q)
    R_L1 = solve_R_L1(rtry, D, q, async_ratio, potential_function);
    pot_L1, ~ = potential_function(R_L1, D, pi/2.0, 0.0, q, async_ratio);
    potS = (pot_L1 + 0.5 * q * q / (1.0 + q)) / fillout - 0.5 * q * q / (1.0 + q)
    rpole = solve_radius(radius_leahy(q), potS, D, 0.0, 0.0, q, async_ratio, potential_function, verbose=false)
    return rpole # Beware, output is the reduced rpole = rpole/a
end

function rl1(roche_parameters; secondary = false)
    secondary == false ? potential_function = compute_potential_primary : potential_function = compute_potential_secondary;
    async_ratio = roche_parameters.rotation_period/roche_parameters.P
    a = roche_parameters.a;
    q = roche_parameters.q;
    D = 1.0 # circular orbit approximation
    secondary == false ? rtry = radius_point_Pathania(q) : rtry = radius_point_Pathania(1/q)
    R_L1 = solve_R_L1(rtry, D, q, async_ratio, potential_function, secondary = secondary)
    return R_L1*a
end

function max_rpole(D, roche_parameters; secondary = false)
    # The maximum rpole will be gotten for a fillout factor of 1.0
    a = roche_parameters.a;
    q = roche_parameters.q;
    async_ratio = roche_parameters.rotation_period/roche_parameters.P
    return filllout_to_rpole(1.0, D, q, async_ratio)*a;
end


function rpole_to_fillout(rpole, D, q, async_ratio; secondary = false)
    # Beware, input is the reduced rpole = rpole/a
    secondary == false ? potential_function = compute_potential_primary : potential_function = compute_potential_secondary;
    # Finds which fillout corresponds to the dimensionless rpole (=rpole/a)
    potS, ~ = potential_function(rpole, D, 0.0, 0.0, q, async_ratio);
    secondary == false ? rtry = radius_point_Pathania(q) : rtry=radius_point_Pathania(1/q)
    R_L1 = solve_R_L1(rtry, D, q, async_ratio, potential_function);
    pot_L1, ~ = potential_function(R_L1, D, Int(-2*(secondary == true)+1)*pi/2.0, 0.0, q, async_ratio);
    fillout =   (pot_L1 + q * q / 2(1.0 + q)) /  (potS + q * q / 2(1.0 + q))
    return fillout
end

# =====================================================================
# Roche Potential Functions
# =====================================================================
# Primary at (0,0,0), secondary at distance D along the x-axis.
# All quantities are dimensionless (divided by semi-major axis a).
# θ = colatitude (0 at pole), φ = azimuth (0 toward companion for primary).
# λ = sin(θ)cos(φ) = direction cosine toward companion.
# async_ratio = ω_rot/ω_orb (1 = synchronous).
#
# Reference: Aufdenberg et al. 2021, ApJ 920 130
# Potential verified against Leahy 2014 (ROCHE.F90) — see comments below.
#
# The secondary potential uses the same q = M2/M1 convention but with the
# potential centered on the secondary. When calling update_roche_radii with
# secondary=true, the caller must pass q = M1/M2 (inverted).
# The D² constant terms in Ω2 are required for correct absolute potential
# values (fillout factor computation); they cancel in root-finding for r(θ,φ).
# =====================================================================

function compute_potential_primary(r, D, θ, ϕ, q, async_ratio) # r and D are dimensionless (were divided by a)
    # Note: async_ratio = ω1rot/ωorb
    λ = sin(θ)*cos(ϕ);
    ν = cos(θ)
    Ω1 = 1/r + q/sqrt( D^2 + r^2 - 2*r*λ*D) - q*r*λ*D + async_ratio^2*(1+q)*r^2*(1-ν^2)/2
    dΩ1 = -1/r^2 - q*(r-λ*D)/sqrt(( D^2 + r^2 - 2*r*λ*D)^3) - q*λ*D + async_ratio^2*(1+q)*r*(1-ν^2)
    ddΩ1 =  2/r^3  + 3*(D*λ - r)^2*q/sqrt((D^2 + r^2 - 2*r*λ*D)^5) - q/sqrt((D^2 + r^2 - 2*r*λ*D)^3) + async_ratio^2*(q + 1)*(1 - ν^2)
    return Ω1, dΩ1, ddΩ1
end

function compute_potential_secondary(r, D, θ, ϕ, q, async_ratio)
    # Note: async_ratio = ω2rot/ωorb
    # Centered on secondary; companion (primary) is at distance D along +x (θ=π/2, φ=0).
    # For asynchronous rotation, the centrifugal term splits into orbital (D-dependent)
    # and spin (async²) contributions. At async=1 they combine to the standard synchronous form.
    λ = cos(ϕ)*sin(θ)
    ν = cos(θ)
    Ω2 = 1/sqrt(D^2+r^2+2*r*λ*D) + q/r + (1+q)/2*(D^2+2*D*r*λ) - q*(D^2+r*λ*D) + async_ratio^2*(1+q)*r^2*(1-ν^2)/2
    dΩ2 = - (D*λ + r)/sqrt((D^2+r^2+2*r*λ*D)^3) - q/r^2 + λ*D + async_ratio^2*(1+q)*r*(1-ν^2)
    ddΩ2 =  2*q/r^3 + 3*(D*λ + r)^2/sqrt((2*D*r*λ + D^2 + r^2)^5) - 1/sqrt((2*D*r*λ + D^2 + r^2)^3) + async_ratio^2*(1+q)*(1-ν^2)
    return Ω2, dΩ2, ddΩ2
end

# =====================================================================
# Root-Finding: solve_radius, Halley, Newton, Brent
# =====================================================================

"""
    solve_radius(r0, pot_surface, D, θ, ϕ, q, async_ratio, potential_function;
                 verbose=true, fillout=0.0, R_L1=0.0, secondary=false)

Find the radius r at direction (θ,ϕ) where the Roche potential equals `pot_surface`.
Uses Halley's method, with a near-L1 interpolation patch for nearly Roche-lobe-filling
stars (fillout > 0.997) to avoid convergence issues where the potential is tangent to zero.
"""
function solve_radius(r0, pot_surface, D, θ, ϕ, q, async_ratio, potential_function;
                      verbose=true, fillout=0.0, R_L1=0.0, secondary=false)
    # Near-L1 interpolation patch (following Leahy 2014, ROCHE.F90 GetPotRoot)
    # The potential is tangent to zero near L1, making root-finders unreliable.
    # For fillout > 0.997 and directions within alim of the L1 axis, linearly
    # interpolate between R_L1 and R(alim).
    if fillout > 0.997 && R_L1 > 0.0
        alim = 0.1 * pi / 180.0  # 0.1 degrees
        # L1 direction: θ=π/2, φ=0 for primary; θ=-π/2, φ=0 for secondary
        θ_L1 = secondary ? -pi/2 : pi/2
        # Angular distance from L1 axis (small angle approximation)
        dθ = θ - θ_L1
        # Handle φ wrapping for secondary (L1 is at φ=0 in both frames)
        dφ = mod(ϕ + pi, 2pi) - pi  # wrap to [-π, π]
        ang = sqrt(dθ^2 + dφ^2)
        if ang < alim
            # Solve at a safe angle offset from L1
            fgh_safe = r -> potential_function(r, D, θ_L1 + alim, 0.0, q, async_ratio)
            R_alim = halley(r0, pot_surface, fgh_safe, verbose=false)
            # Linear interpolation
            return R_L1 * (1.0 - ang/alim) + R_alim * (ang/alim)
        end
    end
    fgh = r->potential_function(r, D, θ, ϕ, q, async_ratio);
    r = halley(r0, pot_surface, fgh, verbose=verbose);
    return r
end


function solve_R_L1(r0, D, q, async_ratio, potential_function; secondary = false, verbose = false, ϵ=1e-6)
    sg = Int(-2*(secondary == true)+1)
    fg = r->potential_function(r, D, sg*pi/2, 0.0, q, async_ratio);
    r = newton_root(r0, fg)
    if r<0
        @warn "Negative R_L1 - unphysical Roche parameters"
    end
    if abs(fg(r)[2])>ϵ
        @warn "R_L1 tol exceeded $(abs(fg(r)[1]))>$ϵ"
    end
    return r
end

"""
    solve_R_L2(D, q, async_ratio) -> Float64

Find the L2 Lagrange point (behind the primary, away from secondary).
Returns the radial distance from the primary along θ=π/2, φ=π.
"""
function solve_R_L2(D, q, async_ratio)
    # L2 is behind the primary (direction away from secondary: θ=π/2, φ=π)
    # dΩ/dr = 0 along this direction
    fg = r -> compute_potential_primary(r, D, pi/2, pi, q, async_ratio)
    # Bracket: R_L2 is between 0 and ~D for typical mass ratios
    r0 = 0.5 * radius_eggleton(q) * D  # initial guess smaller than R_L1
    r = newton_root(r0, fg)
    return r
end

"""
    solve_R_L3(D, q, async_ratio) -> Float64

Find the L3 Lagrange point (behind the secondary, past x=D).
Returns the radial distance from the primary along θ=π/2, φ=0 (r > D).
"""
function solve_R_L3(D, q, async_ratio)
    # L3 is behind the secondary: along θ=π/2, φ=0 but at r > D
    # We use Brent's method because Newton can overshoot past the singularity at r=D
    fg = r -> compute_potential_primary(r, D, pi/2, 0.0, q, async_ratio)
    deriv = r -> fg(r)[2]
    # Bracket: L3 is between D and some large radius
    r = brent_root(deriv, D * 1.001, D * 100.0)
    return r
end

"""
    solve_lagrange_points(D, q, async_ratio) -> (R_L1, R_L2, R_L3)

Find all three collinear Lagrange points. Returns radial distances from the primary:
- R_L1: toward the secondary (0 < R_L1 < D)
- R_L2: away from secondary, behind primary (along φ=π)
- R_L3: behind secondary (R_L3 > D, along φ=0)
"""
function solve_lagrange_points(D, q, async_ratio)
    rtry = radius_point_Pathania(q)
    R_L1 = solve_R_L1(rtry, D, q, async_ratio, compute_potential_primary)
    R_L2 = solve_R_L2(D, q, async_ratio)
    R_L3 = solve_R_L3(D, q, async_ratio)
    return R_L1, R_L2, R_L3
end


function newton_root(x0, fgh; ϵ = 1e-5, verbose=false)
    # Get x so that f'(x) = 0
    n = 1;
    converged = false;
    x = copy(x0)
    while ((converged == false) & (n < 20))
        f, g, h = fgh(x);
        newton_step = g/h;
        if verbose == true
            println("n = $n\t f = $f \t g = $g \t  h = $h \t x = $x \t step = $newton_step");
        end
        x -= newton_step;
        if (abs(newton_step) < ϵ)
            converged = true;
        end
        n += 1;
    end
   return x
end


function halley(x0, f0, fgh; ϵ = 1e-5, verbose=false)
    # Get x so that f(x) = f0
    n = 1;
    converged = false;
    x = copy(x0)
    while ((converged == false) & (n < 20))
        f, g, h = fgh(x);
        halley_step = 2*(f - f0)*g/( 2*g^2 - (f-f0)*h);
        if verbose == true
            println("n = $n\t f = $f,\t f0 = $f0 \t g = $g \t  h = $h \t x = $x \t step = $halley_step");
        end
        x -= halley_step;
        if (abs(halley_step) < ϵ)
            converged = true;
        end
        n += 1;
    end
   return x
end

"""
    brent_root(f, a, b; tol=1e-9, maxiter=100) -> Float64

Find a root of `f(x) = 0` in the interval [a, b] using Brent's method.
Requires `f(a)` and `f(b)` to have opposite signs.
"""
function brent_root(f, a, b; tol=1e-9, maxiter=100)
    fa = f(a); fb = f(b)
    if fa * fb >= 0
        error("brent_root: f(a) and f(b) must have opposite signs (f($a)=$fa, f($b)=$fb)")
    end
    # Ensure |f(a)| >= |f(b)|
    if abs(fa) < abs(fb)
        a, b = b, a
        fa, fb = fb, fa
    end
    c = a; fc = fa
    d = b - a
    mflag = true
    for _ in 1:maxiter
        if abs(b - a) <= tol || fb == 0.0
            return b
        end
        # Inverse quadratic interpolation or secant
        if fa != fc && fb != fc
            s = a*fb*fc/((fa-fb)*(fa-fc)) + b*fa*fc/((fb-fa)*(fb-fc)) + c*fa*fb/((fc-fa)*(fc-fb))
        else
            s = b - fb*(b-a)/(fb-fa)
        end
        # Conditions for bisection fallback
        cond1 = !(((3a+b)/4 < s < b) || (b < s < (3a+b)/4))
        cond2 = mflag && abs(s-b) >= abs(b-c)/2
        cond3 = !mflag && abs(s-b) >= abs(c-d)/2
        cond4 = mflag && abs(b-c) < tol
        cond5 = !mflag && abs(c-d) < tol
        if cond1 || cond2 || cond3 || cond4 || cond5
            s = (a+b)/2
            mflag = true
        else
            mflag = false
        end
        fs = f(s)
        d = c; c = b; fc = fb
        if fa*fs < 0
            b = s; fb = fs
        else
            a = s; fa = fs
        end
        if abs(fa) < abs(fb)
            a, b = b, a
            fa, fb = fb, fa
        end
    end
    @warn "brent_root: maxiter reached"
    return b
end

# =====================================================================
# Gravity Functions (for von Zeipel temperature maps)
# =====================================================================

function compute_gravity_primary(r,θ,ϕ,D,q,async_ratio)
    # r = dimensionless radius
    λ = cos.(ϕ).*sin.(θ)
    μ = sin.(ϕ).*sin.(θ)
    ν = cos.(θ)
    x = r.*λ
    y = r.*μ
    z = r.*ν
    # This means r^2 = x^2 + y^2 + z^2
    ρ = sqrt.(((D.-x).^2+y.^2+z.^2).^(-3)) # faster than .^(-1.5)
    μ = r.^(-3)
    gx = -x.*μ + q*(D.-x).*ρ + async_ratio^2*(1+q)*x.-q*D
    gy = -y.*μ - q*y.*ρ +async_ratio^2*(1+q)*y
    gz = -z.*μ - q*z.*ρ
    return sqrt.(gx.*gx + gy.*gy + gz.*gz);
end

function compute_gravity_secondary(r,θ,ϕ,D,q,async_ratio)
    # r = dimensionless radius
    λ = cos.(ϕ).*sin.(θ)
    μ = sin.(ϕ).*sin.(θ)
    ν = cos.(θ)
    x = D .+ r.*λ
    y = r.*μ
    z = r.*ν
    # This means r^2 = (D-x)^2 + y^2 + z^2
    ρ = sqrt.((x.^2+y.^2+z.^2).^(-3))
    μ = r.^(-3)
    gx = -x.*ρ + q*(D.-x).*μ + (1-async_ratio^2)*(1+q)*(D.-x)+ x*(1+q).-q*D
    gy = -y.*ρ - q*y.*μ + async_ratio^2*(1+q)*y
    gz = -z.*ρ - q*z.*μ
    return sqrt.(gx.*gx + gy.*gy + gz.*gz);
end

# =====================================================================
# Temperature Maps
# =====================================================================

function temperature_map_vonZeipel_roche_single(parameters, star_geom, t::Array{T,1};  secondary = false) where T
    npix = star_geom[1].npix
    nepochs = length(t)
    Tmap = zeros(T, npix, nepochs)
    for i=1:nepochs
        Tmap[:, i] = temperature_map_vonZeipel_roche_single(parameters, star_geom[i], t[i], secondary = secondary, T=T)
    end
    return Tmap
end

function temperature_map_vonZeipel_roche_single(parameters, star_geom, t; secondary = false, T=Float32)
    p = convert_params(T, parameters)
    rpole = p.rpole/p.a
    r = star_geom.vertices_spherical[:, 5 ,1]/p.a
    θ = star_geom.vertices_spherical[:, 5, 2]
    ϕ = star_geom.vertices_spherical[:, 5, 3]
    D = T(compute_separation(p, t))
    async_ratio = p.rotation_period/p.P
    # Compute gravity
    compute_gravity = compute_gravity_primary;
    if secondary == true
        compute_gravity = compute_gravity_secondary;
    end
    g_pole = compute_gravity(rpole,T(0),T(0),D,p.q,async_ratio);
    gravity_map = compute_gravity(r,θ,ϕ,D,p.q,async_ratio);
    # Computes von Zeipel temperature map directly from the gravity map
    Teff = p.tpole*(gravity_map./g_pole).^p.beta;
    return Teff
end

function temperature_map_vonZeipel_roche(binary_parameters, star_geom, t; secondary = false)
    T = eltype(star_geom.vertices_spherical)
    sp = secondary ? binary_parameters.star2 : binary_parameters.star1
    rpole = T(sp.rpole/binary_parameters.a)
    tpole = T(sp.tpole)
    r = star_geom.vertices_spherical[:, 5 ,1]/T(binary_parameters.a)
    θ = star_geom.vertices_spherical[:, 5, 2]
    ϕ = star_geom.vertices_spherical[:, 5, 3]
    D = T(compute_separation(binary_parameters, t))
    async_ratio = T(sp.rotation_period/binary_parameters.P)
    q = T(binary_parameters.q)
    β = T(sp.beta)
    # Compute gravity
    compute_gravity = secondary ? compute_gravity_secondary : compute_gravity_primary;
    g_pole = compute_gravity(rpole,T(0),T(0),D,q,async_ratio);
    gravity_map = compute_gravity(r,θ,ϕ,D,q,async_ratio);
     # Computes von Zeipel temperature map directly from the gravity map
    Teff = tpole*(gravity_map./g_pole).^β;
    return Teff
end

# =====================================================================
# Analytical Roche Radius Approximations
# =====================================================================

# Eggleton (1983) formula for effective Roche lobe radius.
# Input: q = M2/M1 (companion/donor). Returns r_L/a for the primary (M1).
# Matches Leahy 2014 (ROCHE.F90) REg function.
function radius_eggleton(q)
    q1 = 1/q
    return 0.49*q1^(2/3)/(0.6*q1^(2/3)+log(1.0+q1^(1/3)));
end
const radius_equivalent_eggleton = radius_eggleton  # backward compat alias

function radius_point_kopal_polynomial(x)
    # Polynomial = 0 at x=point radius
    return (q + 1)x^5 - (3q + 2)x^4 + (3q + 1)x^3 - x^2 + 2x - 1
end

# Leahy & Leahy (2015) approximation for the L1 radius
function radius_leahy(q)
    a1 = 0.64334;
    a2 = 0.86907;
    a3 = 1.2809;
    a4 = −0.74303;
    a5 = 0.73103;
    return a1*q^a4/(a2*q^a5+log(1+a3*q^(a4+1/3)))
end

# Pathania & Medupe (2012) approximations
function radius_back_Pathania(q)
    return exp(-1.00598q^3 + 2.09674q^2 - 1.69263q - 0.319909)
end

function radius_point_Pathania(q)
    return exp(-0.725742q^3 + 1.53893q^2 - 1.31638q - 0.202505)
end


# =====================================================================
# Roche Lobe Volume and Area (Romberg integration, following Leahy 2014)
# =====================================================================

"""
    romberg_integrate(f, a, b; N=7) -> Float64

Single-variable Romberg integration of `f` over [a, b] with `N` refinement levels.
Uses 2^(N-1) function evaluations. Equivalent to the Fortran ROMBERG subroutines.
"""
function romberg_integrate(f, a, b; N=7)
    R = zeros(N, N)
    h = b - a
    R[1,1] = 0.5 * h * (f(a) + f(b))
    L = 1
    for i in 2:N
        h *= 0.5
        L *= 2
        s = 0.0
        for k in 1:2:(L-1)
            s += f(a + h * k)
        end
        R[i,1] = 0.5 * R[i-1,1] + h * s
        m = 1
        for j in 2:i
            m *= 4
            R[i,j] = R[i,j-1] + (R[i,j-1] - R[i-1,j-1]) / (m - 1)
        end
    end
    return R[N,N]
end

"""
    roche_volume(q, async_ratio, D; fillout=1.0, secondary=false) -> Float64

Compute the dimensionless Roche lobe volume (in units of a³).
Uses double Romberg integration over (cos θ, φ), following Leahy (2014).
"""
function roche_volume(q, async_ratio, D; fillout=1.0, secondary=false)
    potential_function = secondary ? compute_potential_secondary : compute_potential_primary
    secondary ? rtry = radius_point_Pathania(1/q) : rtry = radius_point_Pathania(q)
    R_L1 = solve_R_L1(rtry, D, q, async_ratio, potential_function, secondary=secondary)
    sg = Int(-2*(secondary)+1)
    pot_L1, _ = potential_function(R_L1, D, sg*pi/2, 0.0, q, async_ratio)
    OmegaF = (pot_L1 + q^2 / (2*(1+q))) / fillout - q^2 / (2*(1+q))
    r_init = rtry

    function vol_phi(θ, φ)
        r = solve_radius(r_init, OmegaF, D, θ, φ, q, async_ratio, potential_function,
                         verbose=false, fillout=fillout, R_L1=R_L1, secondary=secondary)
        return (2.0/3.0) * r^3
    end

    function vol_mu(mu)
        θ = acos(mu)
        # Integrate over φ ∈ [0, π], multiply by 2 for φ-symmetry
        return 2.0 * romberg_integrate(φ -> vol_phi(θ, φ), 0.0, pi)
    end

    # Integrate over μ = cos(θ) ∈ [0, 1], multiply by 2 for hemisphere symmetry
    return 2.0 * romberg_integrate(vol_mu, 0.0, 1.0)
end

"""
    roche_area(q, async_ratio, D; fillout=1.0, secondary=false) -> Float64

Compute the dimensionless Roche lobe surface area (in units of a²).
Uses double Romberg integration over (cos θ, φ), following Leahy (2014).
"""
function roche_area(q, async_ratio, D; fillout=1.0, secondary=false)
    potential_function = secondary ? compute_potential_secondary : compute_potential_primary
    secondary ? rtry = radius_point_Pathania(1/q) : rtry = radius_point_Pathania(q)
    R_L1 = solve_R_L1(rtry, D, q, async_ratio, potential_function, secondary=secondary)
    sg = Int(-2*(secondary)+1)
    pot_L1, _ = potential_function(R_L1, D, sg*pi/2, 0.0, q, async_ratio)
    OmegaF = (pot_L1 + q^2 / (2*(1+q))) / fillout - q^2 / (2*(1+q))
    r_init = rtry

    function area_phi(θ, φ)
        r = solve_radius(r_init, OmegaF, D, θ, φ, q, async_ratio, potential_function,
                         verbose=false, fillout=fillout, R_L1=R_L1, secondary=secondary)
        return 2.0 * r^2
    end

    function area_mu(mu)
        θ = acos(mu)
        return 2.0 * romberg_integrate(φ -> area_phi(θ, φ), 0.0, pi)
    end

    return 2.0 * romberg_integrate(area_mu, 0.0, 1.0)
end

"""
    roche_equivalent_radius(q, async_ratio, D; fillout=1.0, secondary=false) -> (R_vol, R_area)

Compute volume-equivalent and area-equivalent radii (dimensionless, in units of a).
"""
function roche_equivalent_radius(q, async_ratio, D; fillout=1.0, secondary=false)
    vol = roche_volume(q, async_ratio, D; fillout=fillout, secondary=secondary)
    area = roche_area(q, async_ratio, D; fillout=fillout, secondary=secondary)
    R_vol = (3*vol / (4*pi))^(1/3)
    R_area = (area / (4*pi))^(1/2)
    return R_vol, R_area
end

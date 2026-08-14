include("orbits.jl")
import Base: Math as math

"""
    synchronicity(p) -> F

Synchronicity parameter `F = ω_rot/ω_orb = P_orb/P_rot` — 1 synchronous, >1
supersynchronous, <1 subsynchronous. This is the quantity that enters the Roche
potential's centrifugal term as `½F²(1+q)r²(1−ν²)`, and it is the same `F` as PHOEBE's
`syncpar`, ELISa's `synchronicity`, and `P` in `RocheLobe.f90`.

Note it is the ratio of *angular rates*, i.e. the RECIPROCAL of the ratio of periods:
`P/rotation_period`, NOT `rotation_period/P`. The inverted form is invisible for a
synchronous binary (both give 1), which is every configuration in the demos and the tests,
so it will not show up in a test run — but for Spica's primary (v sin i = 161 km/s ⇒
F ≈ 1.92) it puts the centrifugal term out by F⁴ ≈ 13.6× and the rotational flattening at
0.5% instead of 9.4%.

Verified against `libphoebe.roche_Omega(q, F, d, [x,y,z])` at F = 0.521, 1.0 and 1.92 —
see the "synchronicity convention" testset in `test/test_reflection.jl`.
"""
@inline synchronicity(p) = p.P / p.rotation_period

# The following function is for single visible roche lobes (= symbiotics or large stars with hidden companions) ONLY
#
function update_roche_radii(tessels::tessellation, roche_parameters, D; use_fillout_factor = false, secondary = false, T=eltype(tessels), omega = nothing)
    # if wanting to call this for secondary=true, invert roche_parameters.q;
    secondary == false ? potential_function = compute_potential_primary : potential_function = compute_potential_secondary;
    secondary == false ? fillout_factor = roche_parameters.fillout_factor_primary : fillout_factor = roche_parameters.fillout_factor_secondary;
    async_ratio = synchronicity(roche_parameters)
    a = roche_parameters.a;
    q = roche_parameters.q;
    rpole = roche_parameters.rpole
    # Compute surface potentials and good init.
    # `omega` (if given) pins the equipotential directly — the volume-conserving path
    # (roche_omega_for_volume) uses it so the star's volume, not its polar radius, is the
    # invariant as D(t) changes around an eccentric orbit.
    if omega === nothing
        pot_surface, r_init = get_surface_potential(rpole/a, D, q, async_ratio, fillout_factor, use_fillout_factor = use_fillout_factor, secondary=secondary);
    else
        pot_surface, r_init = omega, rpole/a
    end
    # Update the radii r(θ,ϕ) to match the surface potential
    npix = tessels.npix
    r = zeros(T, npix, 5);
    for ii = 1:npix
        for jj = 1:5
            r[ii,jj] = a*solve_radius(r_init, pot_surface, D, tessels.unit_spherical[ii,jj,2], tessels.unit_spherical[ii,jj,3], q, async_ratio, potential_function,
                                       verbose=false);
        end
    end
    return r
end


# TODO: update with the new scheme
function update_roche_radii_binary(star1_geom::tessellation, star2_geom::tessellation, binary_parameters, D, use_fillout_factor = false) # updates both primary and secondary
## TEST star1_geom = star1; star2_geom =  star2; use_fillout_factor = false;
    fillout_factor1 = binary_parameters.fillout_factor[1];
    fillout_factor2 = binary_parameters.fillout_factor[1];
    async_ratio1 = binary_parameters.P/binary_parameters.star1.rotation_period
    async_ratio2 = binary_parameters.P/binary_parameters.star2.rotation_period
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

function get_surface_potential(rpole_a, D, q, async_ratio, fillout_factor; secondary = false, use_fillout_factor = false,
                               T = float(promote_type(typeof(rpole_a), typeof(D), typeof(q))))
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

function fillout_to_rpole(fillout, D, q, async_ratio; secondary = false)
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
    async_ratio = synchronicity(roche_parameters)
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
    async_ratio = synchronicity(roche_parameters)
    return fillout_to_rpole(1.0, D, q, async_ratio)*a;
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
# Reference: Aufdenberg et al. 2015 (papers/spica_the_paper_2015.pdf), appendix A;
# cross-checked against /home/baron/SOFTWARE/roche/RocheLobe.f90 (`potential`) and
# PHOEBE 2 (`phoebe/lib/gen_roche.h`).
#
# ---------------------------------------------------------------------------
# The tidal (free-fall) term −q·r·λ/D² — where it comes from
# ---------------------------------------------------------------------------
# The frame is centred on the star itself, NOT on the centre of mass. That origin is in
# free fall, accelerating toward the companion at G·M₂/D². Working in that non-inertial
# frame adds a uniform pseudo-force, whose potential is −G·M₂·x/D² — in units of G·M₁/a,
#
#     −q·x/D²  =  −q·r·λ/D²                                     (the linear term below)
#
# This is a purely translational, tidal term: it involves no assumption whatsoever about
# how fast the frame rotates. The rotational term is then separately the centrifugal
# potential about the star's OWN spin axis, ½F²(1+q)r²(1−ν²) with F = ω_rot/ω_orb
# referenced to the MEAN orbital rate — hence carrying no D dependence.
#
# Getting this framing right matters, because the tempting alternative derivation — put
# the origin at the centre of mass and write the linear term as ω²·x_cm·x with
# x_cm = D·q/(1+q) — is NOT self-consistent with the code. It would force the rotational
# coefficient to be ω² = (1+q)/D³ as well, i.e. D-dependent, which is not what is
# implemented and not what any reference code does. (And (1+q)/D³ is in any case the
# *circular* rate at separation D, not the instantaneous orbital rate, which is larger by
# a factor 1 + e·cos ν.)  Reference: Sepinsky, Willems & Kalogera 2007, ApJ 660, 1624,
# who derive exactly this frame for nonsynchronous eccentric binaries.
#
# Three independent modern implementations agree on both terms:
#   PHOEBE 2  `phoebe/lib/gen_roche.h`:  q[(δ²+ρ²−2ρλδ)^(−1/2) − ρλ/δ²] + ½F²(1+q)ρ²(1−ν²)
#   ELISa     `elisa/binary_system/model.py`, pre_calculate_for_potential_value_primary:
#             b = d²;  d_coef = q·cs/b = q·λ/D²;  e = ½F²(1+q)sin²θ
#             Psi1 = 1/r + q/√(b+r²−c·r) − d_coef·r + e·r²
#   Wilson 1979 (the canonical source both descend from)
#
# NOTE this disagrees with Aufdenberg et al. 2015 eqs. A18/A27/A30, which write the term
# as −q·r·λ·D (internally self-consistent with his A1, but wrong by D³ relative to the
# above). That form overwhelms the tidal attraction past a q- and radius-dependent
# threshold and points the tidal bulge AWAY from the companion: for Spica
# (q = 0.619, rpole/a = 0.290) the threshold is D/a = 1.031, and every modern eccentricity
# determination (0.065–0.123) puts apastron beyond it, so the bulge inverts for roughly
# half of every orbit.
#
# RocheLobe.f90 does not settle the question directly: it works in units of the
# instantaneous separation, so D ≡ 1 and the symbol never appears. But rescaling it to
# units of the semi-major axis (r = D·s) reproduces the /D² form exactly, to machine
# precision — an independent confirmation. Note that in D-scaled lengths the rotational
# coefficient picks up δ³ (PHOEBE writes `b = F*F*delta^3*(1+q)`), so a future switch to
# that convention must carry F² → F²δ³.
#
# The secondary potential is centred on the secondary and normalised by G·M1 (as in
# Aufdenberg A16/A19), so the companion term has coefficient 1 and the self term
# coefficient q. The D² constant in Ω2 is required for correct absolute potential values
# (fillout factor); it cancels in root-finding for r(θ,φ).
# =====================================================================

function compute_potential_primary(r, D, θ, ϕ, q, async_ratio) # r and D are dimensionless (were divided by a)
    # Note: async_ratio = ω1rot/ωorb, referenced to the MEAN orbital rate
    λ = sin(θ)*cos(ϕ);
    ν = cos(θ)
    invD2 = 1/D^2                      # instantaneous-rate centre-of-mass term (see above)
    Ω1 = 1/r + q/sqrt( D^2 + r^2 - 2*r*λ*D) - q*r*λ*invD2 + async_ratio^2*(1+q)*r^2*(1-ν^2)/2
    dΩ1 = -1/r^2 - q*(r-λ*D)/sqrt(( D^2 + r^2 - 2*r*λ*D)^3) - q*λ*invD2 + async_ratio^2*(1+q)*r*(1-ν^2)
    ddΩ1 =  2/r^3  + 3*(D*λ - r)^2*q/sqrt((D^2 + r^2 - 2*r*λ*D)^5) - q/sqrt((D^2 + r^2 - 2*r*λ*D)^3) + async_ratio^2*(q + 1)*(1 - ν^2)
    return Ω1, dΩ1, ddΩ1
end

function compute_potential_secondary(r, D, θ, ϕ, q, async_ratio)
    # Note: async_ratio = ω2rot/ωorb, referenced to the MEAN orbital rate.
    # Centred on the secondary; the companion (primary) lies at distance D along −x,
    # hence the +2rλD in the companion term (and L1 sits at λ = −1, see solve_R_L1).
    λ = cos(ϕ)*sin(θ)
    ν = cos(θ)
    invD2 = 1/D^2
    const_D = (1-q)*D^2/2              # = (1+q)D²/2 − qD²; a pure constant at fixed D
    Ω2 = 1/sqrt(D^2+r^2+2*r*λ*D) + q/r + r*λ*invD2 + const_D + async_ratio^2*(1+q)*r^2*(1-ν^2)/2
    dΩ2 = - (D*λ + r)/sqrt((D^2+r^2+2*r*λ*D)^3) - q/r^2 + λ*invD2 + async_ratio^2*(1+q)*r*(1-ν^2)
    ddΩ2 =  2*q/r^3 + 3*(D*λ + r)^2/sqrt((2*D*r*λ + D^2 + r^2)^5) - 1/sqrt((2*D*r*λ + D^2 + r^2)^3) + async_ratio^2*(1+q)*(1-ν^2)
    return Ω2, dΩ2, ddΩ2
end

# =====================================================================
# Root-Finding: solve_radius, Halley, Newton, Brent
# =====================================================================

"""
    solve_radius(r0, pot_surface, D, θ, ϕ, q, async_ratio, potential_function;
                 verbose=true)

Find the radius r at direction (θ,ϕ) where the Roche potential equals `pot_surface`.
Uses Halley's method (cubic convergence). Converges reliably even near the L1 point
at all fillout levels up to and including 1.0.
"""
function solve_radius(r0, pot_surface, D, θ, ϕ, q, async_ratio, potential_function;
                      verbose=true, fillout=0.0, R_L1=0.0, secondary=false)
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
!!! note "Returns `Float64` on purpose"
    Not an oversight and not a candidate for the element-type-follows-input rule: this is a
    root solve / quadrature whose nested Halley stack needs the headroom. Float32 inputs are
    accepted and widened.
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
!!! note "Returns `Float64` on purpose"
    Not an oversight and not a candidate for the element-type-follows-input rule: this is a
    root solve / quadrature whose nested Halley stack needs the headroom. Float32 inputs are
    accepted and widened.
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

# |∇Ω| at each surface point. These MUST be the exact gradient of the potentials above —
# in particular the constant in gx is the derivative of the centre-of-mass term, so it
# carries the same 1/D² (not ·D) as the potential. `invr3`/`invr3c` are named explicitly
# rather than reusing `μ`, which is already the y direction cosine.

function compute_gravity_primary(r,θ,ϕ,D,q,async_ratio)
    # r = dimensionless radius; body-centred Cartesian, companion at +D x̂
    λ = cos.(ϕ).*sin.(θ)
    μ = sin.(ϕ).*sin.(θ)
    ν = cos.(θ)
    x = r.*λ
    y = r.*μ
    z = r.*ν
    # This means r^2 = x^2 + y^2 + z^2
    invr3c = sqrt.(((D.-x).^2+y.^2+z.^2).^(-3)) # 1/r_companion^3; faster than .^(-1.5)
    invr3  = r.^(-3)
    gx = -x.*invr3 + q*(D.-x).*invr3c + async_ratio^2*(1+q)*x .- q/D^2
    gy = -y.*invr3 - q*y.*invr3c + async_ratio^2*(1+q)*y
    gz = -z.*invr3 - q*z.*invr3c
    return sqrt.(gx.*gx + gy.*gy + gz.*gz);
end

function compute_gravity_secondary(r,θ,ϕ,D,q,async_ratio)
    # r = dimensionless radius; secondary-centred Cartesian (xs,y,z), companion at −D x̂.
    λ = cos.(ϕ).*sin.(θ)
    μ = sin.(ϕ).*sin.(θ)
    ν = cos.(θ)
    xs = r.*λ                 # measured from the SECONDARY's centre
    y  = r.*μ
    z  = r.*ν
    xc = xs .+ D              # measured from the companion (primary) at −D
    invr3c = sqrt.((xc.^2+y.^2+z.^2).^(-3))   # 1/r_companion^3
    invr3  = r.^(-3)
    gx = -xc.*invr3c - q*xs.*invr3 + async_ratio^2*(1+q)*xs .+ 1/D^2
    gy = -y.*invr3c  - q*y.*invr3  + async_ratio^2*(1+q)*y
    gz = -z.*invr3c  - q*z.*invr3
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

function temperature_map_vonZeipel_roche_single(parameters, star_geom, t; secondary = false, T=eltype(star_geom))
    p = convert_params(T, parameters)
    rpole = p.rpole/p.a
    r = star_geom.vertices_spherical[:, 5 ,1]/p.a
    θ = star_geom.vertices_spherical[:, 5, 2]
    ϕ = star_geom.vertices_spherical[:, 5, 3]
    D = T(compute_separation(p, t))
    async_ratio = synchronicity(p)
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
    async_ratio = T(binary_parameters.P/sp.rotation_period)
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
!!! note "Returns `Float64` on purpose"
    Not an oversight and not a candidate for the element-type-follows-input rule: this is a
    root solve / quadrature whose nested Halley stack needs the headroom. Float32 inputs are
    accepted and widened.
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
    roche_polar_radius(OmegaF, D, q, async_ratio, potential_function) -> r_pole

Radius of the `Ω = OmegaF` equipotential along the rotation axis (θ = 0), found by
bracketing outward from the centre and then Brent.

This is the *minimum* radius on the equipotential, which makes it the only universally
safe starting point for the Halley solve in other directions: starting from anything
larger (in particular from the L1 radius, which lies outside the star for an under-filled
lobe) lets Halley step over the saddle and converge to a spurious root beyond L1 —
silently returning a radius roughly twice too large near the companion.
"""
function roche_polar_radius(OmegaF, D, q, async_ratio, potential_function)
    f(r) = potential_function(r, D, 0.0, 0.0, q, async_ratio)[1] - OmegaF
    # Ω → +∞ like 1/r as r → 0 and decreases outward along the pole, so start where the
    # 1/r term dominates OmegaF by a factor 1e6 and march out. The starting radius must
    # scale with 1/Ω, NOT with D: for a widely separated pair (D ≫ 1) a D-scaled guess
    # lands far outside the star and the bracket fails.
    rlo = 1e-6 / max(abs(OmegaF), 1e-12)
    if !(f(rlo) > 0)
        @warn "roche_polar_radius: Ω = $OmegaF unreachable (D=$D, q=$q)"
        return rlo
    end
    rhi = rlo
    for _ in 1:400
        rhi *= 1.3
        f(rhi) < 0 && return brent_root(f, rlo, rhi; tol=1e-14 * rhi)
        rlo = rhi
    end
    @warn "roche_polar_radius: no bracket found for Ω = $OmegaF (D=$D, q=$q)"
    return rlo
end

"""
    roche_volume(q, async_ratio, D; fillout=1.0, secondary=false, omega=nothing) -> Float64

Dimensionless Roche lobe volume (units of a³), by double Romberg quadrature over
(cos θ, φ) following Leahy (2014).

By default the surface equipotential comes from `fillout` (relative to L1); pass `omega`
to give it directly.

Two corrections relative to the original implementation, both verified against a direct
HEALPix mesh sum `(4π/3npix)·Σrᵢ³`:

* the quadrature carried an extra factor of 2 (the μ- and φ-symmetry factors were applied
  on top of an integrand that already included them), so every volume came out exactly 2×
  too large — and hence `R_vol` 2^(1/3) ≈ 1.26× too large;
* the Halley solve was seeded with the L1 radius; see [`roche_polar_radius`](@ref).
!!! note "Returns `Float64` on purpose"
    Not an oversight and not a candidate for the element-type-follows-input rule: this is a
    root solve / quadrature whose nested Halley stack needs the headroom. Float32 inputs are
    accepted and widened.
"""
function roche_volume(q, async_ratio, D; fillout=1.0, secondary=false, omega=nothing)
    potential_function = secondary ? compute_potential_secondary : compute_potential_primary
    secondary ? rtry = radius_point_Pathania(1/q) : rtry = radius_point_Pathania(q)
    if omega === nothing
        R_L1 = solve_R_L1(rtry, D, q, async_ratio, potential_function, secondary=secondary)
        sg = Int(-2*(secondary)+1)
        pot_L1, _ = potential_function(R_L1, D, sg*pi/2, 0.0, q, async_ratio)
        OmegaF = (pot_L1 + q^2 / (2*(1+q))) / fillout - q^2 / (2*(1+q))
    else
        OmegaF = omega
    end
    r_init = roche_polar_radius(OmegaF, D, q, async_ratio, potential_function)

    # V = ∫∫∫ r² dr dΩ = (1/3) ∮ r³ dΩ, and with both symmetries
    #   = (1/3)·2∫₀¹dμ·2∫₀^π dφ  r³  = (4/3) ∫₀¹∫₀^π r³
    function vol_phi(θ, φ)
        # `solve_radius` ignores fillout/R_L1 (accepted only for call compatibility),
        # and R_L1 is not even computed on the `omega` path — so do not pass them.
        r = solve_radius(r_init, OmegaF, D, θ, φ, q, async_ratio, potential_function,
                         verbose=false, secondary=secondary)
        return r^3
    end
    vol_mu(mu) = romberg_integrate(φ -> vol_phi(acos(mu), φ), 0.0, pi)
    return (4.0/3.0) * romberg_integrate(vol_mu, 0.0, 1.0)
end

"""
    roche_area(q, async_ratio, D; fillout=1.0, secondary=false, omega=nothing) -> Float64

Dimensionless Roche lobe surface area (units of a²), by double Romberg quadrature.

Carries the same two fixes as [`roche_volume`](@ref) (spurious factor of 2; L1-seeded
root solve) and additionally the projection factor that the original omitted: for a
star-shaped surface `dA = r²/cos γ dΩ`, where γ is the angle between `r̂` and the surface
normal, so `cos γ = |∂Ω/∂r| / |∇Ω|`. Dropping it underestimates the area of a distorted
lobe.
!!! note "Returns `Float64` on purpose"
    Not an oversight and not a candidate for the element-type-follows-input rule: this is a
    root solve / quadrature whose nested Halley stack needs the headroom. Float32 inputs are
    accepted and widened.
"""
function roche_area(q, async_ratio, D; fillout=1.0, secondary=false, omega=nothing)
    potential_function = secondary ? compute_potential_secondary : compute_potential_primary
    compute_gravity = secondary ? compute_gravity_secondary : compute_gravity_primary
    secondary ? rtry = radius_point_Pathania(1/q) : rtry = radius_point_Pathania(q)
    if omega === nothing
        R_L1 = solve_R_L1(rtry, D, q, async_ratio, potential_function, secondary=secondary)
        sg = Int(-2*(secondary)+1)
        pot_L1, _ = potential_function(R_L1, D, sg*pi/2, 0.0, q, async_ratio)
        OmegaF = (pot_L1 + q^2 / (2*(1+q))) / fillout - q^2 / (2*(1+q))
    else
        OmegaF = omega
    end
    r_init = roche_polar_radius(OmegaF, D, q, async_ratio, potential_function)

    # A = ∮ r²/cos γ dΩ = 2∫₀¹dμ · 2∫₀^π dφ  r²/cos γ
    function area_phi(θ, φ)
        r = solve_radius(r_init, OmegaF, D, θ, φ, q, async_ratio, potential_function,
                         verbose=false, secondary=secondary)
        dΩdr = potential_function(r, D, θ, φ, q, async_ratio)[2]
        g = compute_gravity(r, θ, φ, D, q, async_ratio)
        cosγ = g > 0 ? min(abs(dΩdr) / g, 1.0) : 1.0
        return r^2 / max(cosγ, 1e-3)
    end
    area_mu(mu) = romberg_integrate(φ -> area_phi(acos(mu), φ), 0.0, pi)
    return 4.0 * romberg_integrate(area_mu, 0.0, 1.0)
end

"""
    roche_equivalent_radius(q, async_ratio, D; fillout=1.0, secondary=false, omega=nothing)
        -> (R_vol, R_area)

Volume-equivalent and area-equivalent radii (dimensionless, in units of a).
"""
function roche_equivalent_radius(q, async_ratio, D; fillout=1.0, secondary=false, omega=nothing)
    vol = roche_volume(q, async_ratio, D; fillout=fillout, secondary=secondary, omega=omega)
    area = roche_area(q, async_ratio, D; fillout=fillout, secondary=secondary, omega=omega)
    R_vol = (3*vol / (4*pi))^(1/3)
    R_area = (area / (4*pi))^(1/2)
    return R_vol, R_area
end

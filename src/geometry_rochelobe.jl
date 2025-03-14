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
    # Update the radii r(θ,ϕ) to match the surface potential
    npix = tessels.npix
    r = zeros(T, npix, 5);
    for ii = 1:npix
        for jj = 1:5
            r[ii,jj] = a*solve_radius(r_init, pot_surface, D, tessels.unit_spherical[ii,jj,2], tessels.unit_spherical[ii,jj,3], q, async_ratio, potential_function, verbose=false);
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
    potS = (pot_L1 + q * q / 2(1 + q)) / fillout_factor - q * q / 2(1 + q) # eq 6, Leahy paper
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
    secondary == false ? rtry = radius_point_Pathania(q) : radius_point_Pathania(1/q) 
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

# #
# # Gravitational potential equations from Aufdenberg et al. 2021 + SAGE math
# # "Modeling the Hα Emission Surrounding Spica Using the Lyman Continuum from a Gravity-darkened Central Star"
# # https://iopscience.iop.org/article/10.3847/1538-4357/ac1c0e
# # Primary at x,y,z = (0,0,0) and secondary at x,y,z=(D,0,0) 
# # D is the instantaneous unitless separation unitless (= separation_distance/semimajoraxis)
# # D = (1-e^2)/(1+e*cos.(υ)) υ: true anomaly (upsilon)
# Note: 
# Polar potential λ=0, ν = 1  note that here ν = cos(θ) (yeah, beware υ vs ν)
# Ω1 = 1/r + q/sqrt( D^2 + r^2)
# dΩ1 = -q*r/(D^2 + r^2)^(3/2) - 1/r^2
#
# Example of Sage setup to check derivatives
#
#var('r, D, λ, ν, q, async_ratio')
# f = 1/r + q/sqrt( D^2 + r^2 - 2*r*λ*D) - q*r*λ*D + 1/2*async_ratio^2*(1+q)*r^2*(1-ν^2)
# g = f.derivative(r)
# g.derivative(r)
# -(ν^2 - 1)*async_ratio^2*(q + 1) + 3*(D*λ - r)^2*q/(-2*D*r*λ + D^2 + r^2)^(5/2) - q/(-2*D*r*λ + D^2 + r^2)^(3/2) + 2/r^3

function compute_potential_primary(r, D, θ, ϕ, q, async_ratio) # r and D are dimensionless (were divided by a)
    # Note: async_ratio = ω1rot/ωorb
    λ = sin(θ)*cos(ϕ);
    ν = cos(θ) # WARNING, POTENTIAL SYMBOL CLASH WITH TRUE ANOMALY
    Ω1 = 1/r + q/sqrt( D^2 + r^2 - 2*r*λ*D) - q*r*λ*D + async_ratio^2*(1+q)*r^2*(1-ν^2)/2
    dΩ1 = -1/r^2 - q*(r-λ*D)/sqrt(( D^2 + r^2 - 2*r*λ*D)^3) - q*λ*D + async_ratio^2*(1+q)*r*(1-ν^2)
    ddΩ1 =  2/r^3  + 3*(D*λ - r)^2*q/sqrt((D^2 + r^2 - 2*r*λ*D)^5) - q/sqrt((D^2 + r^2 - 2*r*λ*D)^3) + async_ratio^2*(q + 1)*(1 - ν^2)
    return Ω1, dΩ1, ddΩ1
end

function compute_potential_secondary(r, D, θ, ϕ, q, async_ratio)
    # Note: async_ratio = ω2rot/ωorb
    λ = cos(ϕ)*sin(θ)
    ν = cos(θ)
    Ω2 = 1/sqrt(D^2+r^2+2*r*λ*D) + q/r + (1+q)/2*(D^2+2*D*r*λ) - q*(D^2+r*λ*D) + async_ratio^2*(1+q)*r^2*(1-ν^2)/2
    dΩ2 = - (D*λ + r)/sqrt((D^2+r^2+2*r*λ*D)^3) - q/r^2 + λ*D + async_ratio^2*(1+q)*r*(1-ν^2)
    ddΩ2 =  2*q/r^3 + 3*(D*λ + r)^2/sqrt((2*D*r*λ + D^2 + r^2)^5) - 1/sqrt((2*D*r*λ + D^2 + r^2)^3) + async_ratio^2*(1+q)*(1-ν^2)
    return Ω2, dΩ2, ddΩ2
end


function solve_radius(r0, pot_surface, D, θ, ϕ, q, async_ratio, potential_function; verbose=true)
    # ii=64; jj=1; star = stars[1]
    # r0 = r_init; θ = tessels.unit_spherical[ii,jj,2]; ϕ = tessels.unit_spherical[ii,jj,3]
    # r0 = 0.5923035715675565; θ = 0.0; ϕ = 0.0; pot_surface = 4.45647759841041; potential_function = compute_potential_primary
    # θ=1.5707964f0;ϕ = 357.1875f0
    fgh = r->potential_function(r, D, θ, ϕ, q, async_ratio);
    # close(); plot([fgh(r)[1] for r in range(0.0001, D, length=1000)])
    r = halley(r0, pot_surface, fgh, verbose=verbose);
    # fgh(r)[1]-pot_surface
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


function temperature_map_vonZeipel_roche_single(parameters, star_geom, t::Array{T,1};  secondary = false) where T
    npix = star_geom[1].npix
    nepochs = length(t)
    Tmap = zeros(T, npix, nepochs)
    for i=1:nepochs
        Tmap[:, i] = temperature_map_vonZeipel_roche_single(parameters, star_geom[i], t[i], secondary = false, T=T)
    end
    return Tmap
end

function temperature_map_vonZeipel_roche_single(parameters, star_geom, t; secondary = false, T=Float32)
    rpole = parameters.rpole/parameters.a
    tpole = parameters.tpole
    r = star_geom.vertices_spherical[:, 5 ,1]/parameters.a
    θ = star_geom.vertices_spherical[:, 5, 2]
    ϕ = star_geom.vertices_spherical[:, 5, 3]
    D = compute_separation(parameters, t)
    async_ratio = parameters.rotation_period/parameters.P
    q = parameters.q
    β = parameters.beta
    # Compute gravity
    compute_gravity = compute_gravity_primary;
    if secondary == true
        compute_gravity = compute_gravity_secondary;
    end
    g_pole = compute_gravity(rpole,T(0),T(0),D,q,async_ratio);
    gravity_map = compute_gravity(r,θ,ϕ,D,q,async_ratio);
    # Computes von Zeipel temperature map directly from the gravity map
    Teff = tpole*(gravity_map./g_pole).^β;
    return Teff
end

function temperature_map_vonZeipel_roche(binary_parameters, star_geom, t; secondary = false)
    stellar_parameters = binary_parameters.star1
    if secondary == true
        stellar_parameters = binary_parameters.star2;
    end
    rpole = stellar_parameters.rpole/binary_parameters.a
    tpole = stellar_parameters.tpole
    r = star_geom.vertices_spherical[:, 5 ,1]/binary_parameters.a
    θ = star_geom.vertices_spherical[:, 5, 2]
    ϕ = star_geom.vertices_spherical[:, 5, 3]
    D = compute_separation(binary_parameters, t)
    async_ratio = stellar_parameters.rotation_period/binary_parameters.P
    q = binary_parameters.q
    β = stellar_parameters.beta
    # Compute gravity
    compute_gravity = compute_gravity_primary;
    if secondary == true
        compute_gravity = compute_gravity_secondary;
    end
    g_pole = compute_gravity(rpole,0.0,0.0,D,q,async_ratio);
    gravity_map = compute_gravity(r,θ,ϕ,D,q,async_ratio);
     # Computes von Zeipel temperature map directly from the gravity map
    Teff = tpole*(gravity_map./g_pole).^β;
    return Teff
end


# Formula for Roche radii
# Radius of the Roche equipotential surfaces
# A. Pathania · T. Medupe
# Astrophys Space Sci (2012) 338:127–145
# DOI 10.1007/s10509-011-0928-y
# point radius (λ = 1, ν = 0) 
# back radius  (λ =−1, ν = 0) 
# the side radius (λ = 0, ν = 0)
# the pole radius (λ = 0, ν = 1)

# Eggleton formula for Roche equivalent size (eq 17)
function radius_equivalent_eggleton(q)
    q1=1/q
    return 0.49*q1^(2/3)/(0.6*q1^(2/3)+log(1.0+q1^(1/3)));
end

function radius_point_kopal_polynomial(x)
    # Polynomial = 0 at x=point radius
    return (q + 1)x^5 - (3q + 2)x^4 + (3q + 1)x^3 - x^2 + 2x - 1 
end

function radius_leahy(q) 
    #A calculator for Roche lobe properties, Denis A Leahy & Janet C Leahy   
    #Computational Astrophysics and Cosmology volume 2, Article number: 4 (2015)
    # Eq. 7
    a1 = 0.64334;
    a2 = 0.86907;
    a3 = 1.2809;
    a4 = −0.74303;
    a5 = 0.73103;
    return a1*q^a4/(a2*q^a5+log(1+a3*q^(a4+1/3)))
end

function radius_back_Pathania(q)
    return exp(-1.00598q^3 + 2.09674q^2 - 1.69263q - 0.319909)
end

function radius_point_Pathania(q)
    return exp(-0.725742q^3 + 1.53893q^2 - 1.31638q - 0.202505)
end



# function solve_radius_newton_raphson(r0, pot_surface, D, θ, ϕ, q, async_ratio, potential_function)
#     # Prone to convergence issues unless initialized at good location
#     # Example r0 = 0.5923035715675565; θ = 0.0; ϕ = 0.0; pot_surface = 4.45647759841041; potential_function = compute_potential_primary
#     # ---> will give NaN. Alternative is Interval Newton or Halley's method
#     fg = r->potential_function(r, D, θ, ϕ, q, async_ratio);
#     # close(); plot([fg(r)[1] for r in range(0.0001, D, length=1000)])
#     r = newton_raphson(r0, pot_surface, fg, verbose = true);
#     return r
# end

# function newton_raphson(x0, f0, fg; ϵ = 1e-5, verbose=false)
#     # Newton's method, trying to get x so that f(x) = f0 
#     n = 1; 
#     converged = false;
#     x = copy(x0)
#     while ((converged == false) & (n < 50))
#         f, g = fg(x);
#         newton_step = (f - f0)/g; 
#         if verbose == true
#             println("n = $n\t f = $f,\t f0 = $f0 \t g = $g \t x = $x \t step = $newton_step");
#         end  
#         x -= newton_step;
#         if (abs(newton_step) < ϵ)
#             converged = true;
#         end
#         n += 1;
#     end
#     return x
# end



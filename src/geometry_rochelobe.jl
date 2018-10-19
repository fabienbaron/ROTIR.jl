function latlong_rochelobe(ntheta,nphi,binary_parameters;secondary=false)
    separation = binary_parameters.separation;
    mass_ratio = binary_parameters.mass_ratio; async_ratio = binary_parameters.async_ratio;
    fillout_factor = binary_parameters.fillout_factor; eccentricity = binary_parameters.eccentricity;
    true_anomaly = binary_parameters.true_anomaly;

    # Do limit check for this calculator
    if ((fillout_factor .< 0.1) | (fillout_factor .> 1.0))
        println("Fill factor cannot be less than 0.1 or greater than 1.0.");
        return println("Try a new fill factor")
    elseif ((async_ratio .< 0.1) | (async_ratio .> 2.0))
        println("Rate of stellar rotation divided by binary rotation rate cannot")
        println("be less than 0.01 or greater than 2.")
        return
    end

    npix = ntheta*nphi;
    vertices_spherical = zeros(Float64, npix, 3, 5); # r, theta, phi
    vertices_xyz = zeros(Float64, npix, 3, 5); # x, y, z
    dphi = 2*pi/nphi;
    dtheta = pi/ntheta;

    # get corners
    # calculates theta and phi
    theta = collect(linspace(0,pi-dtheta,ntheta));
    phi = collect(linspace(0,2*pi-dphi,nphi));

    # make vertices -- go counterclockwise
    for i = 1:ntheta
        ilong_range = Int((i-1)*nphi+1):(i*nphi);

        vertices_spherical[ilong_range,2,1], vertices_spherical[ilong_range,3,1] = theta[i], phi;
        vertices_spherical[ilong_range,2,2], vertices_spherical[ilong_range,3,2] = theta[i] + dtheta, phi;
        vertices_spherical[ilong_range,2,3], vertices_spherical[ilong_range,3,3] = theta[i] + dtheta, phi + dphi;
        vertices_spherical[ilong_range,2,4], vertices_spherical[ilong_range,3,4] = theta[i], phi + dphi;
        theta_avg, phi_avg = theta[i] + dtheta*0.5, phi + dphi*0.5;
        vertices_spherical[ilong_range,2,5], vertices_spherical[ilong_range,3,5] = theta_avg, phi_avg;
    end

    if (secondary == true)
        mass_ratio = 1.0/mass_ratio;
    end

    # calculate radius at L1 point_ratio
    # Eggleton formula for approximating radius
    R_L1 = REg_f90(mass_ratio)*separation;
    #R_L1 = separation.*((mass_ratio./3.).^(1./3.)); # Only works if M2 >> M1
    R_L1 = solve_R_L1(R_L1,separation,mass_ratio,async_ratio);

    # Potential at the equator
    # theta & phi are in astro coordinates
    pot_L1, dpot_L1 = compute_potential(R_L1,separation,pi/2.,0.,mass_ratio,async_ratio);

    # New potential given fill out factor
    pot_F = (pot_L1 + 0.5*mass_ratio*mass_ratio/(1. + mass_ratio))/fillout_factor - 0.5*mass_ratio*mass_ratio/(1. + mass_ratio);

    R_L1_f90 = R_L1/separation;

    # calculate radius for roche lobe
    # TBD: Try parallel computing
    for i = 1:npix
        for j = 1:5
            ##vertices_spherical[i,1,j] = compute_radius(R_L1,pot_L1,separation,vertices_spherical[i,2,j],vertices_spherical[i,3,j],mass_ratio,async_ratio);
            #theta_temp = vertices_spherical[i,2,j];
            #phi_temp = vertices_spherical[i,3,j];
            #vertices_spherical[i,1,j] = compute_radius(R_L1,pot_F,separation,theta_temp,phi_temp,mass_ratio,async_ratio,fillout_factor);

            # Multiply by separation to get true radius in milliarcseconds
            vertices_spherical[i,1,j] = getpotroot_f90(mass_ratio,async_ratio,fillout_factor,R_L1_f90,pot_F,
                vertices_spherical[i,2,j],vertices_spherical[i,3,j],eccentricity,true_anomaly)*separation;
        end
    end

    vertices_xyz[:,1,:] = vertices_spherical[:,1,:].*sin.(vertices_spherical[:,2,:]).*cos.(vertices_spherical[:,3,:]); # X
    vertices_xyz[:,2,:] = vertices_spherical[:,1,:].*sin.(vertices_spherical[:,2,:]).*sin.(vertices_spherical[:,3,:]); # Y
    vertices_xyz[:,3,:] = vertices_spherical[:,1,:].*cos.(vertices_spherical[:,2,:]); # Z

    star_base_geom = base_geometry(npix,vertices_xyz, vertices_spherical);
end

# The following functions were adapted from SIMTOI
# https://github.com/bkloppenborg/simtoi
#
# and equations from Leahy & Leahy 2015
# https://link.springer.com/article/10.1186/s40668-015-0008-8


# Initial programming started with Fortran/C++ mentality
# TBD: Once confirmed working, switch to Julia syntax (faster?)
function compute_potential(radius_star,separation,theta,phi,mass_ratio,async_ratio)
    # makes the radius dimensionless
    radius = radius_star/separation;

    # computes the dimensionless potential in the case of non-syncronous rotation
    # (Leahy & Leahy 2015 -- equation 2)
    lambda = sin(theta)*cos(phi);
    sin_squared = sin(theta)*sin(theta);
    sqrt_radius = sqrt(1.0 - 2.0*radius*lambda + radius*radius);

    potential = 1.0/radius + mass_ratio/sqrt_radius - mass_ratio*radius*lambda +
        (mass_ratio + 1.0)*async_ratio*async_ratio*radius*radius*sin_squared/2.;

    # calculate deriviative of potential for Newton's method
    # adapted from SIMTOI
    dpotential = -1.0/(radius*radius) + mass_ratio*(lambda - radius)/(sqrt_radius*sqrt_radius*sqrt_radius) -
        mass_ratio*lambda + (mass_ratio + 1.0)*async_ratio*async_ratio*radius*sin_squared;

    return potential, dpotential
end

function compute_radius(R_init,pot_F,separation,theta,phi,mass_ratio,async_ratio)
    #R_start = deepcopy(R_init);
    radius = R_init/separation;
    # set up integers for guessing
    epsilon = 1.e-5; i = 1; converged = false;

    # Newton's method
    while ((converged == false) & (i <= 500))
        pot, dpot = compute_potential(R_init,separation,theta,phi,mass_ratio,async_ratio);
        println("potential = $pot, dpotential = $dpot");
        newton_step = pot - pot_F; # newton step
        radius = R_init/separation;
        radius = radius - newton_step/dpot;
        R_init = radius*separation;

        println("radius = $radius, i = $i, newton step = $newton_step");
        println(" ");

        if (abs(newton_step) < epsilon)
            converged = true;
        end
        i += 1;
    end

    return radius*separation
end

function solve_R_L1(R_init,separation,mass_ratio,async_ratio)
    # unitless radius
    radius = R_init/separation;
     # set up integers for guessing
    epsilon = 1.e-5; i = 1; converged = false;

    # Newton's method
    while ((converged == false) & (i <= 500))
        # first derivative of potential
        dpot_term1 = -1.0/(radius*radius) + mass_ratio./((1.0 - radius)*(1.0 - radius));
        dpot_term2 = -mass_ratio;
        dpot_term3 = (mass_ratio + 1.0)*async_ratio*async_ratio*radius;
        dpot = dpot_term1 + dpot_term2 + dpot_term3;

        # second deriviative of the potential
        ddpot = 2.0/(radius.*radius.*radius) + mass_ratio.*(2.0/((1.0-radius).*(1.0-radius).*(1.0-radius))) +
            (mass_ratio + 1.0)*async_ratio*async_ratio;

        radius = radius - dpot./ddpot;

        #println("radius = $radius, dpotential = $dpot, i = $i");
        #println(" ")

        if (abs(dpot) < epsilon)
            converged = true;
        end
        i += 1;
    end
    return radius*separation
end

function compute_gravity(radius_star,separation,theta,phi,async_ratio,mass_ratio)
    # makes the radius dimensionless
    radius = radius_star./separation;

    lambda = sin.(theta).*cos.(phi);
    mu = sin.(theta).*sin.(theta);
    nu = sin.(phi).*sin.(theta);
    sqrt_radius = sqrt.(1.0 - 2.0*radius.*lambda + radius.*radius);
    cube_sqrt_radius = sqrt_radius.*sqrt_radius.*sqrt_radius;

    gx = -1.0/(radius.*radius) + mass_ratio*(1.0 - radius.*lambda)./cube_sqrt_radius + async_ratio*async_ratio*(1.0 + mass_ratio).*radius.*lambda - mass_ratio;
    gy = -1.0/(radius.*radius) + mass_ratio.*radius.*nu./cube_sqrt_radius + async_ratio*async_ratio*(1.0 + mass_ratio).*radius.*nu;
    gz = -1.0/(radius.*radius) - mass_ratio.*radius.*cos.(theta)./cube_sqrt_radius;

    gravity = sqrt.(gx.*gx + gy.*gy + gz.*gz);
    return gravity
end

function compute_teff_vonzeipel(R_pole,Teff_pole,radius_star,separation,theta,phi,async_ratio,mass_ratio,beta)
    # computes von Zeipel temperatures directly from gravity
    g_pole = compute_gravity(R_pole,separation,0.0,0.0,async_ratio,mass_ratio);
    gravity = compute_gravity(radius_star,separation,theta,phi,async_ratio,mass_ratio);

    Teff = Teff_pole*((gravity./g_pole).^beta);
    return Teff
end

function getpotroot_f90(mass_ratio,async_ratio,fillout_factor,R_L1,pot_F,theta,phi,eccentricity,true_anomaly)
    #near L1 point the Omega-OmegaF is tangent to zero so rootfinder doesn't work
    #in this case just use R at angle alim off the x axis (set alim=1degree).
    #if F.gt.Flim .and. theta and phi within angle alim of x axis,
    #then use R(theta,phi)=RL1*(1-ang/alim)+R(alim)*ang/alim with ang=angle from x-axis
    #in this case need to find the root for thet2=pi/2.d0;ph2=alim to get R(alim)

    Flim = 0.997;
    alim=0.1*pi/180.;

    # temporary fix for theta = pi/2 and phi = 2*pi by AOM
    # Note: no error at theta = pi/2 and phi = 0.
    if ((theta .>= (pi/2.0-alim)) & (theta .<= (pi/2.0+alim)) & (phi .>= (2.0*pi-alim)))
        phi = 2.0*pi - phi;
    end

    ang2 = phi*phi+(theta-pi/2.0)*(theta-pi/2.0);
    ang = sqrt(ang2);

    if ((fillout_factor .> Flim) & (ang .< alim))
        thet2 = pi/2.0; ph2 = deepcopy(alim);
        flag2 = true;
    else
        thet2 = deepcopy(theta); ph2 = deepcopy(phi);
        flag2 = false;
    end

    # convergence criterion in x
    xcon = 1.e-9;
    # initial guess
    a = 0.05*REg_f90(mass_ratio);
    b = 1.0001*R_L1;
    fa = potential_f90(a,thet2,ph2,mass_ratio,async_ratio)-pot_F;
    fb = potential_f90(b,thet2,ph2,mass_ratio,async_ratio)-pot_F;
    #fa = potential_ellipse_f90(a,thet2,ph2,mass_ratio,async_ratio,eccentricity,true_anomaly)-pot_F;
    #fb = potential_ellipse_f90(b,thet2,ph2,mass_ratio,async_ratio,eccentricity,true_anomaly)-pot_F;

    # test that root is bracketed by a,b
    if (fa*fb .> 0.)
        println("a = $a,b = $b,fa = $fa,fb = $fb");
        println("bad initial guess in potential root")#, i = $i, j = $j");
        println("theta = $thet2, phi = $ph2")
        #return
    end

    # swap a,b if needed
    if (abs(fa) .< abs(fb))
        dum = deepcopy(a); a = deepcopy(b); b = deepcopy(dum)
        dum = deepcopy(fa); fa = deepcopy(fb); fb = deepcopy(dum)
        #write(6,*) "swap: a,b,fa,fb",a,b,fa,fb
    end

    c = deepcopy(a); d = 0.; # AOM added d=c
    fc = potential_f90(c,thet2,ph2,mass_ratio,async_ratio)-pot_F;
    #fc = potential_ellipse_f90(c,thet2,ph2,mass_ratio,async_ratio,eccentricity,true_anomaly)-pot_F;
    #mflag = 1
    mflag = true;

    s = 0.; # AOM
    while ((abs(b-a) .> xcon) & (fb .!= 0.0))   #TODO: change float comparison 
        #println("while loop")
        if ((fa .!= fc) & (fb .!= fc))
            s = a*fb*fc/((fa-fb)*(fa-fc))+b*fa*fc/((fb-fa)*(fb-fc))+c*fa*fb/((fc-fa)*(fc-fb));
        else
            s = b-fb*(b-a)/(fb-fa)
        end
        con1 = (s .< (3.0*a+b)/4.0) .| (s .> b)
        con2 = mflag & (abs(s-b) .> abs(b-c)/2.)
        con3 = (~mflag) & (abs(s-b) .> abs(c-d)/2.)
        con4 = mflag & (abs(b-c) .< abs(xcon))
        con5 = (~mflag) & (abs(c-d) .< abs(xcon))
        if (con1 .| con2 .| con3 .| con4 .| con5)
            s=(a+b)/2.;
            mflag=true;
        else
            mflag=false;
        end

        fs = potential_f90(s,thet2,ph2,mass_ratio,async_ratio)-pot_F;
        #fs = potential_ellipse_f90(s,thet2,ph2,mass_ratio,async_ratio,eccentricity,true_anomaly)-pot_F;
        d=deepcopy(c);
        c=deepcopy(b);
        fc=potential_f90(c,thet2,ph2,mass_ratio,async_ratio)-pot_F;
        #fc = potential_ellipse_f90(c,thet2,ph2,mass_ratio,async_ratio,eccentricity,true_anomaly)-pot_F;

        if (fa*fs .< 0.)
            b=deepcopy(s); fb=deepcopy(fs);
        else
            a=deepcopy(s); fa=deepcopy(fs)
        end

        #swap a,b if needed
        if (abs(fa) .< abs(fb))
            dum = deepcopy(a); a = deepcopy(b); b = deepcopy(dum);
            dum = deepcopy(fa); fa = deepcopy(fb); fb = deepcopy(dum);
        end
    end

    # for case ang<alim then need to calculate R from RL1 and s
    if (flag2)
        # linear interpolate between s and RL1
        xroot=R_L1*(1.0-ang/alim)+s*ang/alim;
    else
        xroot=deepcopy(s);
    end

    return xroot
end

# the dimensionless Roche Potential
function potential_f90(r,theta,phi,q,P)
    dcp=cos(phi);dst=sin(theta)
    term1 = 1.0/r + q/sqrt(1.0-2.0*r*dcp*dst+r*r)
    term2 = -q*r*dcp*dst + (q+1.0)/2.0*P*P*r*r*dst*dst
    potential = term1+term2
    return potential
end

# the dimensionless Roche Potential modified by AOM for elliptical orbit
function potential_ellipse_f90(r,theta,phi,q,P,e,nu)
    dcp=cos(phi);dst=sin(theta);
    eccen_4 = (1.0 + e)^4;
    e_nu_3 = (1.0 + e*cos(nu))^3;
    term1 = 1.0/r + q/sqrt(1.0-2.0*r*dcp*dst+r*r);
    term2 = -q*r*dcp*dst + (q+1.0)/2.0*P*P*r*r*dst*dst*eccen_4/e_nu_3;
    potential = term1+term2
    return potential
end

# Eggleton formula
function REg_f90(q)
    REg = 0.49/(0.6 + q^(2.0/3.0)*log(1.0+q^(-1.0/3.0)));
    return REg
end

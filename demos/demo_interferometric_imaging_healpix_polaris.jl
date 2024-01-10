# using Distances
include("../src/ROTIR.jl"); using Main.ROTIR;
using OITOOLS
# LOAD DATA
oifitsfile = "../../OITOOLS.jl/demos/data/vscale_oiprep_2021Apr02_03_04_MIRCX_HD_8890_AVG15m_all.fits"
nepochs = 1
tepochs = [0.0]
data = [readoifits(oifitsfile, filter_bad_data=true, use_vis=false )[1,1]];

# SETUP STAR MODEL PARAMETERS
stellar_parameters = Array{starparameters}(undef, nepochs);
starparams = [3.18/2,  # milliarcseconds (radius at pole)
     7200.,              # Kelvin (at pole)
     0.,                 # unitless; fractional rotational velocity
     [1, 0.0, 0.0],        # limb darkening,first coefficient is for LD law type, then LD coefficients
     0.0,               # exponent for von Zeipel law
     0.,           # 2nd constant for rotational velocity
     90.,            # degrees; inclination
     0.,            # degrees; position_angle
     0.];           # days; rotation_period

for i=1:nepochs
    stellar_parameters[i]=starparameters(starparams[1],starparams[2],starparams[3],starparams[4],starparams[5],
        starparams[6],starparams[7],starparams[8],0.0,starparams[9]);
end

# SETUP 3D GEOMETRY (HEALPIX)
n=4; star_epoch_geom = create_geometry( healpix_round_star(n,radius=stellar_parameters[1].radius), stellar_parameters);
polyflux, polyft = setup_polygon_ft(data, star_epoch_geom);
x_start = 7200*ones(star_epoch_geom[1].npix);
tvinfo = tv_neighbours_healpix(n);
regularizers = [["tv", 0.01, tvinfo,1:length(x_start)]];
x =  spheroid_oi_reconstruct(x_start, data, polyflux, polyft, regularizers = regularizers, verb = true, maxiter=250);
plot2d_temperature(x, star_epoch_geom[1], plotmesh=false); # Note: Temp = intensity since no LDD

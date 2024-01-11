include("../src/ROTIR.jl"); using Main.ROTIR;
using OITOOLS
oifitsfile = "./data/RWCep_Dec2022.oifits"
oifitsfile = "./data/RWCep_Jul2023.oifits"

nepochs = 1
tepochs = [0.0]
data = [readoifits(oifitsfile, filter_bad_data=true, use_vis=false )[1,1]];

# Determine diameter to use
# 
weights=[1.0,1.0,1.0]
disc = create_component(type="ldpow", name="disc");
disc.vis_params[1].val=2.0
disc.vis_params[1].minval=1.0
disc.vis_params[1].maxval=5.0
model=create_model(disc)
minf, minx, cvis_model, result = fit_model_ultranest(data[1], model, weights=weights); #interesting local minimum with v2 only!
minf, minx, cvis_model, result = fit_model_nlopt(data[1], model, weights=weights);
diameter = minx[1]*1.05
ldpow = minx[2]

# SETUP STAR MODEL PARAMETERS
stellar_parameters = Array{starparameters}(undef, nepochs);
starparams = [diameter/2,  # milliarcseconds (radius at pole)
     4200.,              # Kelvin (at pole)
     0.,                 # unitless; fractional rotational velocity
   # [1,0.0,0.0],
     [3, ldpow],        # limb darkening,first coefficient is for LD law type, then LD coefficients
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
n=3; star_epoch_geom = create_geometry( healpix_round_star(n,radius=stellar_parameters[1].radius), stellar_parameters);
polyflux, polyft = setup_polygon_ft(data, star_epoch_geom);
x_start = 4200*ones(star_epoch_geom[1].npix);
tvinfo = tv_neighbours_healpix(n);
regularizers = [["tv", 0.05, tvinfo,1:length(x_start)]];
regularizers = [["l2", 0.01, tvinfo, 1:length(x_start)]];
#regularizers = [["tv2", 5e-4, tvinfo,1:length(x_start)]];
#x =  spheroid_oi_reconstruct(x_start, data, polyflux, polyft, lower=3000, upper=7000, regularizers = regularizers, verb = true, maxiter=500);
x =  spheroid_oi_reconstruct(x_start, data, polyflux, polyft, regularizers = regularizers, verb = true, maxiter=500);
plot2d_temperature(x, star_epoch_geom[1], plotmesh=false); # Note: Temp = intensity since no LDD

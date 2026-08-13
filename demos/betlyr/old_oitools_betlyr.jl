using OITOOLS
oifitsfiles = ["MIRC_L2.bet_Lyr.2013Jun20_XYZ_2015Oct15.XCHAN.XYZ.AVG15m.oifits",
"MIRC_L2.bet_Lyr.2013Jun21_XYZ_2015Oct14.XCHAN.XYZ.AVG15m.oifits",
"MIRC_L2.bet_Lyr.2013Jun22_XYZ_2015Oct14.XCHAN.XYZ.AVG15m.oifits",
"MIRC_L2.bet_Lyr.2013Jun23_XYZ_2015Oct14.XCHAN.XYZ.AVG15m.oifits",
"MIRC_L2.bet_Lyr.2013Jun24_XYZ_2015Oct16.XCHAN.XYZ.AVG15m.oifits",
"MIRC_L2.bet_Lyr.2013Jun26_XYZ_2015Oct15.XCHAN.XYZ.AVG15m.oifits",
"MIRC_L2.bet_Lyr.2013Jun27_XYZ_2015Oct16.XCHAN.XYZ.AVG15m.oifits",
"MIRC_L2.bet_Lyr.2013Jun28_XYZ_2015Oct20.XCHAN.XYZ.AVG15m.oifits",
"MIRC_L2.bet_Lyr.2013Jun29_tcoh0017ms_FB_2017Apr03.XCHAN.FB.AVG15m.oifits",
"MIRC_L2.bet_Lyr.2013Jun30_tcoh0017ms_FB_2017Apr03.XCHAN.FB.AVG15m.oifits",
"MIRC_L2.bet_Lyr.2013Jul01_tcoh0017ms_XYZ_2017Apr03.XCHAN.FB.AVG15m.oifits"]


titles = ["2013Jun20",
"2013Jun21",
"2013Jun22",
"2013Jun23",
"2013Jun24",
"2013Jun26",
"2013Jun27",
"2013Jun28",
"2013Jun29",
"2013Jun30",
"2013Jul01"]

pixsize = 0.05
nx = 64

# # OITOOLS
for i=1:length(oifitsfiles)
data = readoifits(oifitsfiles[i])[1,1];
ft = setup_nfft(data, nx, pixsize);
#initial image is a simple Gaussian
x_start = gaussian2d(nx,nx,nx/6);
x_start = vec(x_start)/sum(x_start);
regularizers = [["centering", 1e4], ["tv", 7e3]];
x = reconstruct(x_start, data, ft, regularizers = regularizers, verb = true, maxiter=1000);
x = reconstruct(x, data, ft, regularizers = regularizers, verb = true, maxiter=1000);
x = reconstruct(x, data, ft, regularizers = regularizers, verb = true, maxiter=1000);
imdisp(x,pixscale=pixsize)
writefits(reshape(x,nx,nx),"reconstruction_oitools_$i.fits")
end

# Squeeze
for i=1:length(oifitsfiles)
run(`/home/baron/SOFTWARE/squeeze/bin/squeeze $(oifitsfiles[i]) -w 64 -s 0.05 -tv 500 -e 1500 -o reconstruction_squeeze_$i.fits`)
img = readfits("reconstruction_squeeze_$i.fits")[:,:,1]
maxpos = findmax(img)[2]
img = circshift(img, div(nx+1,2) .- Tuple(maxpos))
imdisp(img,pixscale=pixsize);
sleep(0.6)
end

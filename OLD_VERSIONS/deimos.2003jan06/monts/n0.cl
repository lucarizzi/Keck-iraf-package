set dfile="/net/polo/dsk/iraid1/kics/2002jun28/d0628_00"
set mfile="/net/polo/dsk/iraid1/deimos2/monts/f0628_"
qmodel.mu = -18.858
qmodel.o3 = 0.014
qmodel.roll3 = 0.058
qmodel.gmm = 831.
qmodel.norder = 0
qmodel.amap = "amap"
qmodel.bmap = "bmap"

mkmontage.border = 0
mkmontage.background = 1.
mkmontage.fill = -1.
mkmontage.nx = 40
mkmontage.ny = 40
mkmontage.mk_image = yes
mkmontage.nice_form = yes
mkmontage.ncols = 240
mkmontage.nlines = 160
mkmontage.auto_scale = yes
mkmontage.f1 = 0.0
mkmontage.f2 = 1.0
mkmontage.normflux = 1.
mkmontage.verbose = yes

del tmp.coo
qmod goh.6x4 wave=5000  sc=2.0e-3 >> tmp.coo

imdel mfile$13n.fits
mkmont dfile$13 tmp.coo mfile$13n.fits title="test w/ rev qmodel"


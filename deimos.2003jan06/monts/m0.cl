set dfile="/net/polo/dsk/iraid1/kics/2002jun27/d0627_00"
qmodel.mu = 0.376
qmodel.o3 = 0.014
qmodel.roll3 = 0.067
qmodel.gmm = 831.
qmodel.norder = 1
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
mkmontage.nlines = 200
mkmontage.auto_scale = yes
mkmontage.f1 = 0.0
mkmontage.f2 = 1.0
mkmontage.normflux = 1.
mkmontage.verbose = yes

del tmp.coo
qmod loh.six wave=9224.50 sc=0.4e-3 >> tmp.coo
qmod loh.six wave=8521.44 sc=0.4e-3 >> tmp.coo
qmod loh.six wave=7948.18 sc=0.4e-3 >> tmp.coo
qmod loh.six wave=6965.43 sc=0.4e-3 >> tmp.coo
qmod loh.six wave=6266.50 sc=0.4e-3 >> tmp.coo
ty tmp.coo

imdel foc19.fits
mkmont dfile$19 tmp.coo foc19.fits title="test"



del tmp.coo
qmod loh.six wave=9224.50 sc=2.8e-3 >> tmp.coo
qmod loh.six wave=8521.44 sc=2.8e-3 >> tmp.coo
qmod loh.six wave=7948.18 sc=2.8e-3 >> tmp.coo
qmod loh.six wave=6965.43 sc=2.8e-3 >> tmp.coo
qmod loh.six wave=6266.50 sc=2.8e-3 >> tmp.coo
ty tmp.coo

imdel foc26.fits
mkmont dfile$26 tmp.coo foc26.fits title="test"

# imdel foc41.fits
# mkmont dfile$41 tmp.coo foc41.fits title="test"

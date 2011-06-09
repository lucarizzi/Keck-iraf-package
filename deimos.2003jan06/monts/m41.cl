set dfile="/net/polo/dsk/iraid1/kics/2002jun27/d0627_00"
set mfile="/net/polo/dsk/iraid1/deimos2/monts/f0627_"
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

imdel mfile$39.fits
mkmont dfile$39 tmp.coo mfile$39.fits title="test"



del tmp.coo
qmod loh.six wave=9224.50 sc=0.8e-3 >> tmp.coo
qmod loh.six wave=8521.44 sc=0.8e-3 >> tmp.coo
qmod loh.six wave=7948.18 sc=0.8e-3 >> tmp.coo
qmod loh.six wave=6965.43 sc=0.8e-3 >> tmp.coo
qmod loh.six wave=6266.50 sc=0.8e-3 >> tmp.coo
ty tmp.coo

imdel mfile$40.fits
mkmont dfile$40 tmp.coo mfile$40.fits title="test"


del tmp.coo
qmod loh.six wave=9224.50 sc=1.1e-3 >> tmp.coo
qmod loh.six wave=8521.44 sc=1.1e-3 >> tmp.coo
qmod loh.six wave=7948.18 sc=1.1e-3 >> tmp.coo
qmod loh.six wave=6965.43 sc=1.1e-3 >> tmp.coo
qmod loh.six wave=6266.50 sc=1.1e-3 >> tmp.coo
ty tmp.coo

imdel mfile$41.fits
mkmont dfile$41 tmp.coo mfile$41.fits title="test"


del tmp.coo
qmod loh.six wave=9224.50 sc=1.5e-3 >> tmp.coo
qmod loh.six wave=8521.44 sc=1.5e-3 >> tmp.coo
qmod loh.six wave=7948.18 sc=1.5e-3 >> tmp.coo
qmod loh.six wave=6965.43 sc=1.5e-3 >> tmp.coo
qmod loh.six wave=6266.50 sc=1.5e-3 >> tmp.coo
ty tmp.coo

imdel mfile$42.fits
mkmont dfile$42 tmp.coo mfile$42.fits title="test"


del tmp.coo
qmod loh.six wave=9224.50 sc=1.8e-3 >> tmp.coo
qmod loh.six wave=8521.44 sc=1.8e-3 >> tmp.coo
qmod loh.six wave=7948.18 sc=1.8e-3 >> tmp.coo
qmod loh.six wave=6965.43 sc=1.8e-3 >> tmp.coo
qmod loh.six wave=6266.50 sc=1.8e-3 >> tmp.coo
ty tmp.coo

imdel mfile$43.fits
mkmont dfile$43 tmp.coo mfile$43.fits title="test"


del tmp.coo
qmod loh.six wave=9224.50 sc=2.2e-3 >> tmp.coo
qmod loh.six wave=8521.44 sc=2.2e-3 >> tmp.coo
qmod loh.six wave=7948.18 sc=2.2e-3 >> tmp.coo
qmod loh.six wave=6965.43 sc=2.2e-3 >> tmp.coo
qmod loh.six wave=6266.50 sc=2.2e-3 >> tmp.coo
ty tmp.coo

imdel mfile$44.fits
mkmont dfile$44 tmp.coo mfile$44.fits title="test"

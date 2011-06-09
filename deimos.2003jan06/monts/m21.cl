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

imdel mfile$19.fits
mkmont dfile$19 tmp.coo mfile$19.fits title="test"



del tmp.coo
qmod loh.six wave=9224.50 sc=0.8e-3 >> tmp.coo
qmod loh.six wave=8521.44 sc=0.8e-3 >> tmp.coo
qmod loh.six wave=7948.18 sc=0.8e-3 >> tmp.coo
qmod loh.six wave=6965.43 sc=0.8e-3 >> tmp.coo
qmod loh.six wave=6266.50 sc=0.8e-3 >> tmp.coo
ty tmp.coo

imdel mfile$20.fits
mkmont dfile$20 tmp.coo mfile$20.fits title="test"


del tmp.coo
qmod loh.six wave=9224.50 sc=1.1e-3 >> tmp.coo
qmod loh.six wave=8521.44 sc=1.1e-3 >> tmp.coo
qmod loh.six wave=7948.18 sc=1.1e-3 >> tmp.coo
qmod loh.six wave=6965.43 sc=1.1e-3 >> tmp.coo
qmod loh.six wave=6266.50 sc=1.1e-3 >> tmp.coo
ty tmp.coo

imdel mfile$21.fits
mkmont dfile$21 tmp.coo mfile$21.fits title="test"


del tmp.coo
qmod loh.six wave=9224.50 sc=1.5e-3 >> tmp.coo
qmod loh.six wave=8521.44 sc=1.5e-3 >> tmp.coo
qmod loh.six wave=7948.18 sc=1.5e-3 >> tmp.coo
qmod loh.six wave=6965.43 sc=1.5e-3 >> tmp.coo
qmod loh.six wave=6266.50 sc=1.5e-3 >> tmp.coo
ty tmp.coo

imdel mfile$22.fits
mkmont dfile$22 tmp.coo mfile$22.fits title="test"


del tmp.coo
qmod loh.six wave=9224.50 sc=1.8e-3 >> tmp.coo
qmod loh.six wave=8521.44 sc=1.8e-3 >> tmp.coo
qmod loh.six wave=7948.18 sc=1.8e-3 >> tmp.coo
qmod loh.six wave=6965.43 sc=1.8e-3 >> tmp.coo
qmod loh.six wave=6266.50 sc=1.8e-3 >> tmp.coo
ty tmp.coo

imdel mfile$23.fits
mkmont dfile$23 tmp.coo mfile$23.fits title="test"


del tmp.coo
qmod loh.six wave=9224.50 sc=2.2e-3 >> tmp.coo
qmod loh.six wave=8521.44 sc=2.2e-3 >> tmp.coo
qmod loh.six wave=7948.18 sc=2.2e-3 >> tmp.coo
qmod loh.six wave=6965.43 sc=2.2e-3 >> tmp.coo
qmod loh.six wave=6266.50 sc=2.2e-3 >> tmp.coo
ty tmp.coo

imdel mfile$24.fits
mkmont dfile$24 tmp.coo mfile$24.fits title="test"


del tmp.coo
qmod loh.six wave=9224.50 sc=2.5e-3 >> tmp.coo
qmod loh.six wave=8521.44 sc=2.5e-3 >> tmp.coo
qmod loh.six wave=7948.18 sc=2.5e-3 >> tmp.coo
qmod loh.six wave=6965.43 sc=2.5e-3 >> tmp.coo
qmod loh.six wave=6266.50 sc=2.5e-3 >> tmp.coo
ty tmp.coo

imdel mfile$25.fits
mkmont dfile$25 tmp.coo mfile$25.fits title="test"


del tmp.coo
qmod loh.six wave=9224.50 sc=2.8e-3 >> tmp.coo
qmod loh.six wave=8521.44 sc=2.8e-3 >> tmp.coo
qmod loh.six wave=7948.18 sc=2.8e-3 >> tmp.coo
qmod loh.six wave=6965.43 sc=2.8e-3 >> tmp.coo
qmod loh.six wave=6266.50 sc=2.8e-3 >> tmp.coo
ty tmp.coo

imdel mfile$26.fits
mkmont dfile$26 tmp.coo mfile$26.fits title="test"

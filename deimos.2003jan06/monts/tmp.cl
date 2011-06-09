set dfile="/net/polo/dsk/iraid1/kics/2002jun27/d0627_00"
del tmp.coo
qmod loh.six wave=9224.50 sc=2.5e-3 >> tmp.coo
qmod loh.six wave=8521.44 sc=2.5e-3 >> tmp.coo
qmod loh.six wave=7948.18 sc=2.5e-3 >> tmp.coo
qmod loh.six wave=6965.43 sc=2.5e-3 >> tmp.coo
qmod loh.six wave=6266.50 sc=2.5e-3 >> tmp.coo
ty tmp.coo

imdel foc25.fits
mkmont dfile$25 tmp.coo foc25.fits title="test"


del tmp.coo
qmod loh.six wave=9224.50 sc=2.8e-3 >> tmp.coo
qmod loh.six wave=8521.44 sc=2.8e-3 >> tmp.coo
qmod loh.six wave=7948.18 sc=2.8e-3 >> tmp.coo
qmod loh.six wave=6965.43 sc=2.8e-3 >> tmp.coo
qmod loh.six wave=6266.50 sc=2.8e-3 >> tmp.coo
ty tmp.coo

imdel foc26.fits
mkmont dfile$26 tmp.coo foc26.fits title="test"
###################################################################

# imdel foc41.fits
# mkmont dfile$41 tmp.coo foc41.fits title="test"

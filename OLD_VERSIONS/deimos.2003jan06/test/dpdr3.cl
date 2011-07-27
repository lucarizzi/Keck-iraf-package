# THIS IS FOR COLLIMATOR TILTS
del dpdr3.dat
trace pdr.inp  t3=11.67113     gmm=1200 nord=1  > dpdr3.dat
trace pdr.inp  t3=11.67113 p1=90 t1=0.006944      gmm=1200 nord=1 >> dpdr3.dat
trace pdr.inp  t3=11.67113 p1=0  t1=0.006944      gmm=1200 nord=1 >> dpdr3.dat

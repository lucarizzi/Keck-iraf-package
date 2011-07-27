#
sim marc0_obj.sim  marc0_obj.fab  ra0=14.27622 dec0=52.15223 PA0=0 min_s=0 lamb=7500. slit=0.75 plot="marc0_obj.mng"  maskid="Marc0_obj"  surf+
print "Checking marc0_obj:"
!diff marc0_obj.fab ../test/testo.fab
## NB, NOTE RENAME above
!mongo < marc0_obj.mng
#
sim marc0_slit.sim marc0_slit.fab ra0=14.27622 dec0=52.15223 PA0=0 min_s=0 lamb=7500. slit=0.75 mdf=testo.fab plot="marc0_slit.mng" maskid="Marc0_slit" surf+
print "Checking marc0_slit:"
!diff marc0_slit.fab ../marc/test1.fab
## NB, NOTE RENAME above!!
!mongo < marc0_slit.mng
#

sim poh1.sim poh1.fab plot="poh1.mng" maskid="POH_1" surf+
print "Checking poh1:"
!diff poh1.fab ../test/poh1.fab
!mongo < poh1.mng
#
sim gohsp2.sim gohsp2.fab plot="gohsp2.mng" maskid="GOHSP_2" surf+
print "Checking gohsp2:"
!diff gohsp2.fab ../test/gohsp2.fab
!mongo < gohsp2.mng
#
# sim goh1.sim goh1.fab plot="goh1.mng" maskid="GOH_1" surf+
# print "Checking goh1:"	# REQUIRES INTERNAL CHANGES TO SIMTEST
# !diff goh1.fab ../test/goh1.fab
# !mongo < goh1.mng
#
sim goh2.sim goh2.fab plot="goh2.mng" maskid="GOH_2" surf+
print "Checking goh2:"
!diff goh2.fab ../test/goh2.fab
!mongo < goh2.mng
#
sim losh1.sim losh1.fab plot="losh1.mng" maskid="LOSH_1" surf+
print "Checking losh1:"
!diff losh1.fab ../test/loh1.fab
## NB, NOTE RENAME above
!mongo < losh1.mng
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

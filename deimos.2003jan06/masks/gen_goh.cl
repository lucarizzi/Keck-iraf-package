# GRIDS REQUIRE SPECIAL S/W
# sim goh1.sim goh0.fab plot="goh0.mng" maskid="GOH_0" surf+
# print "Checking goh0:"	# REQUIRES INTERNAL CHANGES TO SIMTEST
# !mongo < goh0.mng
#
sim goh1.sim goh1.fab plot="goh1.mng" maskid="GOH_1" surf+
print "Checking goh1:"	# REQUIRES INTERNAL CHANGES TO SIMTEST
!diff goh1.fab ../test/goh1.fab
!mongo < goh1.mng
#
## NB: I see diffs of up to 0.0002 mm between this and the original.

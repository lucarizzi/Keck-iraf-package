# lowd_ir_ohlines.dat / GDW / 2001-Jan-08
# 
# Purpose:
#	Remove insignificant OH lines from linelist, and average
#	together unresolved pairs on the assumption that the line strengths
#	will be equal.
# 
# Author:
#	Gregory D. Wirth, W. M. Keck Observatory
#-----------------------------------------------------------------------
BEGIN{	last_intensity = -1 ;
	min_intensity = 10 ;
	print "# NOTE: this version of the linelist has been modified"
	print "# to (a) exclude all lines with intensity less than", \
		min_intensity
	print "# and (b) combine all P1/P2 pairs of unresolved lines."
	print "# Combined lines are indentified below as 'blend' and"
	print "# the listed wavelengths are the simple average of the"
	print "# two components.  The command used to modify the"
	print "# linelist was:"
	print "#	awk -f lowd_ir_ohlines.awk ir_ohlines.dat > lowd_ir_ohlines.dat"
	print "# GDW 2001-Jan-08"
	print "#-----------------------------------------------------------------------"
}
/^#/{	print ; next }
{	wavelength = $1
	intensity = $2
	if( intensity < min_intensity){ next }
	if( intensity == last_intensity){
		new_wavelength = 0.5*(wavelength + last_wavelength)
		printf "%.3f %9.3e blend\n", new_wavelength, intensity
		last_wavelength = -1
		last_intensity = -1
	} else {
		if( last_intensity > 0){
			printf "%.3f %9.3e solo\n", last_wavelength, last_intensity
		}
		last_wavelength = wavelength
		last_intensity = intensity
	}
}

include <math.h>

#
# T_FLEX_FIT:  Simple tool for playing with fit to flexure data
# Created 29-Jul-97
# Modified 5-Aug-97 to give offsets

define N_ROT	12	# number of rotator positions
define N_EL	 7	# number of elevation positions


procedure t_flex_fit()

real 	x, y
real	rot, el
real	rrot, rel

real	delx1, delx2, dely1, dely2

real	clgetr()
begin
	x = clgetr ("x")
	y = clgetr ("y")
	rot = clgetr ("rot")
	el = clgetr ("el")
	rrot = clgetr ("rrot")
	rel = clgetr ("rel")

	call flex_corr (rot, el, delx1, dely1)
	call flex_corr (rrot, rel, delx2, dely2)

	x = x - (delx1 - delx2)
	y = y - (dely1 - dely2)

	call printf ("x=%7.1f  y=%7.1f \n")
		call pargr (x)
		call pargr (y)
end


#
# FLEX_CORR:  work out the flexure correction to the standand LRIS stow
# position (ROTPPOSN=90, EL=0).
# Currently defined values give 0.10 and 0.17 pixel rms in x and y.
#
# NOW IN MBOXFIND


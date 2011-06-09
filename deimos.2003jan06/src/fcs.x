# FCS: translate FCS coordinates to ICS

include <math.h>

define	PRESCAN 17.
define	DEL1	0.11
define	XCEN1	-4636.6
define	YCEN1	-645.5
define	DEL2	0.52
define	XCEN2	4595.3
define	YCEN2	-616.3

procedure t_fcs ()

real	xobs, yobs			# input mosaic pixel coordinates
real	xoff, yoff			# xoff, yoff are the offsets to mos cen


real	xccd, yccd			# Input CCD values
real	xics, yics			# returned ICS values
int	n
int	i

int	scan()
int	clgeti()
real	clgetr()

begin

# Read in parameters:
	xoff = clgetr ("xoff")
	yoff = clgetr ("yoff")
	n = clgeti ("fcs")

	while (scan() != EOF) {
		call gargr (xccd)
		call gargr (yccd)

		call fcs_to_ics (xccd, yccd, n, xics, yics)

		call printf ("%6.1f,%6.1f,%1d   %8.3f %8.3f    %8.3f %7.3f\n")
			call pargr (xccd)
			call pargr (yccd)
			call pargi (n)
			call pargr (xics)
			call pargr (yics)
			call pargr ((xics+xoff)*0.015)
			call pargr ((yics+yoff)*0.015)
	}

end

procedure	fcs_to_ics (xccd, yccd, n, xics, yics)

real	xccd, yccd
int	n
real	xics, yics

real	x, y
real	xoff, yoff
real	theta, delta
real	cost, sint

begin
	if (n == 1) {
		delta = DEL1
		xoff = XCEN1
		yoff = YCEN1
		theta = DEGTORAD (delta)
		cost = cos (theta)
		sint = sin (theta)

		x = (yccd + 300.) - 300.5
		y = -1. * ((xccd-PRESCAN) - 600.5)

		xics = x * cost - y * sint + xoff
		yics = x * sint + y * cost + yoff

	} else {
		delta = DEL2
		xoff = XCEN2
		yoff = YCEN2
		theta = DEGTORAD (delta)
		cost = cos (theta)
		sint = sin (theta)

		x = -1. * ((yccd + 300.) - 300.5)
		y = -1. * ((xccd-PRESCAN) - 600.5)

		xics = x * cost - y * sint + xoff
		yics = x * sint + y * cost + yoff

	}
end


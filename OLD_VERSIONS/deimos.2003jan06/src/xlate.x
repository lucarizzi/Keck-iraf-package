# XLATE: translate Mosaic pixel coordinates to ICS

include "instrument.h"
include "deimos.h"

procedure t_xlate ()

real	xobs, yobs			# input mosaic pixel coordinates
real	xoff, yoff			# xoff, yoff are the offsets of Chip 3

real	ccd[NCCD,3]				# CCD geometry
double	sys[NPARAM]				# system parameters

real	xccd, yccd			# Input CCD values
real	xics, yics			# returned ICS values
int	n
int	i

int	scan()
real	clgetr()

begin
# Initialize CCD geometry
	call ccd_geom (ccd, sys, YES)

# Read in parameters:
	xoff = clgetr ("xoff")
	yoff = clgetr ("yoff")

	n = 0
	while (scan() != EOF) {
		call gargr (xobs)
		call gargr (yobs)

		call get_ccdnum (xobs, yobs, n, xccd, yccd)
		call ccd_to_ics (xccd, yccd, ccd, n, xics, yics)

		call printf ("%7.3f %7.3f  %6.1f,%6.1f,%1d   %8.3f %8.3f    %8.3f %7.3f\n")
			call pargr (xobs)
			call pargr (yobs)
			call pargr (xccd)
			call pargr (yccd)
			call pargi (n)
			call pargr (xics)
			call pargr (yics)
			call pargr ((xics+xoff)*0.015)
			call pargr ((yics+yoff)*0.015)
	}

#	call phys_sil_corners (ccd, xoff, yoff)
end



procedure	phys_sil_corners (ccd, xoff, yoff)

real	ccd[NCCD,3]				# CCD geometry
real	xoff, yoff			# offsets due to mosaic decenter

int	n
int	i, j, k
real	xccd, yccd, xics, yics

real	xk[4], yk[4], xedge[4], yedge[4]

begin
	xk[1] =    1. ; yk[1] =    1. ; xedge[1] = -CCDXEDG
	xk[2] =    1. ; yk[2] = 4096. ; xedge[2] = -CCDXEDG
	xk[3] = 2048. ; yk[3] = 4096. ; xedge[3] =  CCDXEDG
	xk[4] = 2048. ; yk[4] =    1. ; xedge[4] =  CCDXEDG

	yedge[1] = -CCDYEDG
	yedge[2] =  CCDYEDG
	yedge[3] =  CCDYEDG
	yedge[4] = -CCDYEDG
	
	do j = 1, 2 {
	    do i = 1, 4 {

		n = (j-1)*4 + i
		call printf ("\n## CHIP %d:\n")
			call pargi (n)
		do k = 1, 4 {
			xccd = xk[k] - 1. 
			yccd = yk[k] # - 1. 

			call ccd_to_ics (xccd, yccd, ccd, n, xics, yics)

		call printf ("%6.1f,%6.1f,%1d   %8.3f %8.3f    %8.3f %7.3f\n")
				call pargr (xccd)
				call pargr (yccd)
				call pargi (n)
				call pargr (xics)
				call pargr (yics)
				call pargr (-((xics+xoff)*0.015+xedge[k]))
				call pargr (-((yics+yoff)*0.015+yedge[k]))
		}
	    }
	}
end

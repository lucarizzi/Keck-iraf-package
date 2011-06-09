# VERN: hacked version of XBOX, run only with INPUT

include	<imhdr.h>
include <gset.h>
include <math.h>
include	"instrument.h"
include "align.h"
include	"futil.h"

define	XAMPL	0.		# X-flexure amplitude (px)
define	XPHI	30.		# X-flexure phase offset
define	YAMPL	-7.		# Y-flexure amplitude (px)
define	YPHI	30.		# Y-flexure phase offset

define		CRIT_LEV 100.		# min average level w/in fwhm
define		SATURLEV  64000.	# Saturation level, approx. after bias
define		SZ_FITS	  40		# number of characters in keyword
define		SZ_LAB	40	# size of label in plot

define	DEF_PROP_FILE	"deimos$lib/prop/redprop.txt"	# default properties


procedure t_vern()

char	image[SZ_FNAME]			# Input Image 
char	input[SZ_FNAME]			# Input Box coordinates
char	logfile[SZ_FNAME]		# (opt) log file for moves
int	nxbox, nybox			# Size of subraster
real	xsz, ysz			# Full size of boxes (pix)
real	xfwhm, yfwhm			# FWHM of star (pix)
real	xoff, yoff			# x,y shifts to apply to input box coord
pointer	fda, fdl
int	find_coord			# find box/star coordinate pairs?

char	tchar
int	ncols, nlines
int	npt, ndx, i, j, k, n
int	stat

# new for mosaic images:
pointer	im0			# image pointer to primary HDU
pointer	mmap			# pointer to mosaic vectors

# For auto-grep
char	guiname[SZ_GUINAM]
char	label[SZ_LAB]
real	x0ics, y0ics		# ICS coords of DEIMOS pointing origin
pointer	bdat			# box data struct

int	nx, ny
real 	xs, ys
real	xflex, yflex
real	xb, yb, xstar, ystar
real	zmax, zmin
real	theta
int 	x1, y1, x2, y2
pointer	bufx, bufy, bufzx, bufzy, buftx, bufty
pointer	buf
pointer	xbuf1, ybuf1, xbuf2, ybuf2

real	rotposn					# PA = rotposn+90

int	zip_mode
int	nn, noff
real	pixval, grad

real	saw_xcorr()
int	box_graph(), box_center()
real	vsum1()

int	det_chip()
pointer	chip_sect()

bool	strne()
int	clgeti(), imaccf()
int	fscan(), nscan()
real	clgetr(), imgetr()
pointer	open()

begin
call eprintf ("Welcome to VERN for DEIMOS!\n")
	zip_mode = YES

	call clgstr ("image", image, SZ_FNAME)
	call clgstr ("logfile", logfile, SZ_FNAME)
	nxbox = clgeti ("nxwin")
	nybox = clgeti ("nywin")
	xsz = clgetr ("xsz")
	ysz = clgetr ("ysz")
	xfwhm = clgetr ("fwhm")
	yfwhm = clgetr ("fwhm")
	xoff  = clgetr ("xoff")
	yoff  = clgetr ("yoff")

#	ncols = IM_LEN(im,1)
#	nlines = IM_LEN(im,2)
ncols = 8192		# XXX
nlines = 4096		# XXX

	nx = nxbox		# fixed length now that chip_sect zero-fills
	ny = nybox

	find_coord = YES



# Open mosaic image:
	call mos_init (image, DEF_PROP_FILE, im0, mmap, n)

# Get relevant keywords ...
# ... for the PA and pointing origin
	if (imaccf (im0, "ROTATVAL") == YES) {
		rotposn = imgetr (im0, "ROTATVAL")
	} else {
		call fatal (0, "ROTATVAL missing!\n")
	}


	fdl = open (logfile, APPEND, TEXT_FILE)
	call seek (fdl, EOF)
	call fprintf (fdl, "%s %6.1f")
		call pargstr (image)
		call pargr (rotposn)

# READY: DO WE WANT TO FIND BOX/STAR COORDS ... ?
	if (find_coord == YES) {

## Find flexure correction: #  (XXX NOT NEEDED if FCS is on):
	    if (imaccf (im0, "ROTATVAL") == YES) {
		theta = imgetr (im0, "ROTATVAL")
		xflex = XAMPL * sin (DEGTORAD(theta+XPHI))
		yflex = YAMPL * sin (DEGTORAD(theta+YPHI))
call eprintf ("FLEXURE CORRECTION: ROTATVAL = %4f -->  X,Y flex = %4f  %4f\n")
call pargr (theta)
call pargr (xflex)
call pargr (yflex)
	    } else {
		call eprintf ("ROTATVAL not found -- no flexure applied\n")
	    }

############### FOR TEXT INPUT #######################
# Get and open the input file
	    call clgstr ("input", input, SZ_FNAME)
	    if (strne (input, "")) {
		fda = open (input, READ_ONLY, TEXT_FILE)

# Count entries in input file
		npt = 0
		while (fscan(fda) != EOF)
			npt = npt + 1
		call seek (fda, BOF)

# allocate the box struct
		call box_alloc (bdat, npt)

# Get the input entries
		ndx = 0
		while (fscan (fda) != EOF) {
			call gargwrd (tchar, 1)
			if (tchar == '#' || nscan() == 0) {
				next
			}
			call reset_scan()
			call gargr (xs)
			call gargr (ys)
			if (nscan() < 2) {
				call eprintf ("WARNING: input line skipped\n")
				next
			}

			call gargwrd (OBJNAME(bdat,ndx), SZ_ID-1)
			call gargc (PBAND(bdat,ndx))
			call gargr (MAG(bdat,ndx))

			if (nscan() < 3)
			    call strcpy ("(no ID)", OBJNAME(bdat,ndx), SZ_ID-1)
			if (nscan() < 4 || PBAND(bdat,ndx) == '#')
				PBAND(bdat,ndx) = 'm'
			if (nscan() < 5 || PBAND(bdat,ndx) == '#')
				MAG(bdat,ndx) = INDEF

			XPX(bdat,ndx) = xs
			YPX(bdat,ndx) = ys
			ndx = ndx + 1
		}
		npt = ndx	# npt reduced by bad pairs, comment lines

	    } else {

# get info from Mask Design tables appended to image
		call get_boxes (image, bdat, npt, guiname, x0ics, y0ics)

		call eprintf ("### This is mask  %s\n\n")
			call pargstr (guiname)

	    }

# Allocate arrays for marginal plots
	    call malloc (bufx, nxbox, TY_REAL)
	    call malloc (bufy, nybox, TY_REAL)
	    call malloc (bufzx, nxbox, TY_REAL)
	    call malloc (bufzy, nybox, TY_REAL)
	    call malloc (buftx, nxbox, TY_REAL)
	    call malloc (bufty, nybox, TY_REAL)

# Allocate the coordinate arrays (others deferred until needed)
	    call malloc (xbuf1, npt, TY_REAL)		# ref coord
	    call malloc (ybuf1, npt, TY_REAL)
	    call malloc (xbuf2, npt, TY_REAL)		# input coord
	    call malloc (ybuf2, npt, TY_REAL)



# Loop through expected box positions to measure actual box/star positions
	    ndx = 0
	    do k = 0, npt-1 {
# first, make sure XPX,YPX are OK:
		if (XPX(bdat,k) == INDEF || YPX(bdat,k) == INDEF) {
			call eprintf ("Box off chip -- skipping\n")
			next
		}

# adjust values for offsets:
		xs = XPX(bdat,k) + xflex + xoff
		ys = YPX(bdat,k) + yflex + yoff

# Find the Box; we do a kludge -- run this twice, recentered on box the 2nd time
		xb = xs
		yb = ys
		do j = 1, 2 {
		    x1 = xb - nxbox/2
		    x2 = x1 + nxbox - 1
		    y1 = yb - nybox/2
		    y2 = y1 + nybox - 1
if (xb == INDEF || yb == INDEF){
call eprintf ("xb/yb INDEF on pass %d!\n")
call pargi (j)
next
}
nx = x2 - x1 + 1	# XXX
ny = y2 - y1 + 1	# XXX

### NEED CHIP NUMBER
		    i = det_chip (xb, yb)
# call eprintf ("(%d) %d-%d  %d-%d \n")
# call pargi (i)
# call pargi (x1)
# call pargi (x2)
# call pargi (y1)
# call pargi (y2)

# Get the image section (limit checking and zero-fill included)
		    buf = chip_sect (x1, x2, y1, y2, mmap, i, YES)

# TMP? Check for common error of amp dropout
		    call alimr (Memr[buf], nx*ny, zmax, zmin)
		    if (zmax == zmin) {
			call eprintf ("Subraster contains no useful data!\n")
			stat = ERR
			break
		    }
		    if (zmax > SATURLEV)
			call eprintf ("WARNING! Star may be saturated!!\007\n")

# Fill position vectors
		    do i = 0, nx-1 {
			Memr[bufx+i] = i + x1
		    }
		    do i = 0, ny-1 {
			Memr[bufy+i] = i + y1
		    }
# Get the box position
		    stat = box_center (Memr[bufx], Memr[bufy], Memr[buf], nx,
			ny, Memr[bufzx], Memr[bufzy], xsz, ysz, xb, yb)

		    if (stat != OK)
			break
# On second time through loop recenter on box
		}

		if (stat == OK) {

		    call eprintf ("Box center:  %6.2f %6.2f  (x,yoff:%4.1f,%4.1f) (del:%5.1f,%5.1f)\n")
			call pargr (xb-xoff)
			call pargr (yb-yoff)
			call pargr (xoff)
			call pargr (yoff)
			call pargr (xb-xs)
			call pargr (yb-ys)

# Now get star location
# ... construct label
		    call sprintf (label, SZ_LAB, "%-s: %1s=%4.1f")
			call pargstr (OBJNAME(bdat,k))
			call pargc (PBAND(bdat,k))
			call pargr (MAG(bdat,k))
#		    stat = box_graph (Memr[bufx], Memr[bufy], Memr[bufzx],
#			Memr[bufzy], Memr[buftx], Memr[bufty], nx, ny, label,
#			xb, yb, xsz, ysz, xfwhm, yfwhm, xstar, ystar, zip_mode)

		    noff = xsz / 2. + 0.5
	            xstar = saw_xcorr (Memr[bufx], Memr[bufzx], nx, xfwhm, xfwhm/2., noff, grad)
		    noff = ysz / 2. + 0.5
	            ystar = saw_xcorr (Memr[bufy], Memr[bufzy], ny, yfwhm, yfwhm/2., noff, grad)
#call eprintf ("DEBUG_B %f %f\n"); call pargr (xstar); call pargr (ystar)

# Store box/star coordinates
			x1 = xstar + 0.5 - 0.5 * xfwhm
			x2 = xstar + 0.5 + 0.5 * xfwhm
			y1 = ystar + 0.5 - 0.5 * yfwhm
			y2 = ystar + 0.5 + 0.5 * yfwhm
			nn = (x2-x1+1) * (y2-y1+1)
		    	i = det_chip (xstar, ystar)
		    	buf = chip_sect (x1, x2, y1, y2, mmap, i, YES)
			pixval = vsum1 (Memr[buf], nn) / nn
			if (pixval < CRIT_LEV) {
				xstar = -99.
				ystar = -99.
				xb = -99.
				yb = -99.
		    }
		} else {
			xstar = -99.
                        ystar = -99.
                        xb = -99.
                        yb = -99.
		}
call fprintf (fdl, "%7.1f %6.1f %5.0f ")
call pargr (xstar)
call pargr (ystar)
call pargr (pixval)
	    }

# Done with the box-finding: clean up
	    call mfree (bufty, TY_REAL)
	    call mfree (buftx, TY_REAL)
	    call mfree (bufzy, TY_REAL)
	    call mfree (bufzx, TY_REAL)
	    call mfree (bufy, TY_REAL)
	    call mfree (bufx, TY_REAL)
	    call close (fda)

	}

# Close the mosaic image
	call mos_free (im0, mmap)

call fprintf (fdl, "\n")
		call close (fdl)
end


# INSTALLATION/REVISION NOTES:
# Must change:
#	-- pixel scale [DONE]
#	-- PO Names and pixel locations [Done, except for DREF]
#	-- rotator offsets??  [ROFF -- done; should set up in instrum coord]
#	-- any detector related values
#	-- DCS host name, system name [DONE]

#TBD:	-- review above.
#	-- replace Y-weighting in xb_util.x
#	-- should really use ICS coords, not mosaic coords

include	<imhdr.h>
include <gset.h>
include <math.h>
include	"futil.h"

define	ROFF	90.		# rotational offset ROTPOSN--TO__PA
define DCS_HOST	"polo"		# name of DCS machine
define DCS_SYS	"dcs2"		# name of DCS system

define		ID_CHSZ	9		# Character size of ID string
# define	REL_HT	0.4		# relative height for crossings
# define	LAB_SZ 80		# Char size of title/label strings
# define	X_HWID	1.5		# x half-width of triangle function
# define	Y_HWID	1.5		# y half-width of triangle function
# define	X_CRAD	13		# x centering radius
# define	Y_CRAD	13		# y centering radius
# define	FRAC	0.333		# fractional lev. in sorted list for sky

define		ASECPIX  0.1192		# Arcsec/pixel
define		PA_ACCEPT 0.010		# this is 0.1" at ~4000 px
define		SZ_FITS	  40		# number of characters in keyword

define	DEF_PROP_FILE	"deimos$lib/prop/redprop.txt"	# default properties

#
# XBOX: locate boxes and stars within them, and solve for rotation
# and translation.
# This program is a combination of two previous tasks, mboxfind and maskalign.
#
# Since the boxes are relatively small, read in section and deal with the 
# 2-D image rather than cuts
#
# We currently have a kludge for robustness -- runs box-center twice, the
# second time after recentering image section on initial box center

procedure t_xbox()

char	image[SZ_FNAME]			# Input Image 
char	input[SZ_FNAME]			# Input Box coordinates
char	coordfile[SZ_FNAME]		# Box-Star coordinates (pre-measured)
int	nxbox, nybox			# Size of subraster
real	xsz, ysz			# Full size of boxes (pix)
real	xfwhm, yfwhm			# FWHM of star (pix)
real	xoff, yoff			# x,y shifts to apply to input box coord
pointer	fda, fdb
bool	find_coord			# find box/star coordinate pairs?

char	tchar, idstr[ID_CHSZ]
int	ncols, nlines
int	npt, ndx, i, j
bool	stat

# new for mosaic images:
pointer	im0			# image pointer to primary HDU
pointer	mmap			# pointer to mosaic vectors
int	nccd

int	nx, ny
real 	xs, ys
real	xb, yb, xstar, ystar
int 	x1, y1, x2, y2
pointer	bufx, bufy, bufzx, bufzy, buftx, bufty
pointer	buf

# parameters for aligning
real	xrot, yrot
bool	invert
# int	nx, ny
int	niter
real	bxs, def1err, def2err
bool	dcs, dcspa
bool	no_op

char	poname[SZ_FITS]				# Pointing Origin name
char	cmdline[SZ_LINE]			# Command line for rsh
real	rotposn					# PA = rotposn+90
real	coeff[2,3], err[3]
real	sina, cosa
real	earcsec, narcsec, eerr, nerr, pdeg, perr
pointer	xbuf1, ybuf1, xbuf2, ybuf2
pointer	ebufx1, ebufy1, ebufx2, ebufy2
pointer	rbufx, rbufy, wbuf

int	zip_mode

bool	box_graph()
int	det_chip()
pointer	chip_sect()

bool	clgetb(), imaccf(), streq(), strne()
int	clgeti()
int	fscan(), nscan()
int	oscmd()
real	clgetr(), imgetr()
pointer	open()

begin
call eprintf ("Welcome to XBOX for DEIMOS!\n")
	zip_mode = NO

	call clgstr ("image", image, SZ_FNAME)
	nxbox = clgeti ("nxwin")
	nybox = clgeti ("nywin")
	xsz = clgetr ("xsz") / ASECPIX
	ysz = clgetr ("ysz") / ASECPIX
	xfwhm = clgetr ("fwhm") / ASECPIX
	yfwhm = clgetr ("fwhm") / ASECPIX
	xoff  = clgetr ("xoff")
	yoff  = clgetr ("yoff")
	no_op = clgetb ("practice")

#	ncols = IM_LEN(im,1)
#	nlines = IM_LEN(im,2)
ncols = 8192		# XXX
nlines = 4096		# XXX

	nx = nxbox		# fixed length now that chip_sect zero-fills
	ny = nybox

	call clgstr ("pairs", coordfile, SZ_FNAME)
	if (streq (coordfile, "")) {
		find_coord = YES
	} else {
		fdb = open (coordfile, READ_ONLY, TEXT_FILE)
		find_coord = NO
	}

# Open mosaic image:
	call mos_init (image, DEF_PROP_FILE, im0, mmap, nccd)

# Get relevant keywords ...
# ... for the PA and pointing origin
	if (imaccf (im0, "ROTPOSN")) {
		rotposn = imgetr (im0, "ROTPOSN")
	} else {
		call eprintf ("ROTPOSN missing; -180 (LRIS stow) assumed\n")
		rotposn = -180.
	}

	if (imaccf (im0, "PONAME")) {
		call imgstr (im0, "PONAME", poname, SZ_FITS)
	} else {
		call strcpy ("UNKNOWN", poname, SZ_FITS)
	}
	
# Assign the proper rotation center for the PO
	if (streq (poname, "REF")) {
		xrot = 4884.9
		yrot = 3792.5
	} else if (streq (poname, "DEIMOS")) {
		xrot = 4096.
		yrot = 3149
	} else if (streq (poname, "Slit")) {
		xrot = 5132.
		yrot = 3064
	} else if (streq (poname, "Image")) {
		xrot = 5132.
		yrot = 3148
	} else {
		if (strne (poname, ""))
			call eprintf ("Unknown Pointing Origin!\n")
		xrot = 4096.5
		yrot = 3150
	}

# READY: DO WE WANT TO FIND BOX/STAR COORDS?
	if (find_coord == YES) {

# Get and open the input file
	    call clgstr ("input", input, SZ_FNAME)
	    fda = open (input, READ_ONLY, TEXT_FILE)

# Allocate arrays for marginal plots
	    call malloc (bufx, nxbox, TY_REAL)
	    call malloc (bufy, nybox, TY_REAL)
	    call malloc (bufzx, nxbox, TY_REAL)
	    call malloc (bufzy, nybox, TY_REAL)
	    call malloc (buftx, nxbox, TY_REAL)
	    call malloc (bufty, nybox, TY_REAL)

# Count entries in input file
	    npt = 0
	    while (fscan(fda) != EOF)
		npt = npt + 1
	    call seek (fda, BOF)

# Allocate the coordinate arrays (others deferred until needed)
	    call malloc (xbuf1, npt, TY_REAL)		# ref coord
	    call malloc (ybuf1, npt, TY_REAL)
	    call malloc (xbuf2, npt, TY_REAL)		# input coord
	    call malloc (ybuf2, npt, TY_REAL)


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

# adjust values for offsets:
		xs = xs + xoff
		ys = ys + yoff

		call gargwrd (idstr, ID_CHSZ)
		if (nscan() < 3)
			call strcpy ("(no ID)", idstr, ID_CHSZ)

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
call eprintf ("(%d) %d-%d  %d-%d \n")
call pargi (i)
call pargi (x1)
call pargi (x2)
call pargi (y1)
call pargi (y2)

# Get the image section (limit checking and zero-fill included)
		    buf = chip_sect (x1, x2, y1, y2, mmap, i, YES)

# Fill position vectors
		    do i = 0, nx-1 {
			Memr[bufx+i] = i + x1
		    }
		    do i = 0, ny-1 {
			Memr[bufy+i] = i + y1
		    }
# Get the box position
		    call box_center (Memr[bufx], Memr[bufy], Memr[buf], nx, ny,
			Memr[bufzx], Memr[bufzy], xsz, ysz, xb, yb)

# On second time through loop recenter on box
		}

		call printf ("Box center:  %6.2f %6.2f  (%4.1fx,%4.1fy removed) (del:%4.1f,%4.1f)\n")
			call pargr (xb-xoff)
			call pargr (yb-yoff)
			call pargr (xoff)
			call pargr (yoff)
			call pargr (xb-xs)
			call pargr (yb-ys)

# Now get star location
		stat = box_graph (Memr[bufx], Memr[bufy], Memr[bufzx],
			Memr[bufzy], Memr[buftx], Memr[bufty], nx, ny, idstr,
			xb, yb, xsz, ysz, xfwhm, yfwhm, xstar, ystar, zip_mode)

# Store box/star coordinates
		if (stat == OK) {
			Memr[xbuf1+ndx] = xstar
			Memr[ybuf1+ndx] = ystar
			Memr[xbuf2+ndx] = xb
			Memr[ybuf2+ndx] = yb
			ndx = ndx + 1
####### THIS SECTION FOR RAJA:
# i = det_chip (xstar, ystar)
# call printf ("RAJA: %1d  %7.2f %7.2f %-s\n")
# call pargi (i)
# call pargr (xstar - 2048.*(i-1))
# call pargr (ystar)
# call pargstr (idstr)
		}
	    }
	    npt = ndx		# npt reduced by bad pairs, comment lines

# Done with the box-finding: clean up
	    call mfree (bufty, TY_REAL)
	    call mfree (buftx, TY_REAL)
	    call mfree (bufzy, TY_REAL)
	    call mfree (bufzx, TY_REAL)
	    call mfree (bufy, TY_REAL)
	    call mfree (bufx, TY_REAL)
	    call close (fda)

	} else {
# ... OR READ IN COORDINATE PAIRS?

# Count entries in input file
	    npt = 0
	    while (fscan(fdb) != EOF)
		npt = npt + 1
	    call seek (fdb, BOF)

# Allocate the coordinate arrays (others deferred until needed)
	    call malloc (xbuf1, npt, TY_REAL)		# ref coord
	    call malloc (ybuf1, npt, TY_REAL)
	    call malloc (xbuf2, npt, TY_REAL)		# input coord
	    call malloc (ybuf2, npt, TY_REAL)

# Get the input entries
	    ndx = 0
	    while (fscan (fdb) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#' || nscan() == 0) {
			next
		}
		call reset_scan()
		call gargr (xb)
		call gargr (yb)
		call gargr (xstar)
		call gargr (ystar)
		if (nscan() < 4) {
			call eprintf ("WARNING: input line skipped\n")
			next
		}
		Memr[xbuf1+ndx] = xstar
		Memr[ybuf1+ndx] = ystar
		Memr[xbuf2+ndx] = xb
		Memr[ybuf2+ndx] = yb
		ndx = ndx + 1
	    }
	    npt = ndx		# npt reduced by bad pairs, comment lines
	    call close (fdb)
	}

# Close the mosaic image
	call mos_free (im0, mmap)

# ON TO ALIGNING: setup
	nx = ncols			# XXX TMP: should fix; also symmetrize
	ny = nlines			# XXX TMP   ... maskalign
	bxs = -2. * clgetr ("box_size")
	niter = clgeti ("niter")
	def1err = clgetr ("def_ref_err")
	def2err = clgetr ("def_err")

# Allocate the remaining arrays needed
	call malloc (ebufx1, npt, TY_REAL)		# error arrays
	call malloc (ebufy1, npt, TY_REAL)
	call malloc (ebufx2, npt, TY_REAL)
	call malloc (ebufy2, npt, TY_REAL)
	call malloc (rbufx, npt, TY_REAL)		# residuals
	call malloc (rbufy, npt, TY_REAL)
	call malloc (wbuf, npt, TY_REAL)		# weights

	call amovkr (1., Memr[wbuf], npt)

# Store default error estimates (someday could be individual...)
	call amovkr (def1err, Memr[ebufx1], npt)
	call amovkr (def1err, Memr[ebufy1], npt)
	call amovkr (def2err, Memr[ebufx2], npt)
	call amovkr (def2err, Memr[ebufy2], npt)

# square the errors, as required for get_lsqf:
	call amulr (Memr[ebufx1], Memr[ebufx1], Memr[ebufx1], npt)
	call amulr (Memr[ebufy1], Memr[ebufy1], Memr[ebufy1], npt)
	call amulr (Memr[ebufx2], Memr[ebufx2], Memr[ebufx2], npt)
	call amulr (Memr[ebufy2], Memr[ebufy2], Memr[ebufy2], npt)

# Reference coordinates to center of rotation axis
	call aaddkr (Memr[xbuf1], -xrot, Memr[xbuf1], npt)
	call aaddkr (Memr[ybuf1], -yrot, Memr[ybuf1], npt)
	call aaddkr (Memr[xbuf2], -xrot, Memr[xbuf2], npt)
	call aaddkr (Memr[ybuf2], -yrot, Memr[ybuf2], npt)


	call maskalign (Memr[xbuf1], Memr[ybuf1], Memr[xbuf2], Memr[ybuf2],
		Memr[ebufx1], Memr[ebufy1], Memr[ebufx2], Memr[ebufy2],
		Memr[rbufx], Memr[rbufy], Memr[wbuf], npt, nx, ny, niter,
						xrot, yrot, bxs, coeff, err)
# Print out results
	call printf ("#   obj-x   obj-y   slit-x   xres   slit-y   yres     w \n")
	do i = 0, npt-1 {
		call printf ("#%8.2f%8.2f  %7.2f (%5.3f) %7.2f (%5.3f) %6.2f\n")
			call pargr(Memr[xbuf1+i])
			call pargr(Memr[ybuf1+i])
			call pargr(Memr[xbuf2+i])
			call pargr(Memr[rbufx+i])
			call pargr(Memr[ybuf2+i])
			call pargr(Memr[rbufy+i])
			call pargr(Memr[wbuf+i])
	}

	call printf ("\n\n#  x-xform: %7.5fx + %7.5fy  + %7.3f (%5.3f)\n")
		call pargr (coeff[1,1])
		call pargr (coeff[1,2])
		call pargr (coeff[1,3])
		call pargr (err[2])
	call printf ("#  y-xform: %7.5fx + %7.5fy  + %7.3f (%5.3f)\n")
		call pargr (coeff[2,1])
		call pargr (coeff[2,2])
		call pargr (coeff[2,3])
		call pargr (err[3])

# Now adjust for existing PA (convert x,y into E,N):
# Note that "angle" = (angle of "PA" wrt x-axis) - (PA = rotposn+90)
#        ...        = 90 - (rotposn+90)
	invert = NO			# Left for historical reasons
	if (invert) {
		call eprintf ("Sorry -- not yet tested\n")
		pdeg = RADTODEG( atan2(coeff[1,2], coeff[1,1]))
		perr = RADTODEG(asin(err[1]))
		cosa = cos (DEGTORAD(-rotposn-ROFF))
		sina = sin (DEGTORAD(-rotposn-ROFF))
		narcsec = ASECPIX * (-cosa * coeff[1,3] + sina * coeff[2,3])
		earcsec = ASECPIX * ( sina * coeff[1,3] + cosa * coeff[2,3])
		nerr = ASECPIX * sqrt ((cosa * err[2])**2 + (sina * err[3])**2)
		eerr = ASECPIX * sqrt ((sina * err[2])**2 + (cosa * err[3])**2)
	} else {
		pdeg = RADTODEG(-atan2(coeff[1,2], coeff[1,1]))
		perr = RADTODEG(asin(err[1]))
		cosa = cos (DEGTORAD(-rotposn-ROFF))
		sina = sin (DEGTORAD(-rotposn-ROFF))
		narcsec = ASECPIX * ( cosa * coeff[1,3] + sina * coeff[2,3])
		earcsec = ASECPIX * (-sina * coeff[1,3] + cosa * coeff[2,3])
		nerr = ASECPIX * sqrt ((cosa * err[2])**2 + (sina * err[3])**2)
		eerr = ASECPIX * sqrt ((sina * err[2])**2 + (cosa * err[3])**2)
	}
	
# PA offset is backwards from rotation:
	call printf ("\n====================================================\n")
	call printf ("\n*** MOVE TELESCOPE/ROTATOR by the following offsets:\n")
	call printf ("\n Offset PA by %6.3f (%5.3f) degree\n")
		call pargr (pdeg)
		call pargr (perr)
	call printf (" Offsets: %6.2f\" EAST (%4.2f)\t %6.2f\" NORTH (%4.2f)\n\n")
		call pargr (earcsec)
		call pargr (eerr)
		call pargr (narcsec)
		call pargr (nerr)
	call printf ("====================================================\n\n")
	call flush (STDOUT)

## Is this just a practice run?
	if (no_op)
		call fatal (0, "ignore error -- terminated because PRACTICE set")
	
##  Apply offsets?
	dcs = clgetb ("dcs")

	if (dcs) {
# Check rotator error for magnitude
	    if (abs (pdeg) < PA_ACCEPT) {
		call printf ("PA error is small...  ")
		dcspa = clgetb ("dcs_rot")
	    } else {
		dcspa = YES
	    }
# Apply rotation
	    if (dcspa) {
		call sprintf (cmdline, SZ_LINE,
		  "rsh %s modify -s %s ROTDEST=%.3f ROTMODE=1")	# (1=pos.angle)
			call pargstr (DCS_HOST)
			call pargstr (DCS_SYS)
			call pargr (rotposn+pdeg)
		call eprintf ("\n(sending)  %s \n")
			call pargstr (cmdline)
		stat = oscmd (cmdline)
		if (stat != OK) {
			call eprintf ("command failed!  (%d)\n")
				call pargi (stat)
		}

		call sprintf (cmdline, SZ_LINE,
			"rsh %s waitfor -s %s ROTSTAT=8")	# (8=tracking)
			call pargstr (DCS_HOST)
			call pargstr (DCS_SYS)
		call eprintf ("\n(sending)  %s \n")
			call pargstr (cmdline)
		stat = oscmd (cmdline)
		if (stat != OK) {
			call eprintf ("command failed!  (%d)\n")
				call pargi (stat)
		}
	    }

# ... and apply translation
		call sprintf (cmdline, SZ_LINE,
		  "rsh %s modify -s %s RAOFF=%.2f DECOFF=%.2f REL2CURR=1")
			call pargstr (DCS_HOST)
			call pargstr (DCS_SYS)
			call pargr (earcsec)
			call pargr (narcsec)
		call eprintf ("\n(sending)  %s \n")
			call pargstr (cmdline)
		stat = oscmd (cmdline)
		if (stat != OK) {
			call eprintf ("command failed!  (%d)\n")
				call pargi (stat)
		}

		call sprintf (cmdline, SZ_LINE,
			"rsh %s waitfor -s %s AXESTAT=64")	# (64=tracking)
			call pargstr (DCS_HOST)
			call pargstr (DCS_SYS)
		call eprintf ("\n(sending)  %s \n")
			call pargstr (cmdline)
		stat = oscmd (cmdline)
		if (stat != OK) {
			call eprintf ("command failed!  (%d)\n")
				call pargi (stat)
		}

		call eprintf ("... done! \007 \n")
	}

end
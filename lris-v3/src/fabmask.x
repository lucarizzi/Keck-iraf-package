include	<math.h>
include	<error.h>
include <time.h>

define		MM_INCH		25.4	# mm per inch
define		PUNCH_HT	3.0	# Keck punch height (mm)
define		MAX_PUNCH	1000	# Limit to punches (1000 for AUTOSLIT)
define		Y_STEP		2.5	# Keck punch y-step size (mm)
define		X_STEPFCT	1.5  	# Keck punch x-step factor
define		NPATTERN	20	# max number of x-steps for boxes (>3)
define		MM_ASEC		0.7253	# mm/arcsec, for slit width conversion

define		SZ_FMT		8	# Char length of file format string
define		SZ_LAB		25	# Char length of mask name, owner

define		CUT_SLIT	0	# Code to mill/punch a slit
define		CUT_BOX		1	# Code to mill/punch a box
define		CODE_GS		-1	# Priority code for guide star   XXX
define		CODE_AS		-2	# Priority code for alignment star   XXX

define		PUNCH_XOFF	0.0	# X-offset for Keck Punch
define		PUNCH_YOFF	0.0	# Y-offset for Keck Punch
define		MILL_XOFF     -240.285	# X-offset for Mill (from phys. edge)
define		MILL_YOFF	10.0	# Y-offset for Mill (from phys. edge)

#
# FABMASK:  Translate autoslit "setup" file into  coords that can be used by
# the Lick Shops or Keck punch machine to fabricate a slit mask.
#
# REVISIONS:  for Lick shops, now creates "ACSII" format files, which are
# then "digitized" in "line" mode
# Rev: 27-jun-97 --- restored first corner as final corner. Changed all
# "keck" references to "punch"; all "lick" to "mill"
# Rev:  2-jul-97 --- removed slitsep, put in mapmask
# Rev: 16-jul-97 --- changed sign on Y-axis to conform to Keck preferences
# Rev: 18-jul-97 --- changed origin of Y-axis to conform to Keck preferences

procedure t_fabmask()

char	input[SZ_FNAME]
char	output[SZ_FNAME]
char	outfmt[SZ_FMT]
real	xzpt, yzpt
real	slitwid, punchwid
pointer	fda, fdb

char	tchar
int	npts, ndx
int	priority
real	arg1, arg2, arg3, arg4, arg5, arg6, arg7
pointer	bufx, bufy, bufy1, bufy2, buftype, bufpa

bool	streq()
int	fscan(), nscan()
real	clgetr()
pointer	open()

begin
	call clgstr ("input", input, SZ_FNAME)
	fda = open (input, READ_ONLY, TEXT_FILE)
	call clgstr ("output", output, SZ_FNAME)
	fdb = open (output, NEW_FILE, TEXT_FILE)
	call clgstr ("format", outfmt, SZ_FMT)

	xzpt = clgetr ("x_zpt")
	yzpt = clgetr ("y_zpt")
	if (streq (outfmt, "punch")) {
		xzpt = xzpt + PUNCH_XOFF
		yzpt = yzpt + PUNCH_YOFF
		punchwid = clgetr ("punch_width")	
	} else {
		xzpt = xzpt + MILL_XOFF
		yzpt = yzpt + MILL_YOFF
		slitwid = clgetr ("slit_width") * MM_ASEC
	}

# Count the entries
	ndx = 0
	while (fscan (fda) != EOF) {
		call gargwrd (tchar, 1)
		call reset_scan()
		if (tchar == '#') {
			next
		}
		ndx = ndx + 1
	}
	call seek (fda, BOF)

# Allocate arrays
	call malloc (bufx, ndx, TY_REAL)
	call malloc (bufy, ndx, TY_REAL)
	call malloc (bufy1, ndx, TY_REAL)
	call malloc (bufy2, ndx, TY_REAL)
	call malloc (bufpa, ndx, TY_REAL)
	call calloc (buftype, ndx, TY_INT)

# Get the input entries
	ndx = 0
	while (fscan (fda) != EOF) {
		call gargwrd (tchar, 1)
		call reset_scan()
		if (tchar == '#') {
			next
		}
		call gargr (arg1)	# xobj
		call gargr (arg2)	# yobj
		call gargr (arg3)	# ymin
		call gargr (arg4)	# ymax
		call gargr (arg5)	# rel. PA
		call gargr (arg6)	# ccdx
		call gargr (arg7)	# ccdy
		call gargi (priority)
		if (nscan() < 8) {
			call eprintf ("WARNING: input line skipped\n")
			next
		}
		Memr[bufx+ndx]   = arg1 + xzpt
		Memr[bufy+ndx]   = arg2 + yzpt
		Memr[bufy1+ndx]  = arg4 + yzpt
		Memr[bufy2+ndx]  = arg3 + yzpt

#	Work out PA
		if (arg5 == INDEF)
			arg5 = 0.
		if (arg5 != 0.) {
			if (arg5 < -90.) {
				arg5 = arg5 + 180.
			} else if (arg5 > 90.) {
				arg5 = arg5 - 180.
			}
		}
		Memr[bufpa+ndx] = arg5

#	Determine type
		if (priority >= 1 && priority <= 9999) {
			Memi[buftype+ndx] = CUT_SLIT
		} else if (priority == CODE_GS) {
			Memi[buftype+ndx] = CUT_SLIT
			Memr[bufy1+ndx]  = Memr[bufy+ndx] - 0.5
			Memr[bufy2+ndx]  = Memr[bufy+ndx] + 0.5
		} else if (priority == CODE_AS) {
			Memi[buftype+ndx] = CUT_BOX
		} else {
			call eprintf ("Unused priority code -- skipped!\n")
			next
		}
			
		ndx = ndx + 1
	}
	npts = ndx
	call close (fda)

# Write out the slit positions:
	if (streq (outfmt, "punch")) {
		call punchmask (fdb, Memr[bufx], Memr[bufy1], Memr[bufy2],
				Memi[buftype], npts, punchwid)
	} else if (streq (outfmt, "mill")) {
		call millmask (fdb, Memr[bufx], Memr[bufy], Memr[bufy1],
			Memr[bufy2], Memr[bufpa], Memi[buftype], npts, slitwid)
	} else {
		call eprintf ("Unrecognized mask format!\n")
	}

# Close up
	call close (fdb)

	call mfree (buftype, TY_INT)
	call mfree (bufpa, TY_REAL)
	call mfree (bufy2, TY_REAL)
	call mfree (bufy1, TY_REAL)
	call mfree (bufx, TY_REAL)
end

#
# PUNCHMASK:  makes a Keck PUNCH file
#
# Modified 25 sep to remove double-punch option, but to allow a scheme for
# multi-punching boxes
# Modified 27-jun-97: renamed keck --> punch to avoid confusion with keck mill

procedure punchmask (fdb, x, ymin, ymax, type, npts, punchwid)

pointer	fdb				# Output file
real	x[npts], ymin[npts], ymax[npts]	# Coords (mm)
int	type[npts]			# Type 1=slit; 4=box
int	npts
real	punchwid			# punch width (mm)

int	ndx
real	xsz, ysz			# align box size (mm)
real	xpat[NPATTERN]			# x-punch pattern for box
int	npat				# number of steps in box pattern
pointer bufx, bufy			# punch x,y vectors

int	k
int	i, j, j2
real	xpunch, ypunch
real	xstep
real	ystart, yend, yoff 
real	py				# number of punches in y
real	pxb, pyb			# number of punches in x,y for box
int	np
real	xoffbox, yoffbox		# offsets to first punch for box

begin
# Setup for slits: (none)
	
# Setup for boxes:
	ysz = PUNCH_HT			# TMP FIX -- force size to PUNCH_HT
	xsz = ysz
	yoffbox = max (0.5 * (ysz - PUNCH_HT), 0.)
	pyb = (ysz - PUNCH_HT) / Y_STEP + 1.
	if (pyb < 1.0) {
		call eprintf ("Box shorter than punch height\n")
		pyb = 1.0
	}
# Set up box pattern here: we don't want punches with uneven pressure
	xoffbox = max (0.5 * (xsz - punchwid), 0.)
#  ... we want the edges to be clean -- start there
	xpat[1] =  xoffbox
	xpat[2] = -xoffbox
	npat = 2
	pxb = xsz / punchwid - 2	# remaining width, punch units
	if (pxb > 1.0) {
		np = (pxb-1.) / 1.6 + 1		# number of punches in space
		xstep = (pxb + 1.) / (np + 1.)  # step factor between punches
		xpunch = -xoffbox
		do i = 1, np {
			npat = npat + 1
			if (npat > NPATTERN)
				call fatal (0, "Increase NPATTERN\n")
			xpunch = xpunch + xstep * punchwid
			xpat[npat] = xpunch
		}
		xpunch = 0.5 * (-xoffbox + xpat[3]) - xstep * punchwid
		do i = 1, np+1 {
			npat = npat + 1
			if (npat > NPATTERN)
				call fatal (0, "Increase NPATTERN\n")
			xpunch = xpunch + xstep * punchwid
			xpat[npat] = xpunch
		}
	} else {
		npat = 3
		xpat[npat] = 0.
	}

	call malloc (bufx, MAX_PUNCH, TY_REAL)
	call malloc (bufy, MAX_PUNCH, TY_REAL)

# Work out the punch positions:
	ndx = 0
	yoff = 0.5 * PUNCH_HT
	do k = 1, npts {
# SLITS:
	    if (type[k] == CUT_SLIT) {
		py = (ymax[k] - ymin[k] - PUNCH_HT) / Y_STEP + 1.
		ystart = ymax[k] - yoff
		yend   = ymin[k] + yoff
		if (py < 1.0) {
			call eprintf ("Slit %d shorter than punch height\n")
				call pargi (k)
			py = 1.0
			ystart = 0.5 * (ymax[k] + ymin[k])
			yend = ystart
		}

# punch the slit (single punches only)
		j2 = py
		if ((py - j2) > 0.)
			j2 = j2 + 1

		ypunch = ystart
		do j = 1, j2 {
			Memr[bufx+ndx] = x[k]
			Memr[bufy+ndx] = max (ypunch, yend)
			ndx = ndx + 1
			ypunch = ypunch - Y_STEP
		}
# BOXES:
	    } else if (type[k] == CUT_BOX) {
		py = pyb
		ystart = 0.5 * (ymax[k] + ymin[k]) + yoffbox
		yend   = 0.5 * (ymax[k] + ymin[k]) - yoffbox

# punch the box (complex pattern)
		j2 = pyb
		if ((py - j2) > 0.)
			j2 = j2 + 1

		ypunch = ystart
		do j = 1, j2 {
		    ypunch = max (ypunch, yend)
		    do i = 1, npat {
			Memr[bufx+ndx] = x[k] + xpat[i]
			Memr[bufy+ndx] = ypunch
			ndx = ndx + 1
		    }
		    ypunch = ypunch - Y_STEP
		}
	    }
	}

# Sort the final list in Y
	call r2sortd (Memr[bufy], ndx, Memr[bufx])

# HEADER:
# there is some initialization that goes here...
	call fprintf (fdb, "-9999,-9999\n")

	do i = 0, ndx-1 {
		call fprintf (fdb, "  %9.5f, %9.5f\n")
			call pargr (Memr[bufx+i])
			call pargr (Memr[bufy+i])
	}

	call mfree (bufy, TY_REAL)
	call mfree (bufx, TY_REAL)
end

#
# MILLMASK:  makes a Lick shop MILL file, in ASCII format (now with tilts)
# Rev:  8 oct 96 -- specify corners
# Mod: 27-jun-97: renamed lick --> mill
#

procedure millmask (fdb, x, y, ymin, ymax, theta, type, npts, slitwid)

pointer	fdb						# Output file
real	x[npts], y[npts], ymin[npts], ymax[npts]	# Coords (mm)
real	theta[npts]					# tilt-angle of slits
int	type[npts]					# Type 1=slit; 4=box
int	npts
real	slitwid						# slit width in mm

int	i
real	boxsz						# align box size (mm)
real	x1, x2, y1, y2					# ends of slit
real	x1a, x1b, x2a, x2b				# x corners of slit/box
real	maxtheta, mintheta				# extremes of theta

char	idmask[SZ_LAB], owner[SZ_LAB]			# mask ID and owner
char	date[SZ_TIME]					# time/date of file
long	time

long	clktime()

begin
# Write comments at beginning of file
	call clgstr ("idmask", idmask, SZ_LAB)
	call clgstr ("owner", owner, SZ_LAB)

	call fprintf (fdb, "# (     MaskName:    %-25s )\n")
		call pargstr (idmask, SZ_LAB)

	call fprintf (fdb, "# (     MaskOwner:   %-25s )\n")
		call pargstr (owner, SZ_LAB)

# Stamp with time
	time = clktime (long (0))
	call cnvtime (time, date, SZ_TIME)
	call fprintf (fdb, "# (     Date:        %-25s )\n")
		call pargstr (date, SZ_TIME)

	call fprintf (fdb, "# (     Slit width:  %5.3f mm  = %5.2f arcsec  )\n")
		call pargr (slitwid)
		call pargr (slitwid / MM_ASEC)

# work out max tool size:
	call alimr (theta, npts, mintheta, maxtheta)
	maxtheta = DEGTORAD (max (abs (mintheta), abs (maxtheta)))

	call fprintf (fdb, "# (     Max Tool Sz: %5.3f mm  = %5.1f thou    )\n")
		call pargr (slitwid * cos (maxtheta))
		call pargr (slitwid * cos (maxtheta) / 25.4 * 1000.)

# Stamp with internal units
	call fprintf (fdb, "# (     Units are:   MM                        )\n")
	call fprintf (fdb, "# (                                            )\n")

# Now loop through entities
	do i = 1, npts {
		y1 = ymin[i]
		y2 = ymax[i]
		if (type[i] == CUT_SLIT) {
			x1 = x[i] + (y1 - y[i]) * tan (DEGTORAD(theta[i]))
			x2 = x[i] + (y2 - y[i]) * tan (DEGTORAD(theta[i]))
			x1a = x1 - slitwid/2.
			x1b = x1 + slitwid/2.
			x2a = x2 - slitwid/2.
			x2b = x2 + slitwid/2.
		} else if (type[i] == CUT_BOX) {
			boxsz = y2 - y1
			x1a = x[i] - boxsz/2.
			x1b = x[i] + boxsz/2.
			x2a = x[i] - boxsz/2.
			x2b = x[i] + boxsz/2.
			y1 = y[i] - boxsz/2.
			y2 = y[i] + boxsz/2.
		}
# Note: the order x,y are reversed, as the milling stage has a RHS long x-axis
# and short y-axis. The masks are then cut from the front.
		call fprintf (fdb, "%10.4f %10.4f 0.00000\n")
			call pargr (y1)
			call pargr (-x1a)
		call fprintf (fdb, "%10.4f %10.4f 0.00000\n")
			call pargr (y1)
			call pargr (-x1b)
		call fprintf (fdb, "%10.4f %10.4f 0.00000\n")
			call pargr (y2)
			call pargr (-x2b)
		call fprintf (fdb, "%10.4f %10.4f 0.00000\n")
			call pargr (y2)
			call pargr (-x2a)
		call fprintf (fdb, "%10.4f %10.4f 0.00000\n")
			call pargr (y1)
			call pargr (-x1a)
		call fprintf (fdb, "newrow\n")
	}
	call fprintf (fdb, "# (     -- end of file --     )\n")
end

#
# R3SORT: Sort a list with 3 real items
#

procedure	r3sort (r, npt, r2, r3)

real	r[npt]				# Real vector to sort
int	npt				# number of points
real	r2[npt]				# Additional real vector
real	r3[npt]				# Additional real vector

int	i, j
real	h, h2, h3

begin
# Sort the list (low-to-high)
	do i = 1, npt-1 {
	    do j = 1, npt-i {
		if (r[j] > r[j+1]) {
			h = r[j+1]
			r[j+1] = r[j]
			r[j] = h

			h2 = r2[j+1]
			r2[j+1] = r2[j]
			r2[j] = h2

			h3 = r3[j+1]
			r3[j+1] = r3[j]
			r3[j] = h3
		}
	    }
	}
end

#
# R2SORTD: Sort DOWNWARD a list with 2 real items
#

procedure	r2sortd (r, npt, r2)

real	r[npt]				# Real vector to sort
int	npt				# number of points
real	r2[npt]				# Additional real vector

int	i, j
real	h, h2

begin
# Sort the list (high-to-low)
	do i = 1, npt-1 {
	    do j = 1, npt-i {
		if (r[j] < r[j+1]) {
			h = r[j+1]
			r[j+1] = r[j]
			r[j] = h

			h2 = r2[j+1]
			r2[j+1] = r2[j]
			r2[j] = h2
		}
	    }
	}
end

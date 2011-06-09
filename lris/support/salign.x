# Maskalign, but using simpler file format

include <gset.h>
include <math.h>
include	"futil.h"

define	KEYSFILE	"ucsclris$maskalign.key"
define	PA_ACCEPT 0.05

#
# T_SALIGN -- work out shifts and rotation between two sets of coordinates
# (This is a hacked-up version of t_geoscale constrained to rotations for
# slit-mask alignment)
#
# Created 22-Jul-97 from MASKALIGN



procedure t_salign()

char	coordfile[SZ_FNAME]
real	xrot, yrot
real	pa, pixscale
bool	invert
int	nx, ny
int	niter
real	bxs, def1err, def2err
bool	dcs, dcspa
pointer	fd

char	tchar
char	cmdline[SZ_LINE]			# command string for DCS
int	stat					# command stat
int	npt, ndx, i
real	coeff[2,3], err[3]
#real	x1, x2, xfract
#real	y1, y2, yfract
real	rotposn
pointer	xbuf1, ybuf1, xbuf2, ybuf2
pointer	ebufx1, ebufy1, ebufx2, ebufy2
pointer	rbufx, rbufy, wbuf

real	sina, cosa
real	earcsec, narcsec, eerr, nerr, pdeg, perr
#real	sint, cost, rrot

bool	clgetb()
int	clgeti(), fscan(), nscan(), oscmd()
real	clgetr()
pointer	open()

begin
	call clgstr ("coordfile", coordfile, SZ_FNAME)
	xrot = clgetr ("xrot")
	yrot = clgetr ("yrot")
	nx = clgeti ("nx")
	ny = clgeti ("ny")
	bxs = -2. * clgetr ("box_size")
	niter = clgeti ("niter")
	def1err = clgetr ("def_ref_err")
	def2err = clgetr ("def_err")

        fd = open (coordfile, READ_ONLY, TEXT_FILE)
	npt = 0
	while (fscan(fd) != EOF)
		npt = npt + 1
	call seek (fd, BOF)

	call malloc (xbuf1, npt, TY_REAL)		# ref coord
	call malloc (ybuf1, npt, TY_REAL)
	call malloc (xbuf2, npt, TY_REAL)		# input coord
	call malloc (ybuf2, npt, TY_REAL)
	call malloc (ebufx1, npt, TY_REAL)		# error arrays
	call malloc (ebufy1, npt, TY_REAL)
	call malloc (ebufx2, npt, TY_REAL)
	call malloc (ebufy2, npt, TY_REAL)
	call malloc (rbufx, npt, TY_REAL)		# residuals
	call malloc (rbufy, npt, TY_REAL)
	call malloc (wbuf, npt, TY_REAL)		# weights

	call amovkr (1., Memr[wbuf], npt)

# Read in list of coord pairs and errors (and weights if present)
	ndx = 0
	while (fscan(fd) != EOF) {
		call gargwrd (tchar, 1)
		call reset_scan()
		if (tchar == '#') {
			next
		}
		call gargr (Memr[xbuf2+ndx])
		call gargr (Memr[ybuf2+ndx])
		call gargr (Memr[xbuf1+ndx])
		call gargr (Memr[ybuf1+ndx])
		call gargr (Memr[ebufx1+ndx])
		call gargr (Memr[ebufy1+ndx])
		call gargr (Memr[ebufx2+ndx])
		call gargr (Memr[ebufy2+ndx])
		if (nscan() < 4)
			next
		if (nscan() < 8) {
			Memr[ebufx1+ndx] = INDEF
			Memr[ebufy1+ndx] = INDEF
			Memr[ebufx2+ndx] = INDEF
			Memr[ebufy2+ndx] = INDEF
		}
		ndx = ndx + 1
	}

	call close (fd)
	npt = ndx		# npt reduced by rejection of bad pairs

	do i = 0, npt-1 {
		if (Memr[ebufx1+i] == INDEF || Memr[ebufx1+i] <= 0.)
			Memr[ebufx1+i] = def1err
		if (Memr[ebufy1+i] == INDEF || Memr[ebufy1+i] <= 0.)
			Memr[ebufy1+i] = def1err
		if (Memr[ebufx2+i] == INDEF || Memr[ebufx2+i] <= 0.)
			Memr[ebufx2+i] = def2err
		if (Memr[ebufy2+i] == INDEF || Memr[ebufy2+i] <= 0.)
			Memr[ebufy2+i] = def2err
	}

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
		Memr[rbufx], Memr[rbufy], Memr[wbuf], npt, nx, ny, niter, bxs,
								coeff, err)

	call printf ("#  obj-x  obj-y  slit-x   xres   slit-y   yres     w \n")
	do i = 0, npt-1 {
		call printf ("#%7.2f%7.2f  %7.2f (%5.3f) %7.2f (%5.3f) %6.2f\n")
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

# Now adjust for existing PA:
	pa = clgetr ("zpa") - clgetr ("pa")
	invert = clgetb ("invert")
	pixscale = clgetr ("pixscale")
	if (invert) {
		call eprintf ("Sorry -- not yet tested\n")
		pdeg = RADTODEG( atan2(coeff[1,2], coeff[1,1]))
		perr = RADTODEG(asin(err[1]))
		pa = pa - pdeg			# Small, ...note sign: real PA
		cosa = cos (DEGTORAD(pa))
		sina = sin (DEGTORAD(pa))
		narcsec = pixscale * (-cosa * coeff[1,3] + sina * coeff[2,3])
		earcsec = pixscale * ( sina * coeff[1,3] + cosa * coeff[2,3])
		nerr = pixscale * sqrt ((cosa * err[2])**2 + (sina * err[3])**2)
		eerr = pixscale * sqrt ((sina * err[2])**2 + (cosa * err[3])**2)
	} else {
		pdeg = RADTODEG(-atan2(coeff[1,2], coeff[1,1]))
		perr = RADTODEG(asin(err[1]))
		pa = pa - pdeg			# Small, but...
		cosa = cos (DEGTORAD(pa))
		sina = sin (DEGTORAD(pa))
		narcsec = pixscale * ( cosa * coeff[1,3] + sina * coeff[2,3])
		earcsec = pixscale * (-sina * coeff[1,3] + cosa * coeff[2,3])
		nerr = pixscale * sqrt ((cosa * err[2])**2 + (sina * err[3])**2)
		eerr = pixscale * sqrt ((sina * err[2])**2 + (cosa * err[3])**2)
	}
	
# PA offset is backwards from rotation:
	call printf ("\n====================================================\n")
	Call printf ("\n*** MOVE TELESCOPE/ROTATOR by the following offsets:\n")
	call printf ("\n Offset PA by %6.3f (%5.3f) degree\n")
		call pargr (pdeg)
		call pargr (perr)
	call printf (" Offsets: %7.3f\" EAST (%5.3f)\t %7.3f\" NORTH (%5.3f)\n\n")
		call pargr (earcsec)
		call pargr (eerr)
		call pargr (narcsec)
		call pargr (nerr)
	call printf ("====================================================\n\n")

	dcs = NO
#	dcs = clgetb ("dcs")
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
		rotposn = clgetr ("rotposn")
		call sprintf (cmdline, SZ_LINE,
		  "%s modify -s dcs ROTDEST=%.2f ROTMODE=1")	# (1=pos.angle)
			call pargstr ("rsh manuka")
			call pargr (rotposn+pdeg)
		call eprintf ("\n(sending)  %s \n")
			call pargstr (cmdline)
		stat = oscmd (cmdline)
		if (stat != OK) {
			call eprintf ("command failed!  (%d)\n")
				call pargi (stat)
		}

		call sprintf (cmdline, SZ_LINE, "%s waitfor -s dcs ROTSTAT=8")
			call pargstr ("rsh manuka")	# (8=tracking)
		call eprintf ("\n(sending)  %s \n")
			call pargstr (cmdline)
		stat = oscmd (cmdline)
		if (stat != OK) {
			call eprintf ("command failed!  (%d)\n")
				call pargi (stat)
		}
	    }

# Apply translation
		call sprintf (cmdline, SZ_LINE,
		  "%s modify -s dcs RAOFF=%.2f DECOFF=%.2f REL2CURR=1")
			call pargstr ("rsh manuka")
			call pargr (earcsec)
			call pargr (narcsec)
		call eprintf ("\n(sending)  %s \n")
			call pargstr (cmdline)
		stat = oscmd (cmdline)
		if (stat != OK) {
			call eprintf ("command failed!  (%d)\n")
				call pargi (stat)
		}

		call sprintf (cmdline, SZ_LINE, "%s waitfor -s dcs AXESTAT=64")
			call pargstr ("rsh manuka")	# (64=tracking)
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

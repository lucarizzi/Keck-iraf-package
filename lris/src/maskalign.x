# TBD: work out errors on angle based on RESIDUALS! ditto position errors

include <gset.h>
include <math.h>
include	"futil.h"

define	KEYSFILE	"lris$lib/scr/maskalign.key"
define	PA_ACCEPT 0.027			# this is 0.1" at 1000 px

#
# T_MASKALIGN -- work out shifts and rotation between two sets of coordinates
# (This is a hacked-up version of t_geoscale constrained to rotations for
# slit-mask alignment)
#
# Rev 07-Jul-97:  put in rsh manuka commands to apply offsets directly
# Rev 16-Aug-2000 (GDW) Changed 'manuka' to 'punaluu'


procedure t_maskalign()

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
real	x1, x2, xfract
real	y1, y2, yfract
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
		call gargr (x1)
		call gargr (x2)
		call gargr (xfract)
		Memr[xbuf2+ndx] = x1 + (x2 - x1) * xfract
		call gargr (y1)
		call gargr (y2)
		call gargr (yfract)
		Memr[ybuf2+ndx] = y1 + (y2 - y1) * yfract
		call gargr (x1)
		call gargr (y1)
		if (x1 == INDEF || y1 == INDEF)		# skip indefs
			next
		Memr[xbuf1+ndx] = x1
		Memr[ybuf1+ndx] = y1
		call gargr (Memr[ebufx1+ndx])
		call gargr (Memr[ebufy1+ndx])
		call gargr (Memr[ebufx2+ndx])
		call gargr (Memr[ebufy2+ndx])
		if (nscan() < 8)
			next
		if (nscan() < 12) {
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
		Memr[rbufx], Memr[rbufy], Memr[wbuf], npt, nx, ny, niter,
						xrot, yrot, bxs, coeff, err)

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

# Now adjust for existing PA (convert x,y into E,N):
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
	call printf ("\n Offset PA by %5.2f (%4.2f) degree\n")
		call pargr (pdeg)
		call pargr (perr)
	call printf (" Offsets: %6.2f\" EAST (%4.2f)\t %6.2f\" NORTH (%4.2f)\n\n ")
		call pargr (earcsec)
		call pargr (eerr)
		call pargr (narcsec)
		call pargr (nerr)
	call printf ("====================================================\n\n")

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
		rotposn = clgetr ("rotposn")
		call sprintf (cmdline, SZ_LINE,
		  "%s modify -s dcs ROTDEST=%.2f ROTMODE=1")	# (1=pos.angle)
			call pargstr ("rsh punaluu")
			call pargr (rotposn+pdeg)
		call eprintf ("\n(sending)  %s \n")
			call pargstr (cmdline)
		stat = oscmd (cmdline)
		if (stat != OK) {
			call eprintf ("command failed!  (%d)\n")
				call pargi (stat)
		}

		call sprintf (cmdline, SZ_LINE, "%s waitfor -s dcs ROTSTAT=8")
			call pargstr ("rsh punaluu")	# (8=tracking)
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
			call pargstr ("rsh punaluu")
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
			call pargstr ("rsh punaluu")	# (64=tracking)
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

#
# MASKALIGN -- solve for rotation/translation of two sets of coordinates
# The solution for PA has been modified significantly, 26mar95. It now
# relies on reasonable approximations (eg small angles), and assumes that
# the delta-y are generally significant. There is also a trap for duplicates.
#

procedure maskalign (xref, yref, x, y, xrerr, yrerr, xerr, yerr, xres, yres, w,
				npt, nx, ny, niter, xrot, yrot, bxs, coeff, err)

real	xref[npt], yref[npt]			# reference coord.
real	x[npt], y[npt]				# input coord.
real	xrerr[npt], yrerr[npt], xerr[npt], yerr[npt]	# error vectors
real	xres[npt], yres[npt]			# residual vectors
real	w[npt]					# weights
int	nx, ny					# image size/plot dimensions
int	npt					# vector lengths
int	niter					# no. of iterations in fit
real	xrot, yrot				# im. coord of origin
real	bxs					# box size for unit weight
real	coeff[2,3]				# geotran coeffs
real	err[3]					# errors: ang, xoff, yoff

real	stats[NFITPAR]				# fit info struct

int	i, j, ndx
real	axx, axy, ayx, ayy, bx, by
real	me1x, me1y, bxerr, byerr
real	vmag
real	xf, yf				# fractional location of origin

int	npr
real	delx, dely, delxref, delyref
real	scale, delxlen, rsq, sigysq
real	a, aerr
#real	angerr
pointer	bufr, bufrerr, bufw, buft, buferr
pointer	gp

char	command[32]			# not sure if 32 is good
int	wcs, key
real	wx, wy

int	clgcur(), get_nearest()
real	clgetr()
pointer	gopen()

begin
	npr = npt * (npt - 1) / 2
	call malloc (bufr, npr, TY_REAL)
	call malloc (bufrerr, npr, TY_REAL)
	call malloc (bufw, npr, TY_REAL)

	call malloc (buft, npt, TY_REAL)
	call malloc (buferr, npt, TY_REAL)

	ndx = 0
	do i = 1, npt-1 {
	    do j = i+1, npt {
		delx = x[j] - x[i]
		dely = y[j] - y[i]
		delxref = xref[j] - xref[i]
		delyref = yref[j] - yref[i]
		rsq = delx * delx + dely * dely
		if (rsq == 0.) {
			call eprintf ("WARNING: Duplicate in list!\n")
# Assign zero weight
		    Memr[bufr+ndx] = 0.
		    Memr[bufrerr+ndx] = 1.
		    Memr[bufw+ndx] = 0.
		} else if (dely == 0.) {
			call eprintf ("WARNING: Degenerate case!\n")
# Assign zero weight
		    Memr[bufr+ndx] = 0.
		    Memr[bufrerr+ndx] = 1.
		    Memr[bufw+ndx] = 0.
		} else {
		    scale = sqrt ((delxref*delxref + delyref*delyref) / rsq)
		    delxlen = delx - delxref / scale
		    sigysq = yerr[i] + yerr[j]	# should put in yrerr, too
		    Memr[bufr+ndx] = delxlen / dely
		    Memr[bufrerr+ndx] = (delxlen*delxlen/(dely*dely) * sigysq +
			xerr[i] + xerr[j] + xrerr[i] + xrerr[j]) / (dely*dely)
		    Memr[bufw+ndx] = w[j] * w[i]
		}
		ndx = ndx + 1
	    }
	}

	call get_0lsqf0 (Memr[bufr], Memr[bufrerr], Memr[bufw], npr, stats)
	a = OFFSET[stats]		# tan(pa), radians
	aerr = EOFFSET[stats]		# err (tan(pa)), radians
	a = atan (a)
	axx = cos (a)
	axy = -sin (a)
	ayx = -axy
	ayy = axx

	ndx = 0
#	angerr = aerr * aerr
	do i = 1, npt {
		Memr[buft+ndx] = axx * x[i] + axy * y[i]
		Memr[buferr+ndx] = xerr[i] + axy*axy*yerr[i]
#		Memr[buferr+ndx] = y[i]*y[i]*angerr + xerr[i] + axy*axy*yerr[i]
		ndx = ndx + 1
	}
	call get_lsqf0 (Memr[buft], xref, Memr[buferr], xrerr, w, npt, stats)
	bx = OFFSET[stats]
	bxerr = EOFFSET[stats]
	me1x = ME1[stats]

	ndx = 0
	do i = 1, npt {
		Memr[buft+ndx] = ayx * x[i] + ayy * y[i]
		Memr[buferr+ndx] = yerr[i] + ayx*ayx*xerr[i]
#		Memr[buferr+ndx] = x[i]*x[i]*angerr + yerr[i] + ayx*ayx*xerr[i]
		ndx = ndx + 1
	}
	call get_lsqf0 (Memr[buft], yref, Memr[buferr], yrerr, w, npt, stats)
	by = OFFSET[stats]
	byerr = EOFFSET[stats]
	me1y = ME1[stats]

	do i = 1, npt {
		xres[i] = xref[i] - (axx * x[i] + axy * y[i] + bx)
		yres[i] = yref[i] - (ayx * x[i] + ayy * y[i] + by)
	}

	vmag = 100.
# Open the graphics stream
	gp = gopen ("stdgraph", NEW_FILE, STDGRAPH)

# Fractional location of origin:
	xf = xrot / nx
	yf = yrot / ny

# Plot the data
	call plot_res2d (gp, x, y, xres, yres, w, npt, nx, ny, xf, yf, bxs, vmag)

	while ( clgcur("coord", wx, wy, wcs, key, command, 32) != EOF ) {

	if (key == 'q')
		break

	switch (key) {

	case 'c':
		i = get_nearest (gp, x, y, npt, wx, wy, wcs)
		w[i] = clgetr ("wt")
		if (w[i] == 0.)
			call gmark (gp, x[i], y[i], GM_CROSS, 1., 1.)
		else
			call gmark (gp, x[i], y[i], GM_BOX, w[i]*bxs, w[i]*bxs)

	case 'd':
		i = get_nearest (gp, x, y, npt, wx, wy, wcs)
		w[i] = 0.
		call gmark (gp, x[i], y[i], GM_CROSS, 1., 1.)

	case 'f':
		ndx = 0
		do i = 1, npt-1 {
		    do j = i+1, npt {
			delx = x[j] - x[i]
			dely = y[j] - y[i]
			delxref = xref[j] - xref[i]
			delyref = yref[j] - yref[i]
			rsq = delx * delx + dely * dely
			if (rsq == 0.) {
				call eprintf ("WARNING: Duplicate in list!\n")
# Assign zero weight
			    Memr[bufr+ndx] = 0.
			    Memr[bufrerr+ndx] = 1.
			    Memr[bufw+ndx] = 0.
			} else if (dely == 0.) {
				call eprintf ("WARNING: Degenerate case!\n")
# Assign zero weight
			    Memr[bufr+ndx] = 0.
			    Memr[bufrerr+ndx] = 1.
			    Memr[bufw+ndx] = 0.
			} else {
			    scale = sqrt (
				(delxref*delxref + delyref*delyref) / rsq)
			    delxlen = delx - delxref / scale
			    sigysq = yerr[i] + yerr[j]
			    Memr[bufr+ndx] = delxlen / dely
			    Memr[bufrerr+ndx] = (delxlen*delxlen/(dely*dely) *
				sigysq + xerr[i] + xerr[j] + xrerr[i] +
				xrerr[j]) / (dely*dely)
			    Memr[bufw+ndx] = w[j] * w[i]
			}
			ndx = ndx + 1
		    }
		}

		call get_0lsqf0 (Memr[bufr], Memr[bufrerr], Memr[bufw], npr, stats)
		a = OFFSET[stats]		# tan(pa), radians
		aerr = EOFFSET[stats]		# err (tan(pa)), radians
		a = atan (a)
		axx = cos (a)
		axy = -sin (a)
		ayx = -axy
		ayy = axx

		ndx = 0
		do i = 1, npt {
			Memr[buft+ndx] = axx * x[i] + axy * y[i]
			Memr[buferr+ndx] = xerr[i] + axy*axy*yerr[i]
			ndx = ndx + 1
		}
		call get_lsqf0 (Memr[buft], xref, Memr[buferr], xrerr, w, npt, stats)
		bx = OFFSET[stats]
		bxerr = EOFFSET[stats]
		me1x = ME1[stats]

		ndx = 0
		do i = 1, npt {
			Memr[buft+ndx] = ayx * x[i] + ayy * y[i]
			Memr[buferr+ndx] = yerr[i] + ayx*ayx*xerr[i]
			ndx = ndx + 1
		}
		call get_lsqf0 (Memr[buft], yref, Memr[buferr], yrerr, w, npt, stats)
		by = OFFSET[stats]
		byerr = EOFFSET[stats]
		me1y = ME1[stats]

		do i = 1, npt {
			xres[i] = xref[i] - (axx * x[i] + axy * y[i] + bx)
			yres[i] = yref[i] - (ayx * x[i] + ayy * y[i] + by)
		}
		call gclear (gp)
		call plot_res2d (gp, x, y, xres, yres, w, npt, nx, ny, xf, yf, bxs, vmag)

	case 'm':
		vmag = clgetr ("vmag")
		call gclear (gp)
		call plot_res2d (gp, x, y, xres, yres, w, npt, nx, ny, xf, yf, bxs, vmag)

	case 'I':
		call fatal (0, "INTERRUPT")

        case '?':
                call gpagefile (gp, KEYSFILE, "maskalign cursor commands")

	}
	}
	
	call gclose (gp)

# Print stats again on text plane
	call printf ("#  x-xform: %7.5fx + %7.5fy  + %7.3f  (%5.3f) [%5.2f]\n")
		call pargr (axx)
		call pargr (axy)
		call pargr (bx)
		call pargr (bxerr)
		call pargr (me1x)
	call printf ("#  y-xform: %7.5fx + %7.5fy  + %7.3f  (%5.3f) [%5.2f]\n")
		call pargr (ayx)
		call pargr (ayy)
		call pargr (by)
		call pargr (byerr)
		call pargr (me1y)

	coeff[1,1] = axx
	coeff[1,2] = axy
	coeff[1,3] = bx
	coeff[2,1] = ayx
	coeff[2,2] = ayy
	coeff[2,3] = by
	err[1] = aerr
	err[2] = bxerr
	err[3] = byerr

	call mfree (buferr, TY_REAL)
	call mfree (buft, TY_REAL)
	call mfree (bufw, TY_REAL)
	call mfree (bufrerr, TY_REAL)
	call mfree (bufr, TY_REAL)
	return
end


#
# GET_0LSQF0 -- calculate the weighted mean
# fit equation  <x> = b

	procedure get_0lsqf0 (x, xerr, w, npt, stats)

int	npt				# vector length
real	x[npt]				# input vector
real	xerr[npt]			# error vector
real	w[npt]				# weight vector
real	stats[NFITPAR]			# fit info struct

double	sumxx, sumx, sumw
pointer	bufw

double	dsum1(), dsum2(), dsum3()

begin
	call malloc (bufw, npt, TY_REAL)

	call adivr (w, xerr, Memr[bufw], npt)

	sumxx = dsum3 (x, x, Memr[bufw], npt)
	sumx = dsum2 (x, Memr[bufw], npt)
	sumw = dsum1 (Memr[bufw], npt)

	OFFSET[stats] = sumx / sumw
	if (npt - 1 > 0) {
		ME1[stats] = sqrt (real ((sumxx - sumx*sumx/sumw) / (npt - 1)))
		EOFFSET[stats] = ME1[stats] / sqrt (real (sumw))
	} else {
		ME1[stats] = -1.
		EOFFSET[stats] = -1.
	}

	call mfree (bufw, TY_REAL)
	return
end

int procedure get_nearest (gp, xdata, ydata, ndata, wx, wy, wcs)

pointer	gp
real 	xdata[ARB], ydata[ARB]
int	ndata
real	wx, wy
int	wcs

int	nearest, i
real	ycorr
real	rsq, rsq_min
real	xndc, yndc, xgcndc, ygcndc

real	ggetr()

begin

# need to put in INDEF check 

# Get aspect ratio 
	ycorr = ggetr (gp, "ar")
	if (ycorr == 0.)
		ycorr = 1.

	rsq_min = 2.			# by def'n larger than NDC possible
	nearest = 0

	call gctran (gp, wx, wy, xgcndc, ygcndc, wcs, 0)

	do i = 1, ndata {
		call gctran (gp, xdata[i], ydata[i], xndc, yndc, wcs, 0)
		rsq = (xndc - xgcndc) ** 2 + ( (yndc - ygcndc) * ycorr) ** 2
		if (rsq < rsq_min) {
			rsq_min = rsq
			nearest = i
		}
	}

	if (nearest != 0)
		call gscur (gp, xdata[nearest], ydata[nearest])
	
	return (nearest)
end

#
# PLOT_RES2D -- plot residuals as vector on 2-d surface
#

procedure	 plot_res2d (gp, x, y, xres, yres, w, npt, nx, ny, xf, yf, bxs, 
vmag)

pointer	gp
real	x[npt], y[npt]				# position vectors
real	xres[npt], yres[npt]			# residual vectors
real	w[npt]					# weight vector
int	npt					# length of vectors
int	nx, ny					# plot dimension
real	xf, yf					# fractional loc. of origin
real	bxs					# box size (neg for WCS)
real	vmag					# residual magnification

char	title[SZ_LINE]
int	i
real	wfactor					# box size for full weight

begin
	call gswind (gp, -nx*xf, nx*(1.-xf), -ny*yf, ny*(1.-yf))
	call sprintf (title, SZ_LINE, "Residual Vectors; Box Half-Length= %f pix")
		call pargr (-bxs/2.)
	call glabax (gp, title, "X (pix)", "Y (pix)")

	wfactor = bxs * vmag
	for (i = 1; i <= npt; i = i + 1) {
		if (w[i] == 0.)
			call gmark (gp, x[i], y[i], GM_CROSS, 1., 1.)
		else
			call gmark (gp, x[i], y[i], GM_BOX, w[i]*wfactor,
								w[i]*wfactor)
		call grdraw (gp, xres[i]*vmag, yres[i]*vmag)
	}

	call gflush (gp)

	return
end

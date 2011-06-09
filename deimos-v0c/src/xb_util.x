# Copied from ucsclris/maskalign.x and ucsclris/mboxfind.x
#
# XXX saw_xcorr() has grad added
# XXX TBD: review weighting -- should now be in X
# TBD:  review angle calc
# TBD:	move fitting into separate routine; call initially and with 'f'


include <gset.h>
include <math.h>
include	"futil.h"
include	"instrument.h"		# needed only for det_to_ics
include	"deimos.h"		# needed only for det_to_ics

# define	REL_HT	0.4		# relative height for crossings
# define	ID_CHSZ	9		# Character size of ID string
# define	LAB_SZ 80		# Char size of title/label strings
define		X_HWID	1.5		# x half-width of triangle function
define		Y_HWID	1.5		# y half-width of triangle function
define		FRAC	0.333		# fractional lev. in sorted list for sky

define		KEYSFILE1	"deimos$lib/keys/boxfind.key"
define		KEYSFILE2	"deimos$lib/keys/maskalign.key"

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
#real	scale, delxlen, rsq, sigysq	# variables for Y emphasis
real	scale, delylen, rsq, sigxsq	# variables for X emphasis
real	a, aerr
#real	angerr
pointer	bufr, bufrerr, bufw, buft, buferr
pointer	gp

char	command[32]			# not sure if 32 is good
int	wcs, key
real	wx, wy

int	req_valr(), get_nearest()
int	clgcur()
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
#		} else if (dely == 0.) {
 		} else if (delx == 0.) {
			call eprintf ("WARNING: Degenerate case!\n")
# Assign zero weight
		    Memr[bufr+ndx] = 0.
		    Memr[bufrerr+ndx] = 1.
		    Memr[bufw+ndx] = 0.
		} else {
		    scale = sqrt ((delxref*delxref + delyref*delyref) / rsq)
# emphasize Y spacing:
#		    delxlen = delx - delxref / scale
#		    sigysq = yerr[i] + yerr[j]	# should put in yrerr, too
#		    Memr[bufr+ndx] = delxlen / dely
#		    Memr[bufrerr+ndx] = (delxlen*delxlen/(dely*dely) * sigysq +
#			xerr[i] + xerr[j] + xrerr[i] + xrerr[j]) / (dely*dely)
# emphasize X spacing:
		    delylen = dely - delyref / scale
		    sigxsq = xerr[i] + xerr[j]	# should put in xrerr, too
		    Memr[bufr+ndx] = -delylen / delx
		    Memr[bufrerr+ndx] = (delylen*delylen/(delx*delx) * sigxsq +
			yerr[i] + yerr[j] + yrerr[i] + yrerr[j]) / (delx*delx)
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
#	xf = xrot / nx
#	yf = yrot / ny
	xf = 0.5		# XXX TMP, for DEIMOS
	yf = 0.5		# XXX TMP, for DEIMOS


# Plot the data
	call plot_res2d (gp, x, y, xres, yres, w, npt, nx, ny, xf, yf, bxs, vmag)

	while ( clgcur("coord", wx, wy, wcs, key, command, 32) != EOF ) {

	if (key == 'q')
		break

	switch (key) {

	case 'c':
		i = get_nearest (gp, x, y, npt, wx, wy, wcs)
		if (req_valr ("weight", w[i], w[i], 0., 1.) != OK)
			next
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
# emphasize Y spacing:
#			    delxlen = delx - delxref / scale
#			    sigysq = yerr[i] + yerr[j]
#			    Memr[bufr+ndx] = delxlen / dely
#			    Memr[bufrerr+ndx] = (delxlen*delxlen/(dely*dely) *
#				sigysq + xerr[i] + xerr[j] + xrerr[i] +
#				xrerr[j]) / (dely*dely)
# emphasize X spacing:
			    delylen = dely - delyref / scale
			    sigxsq = xerr[i] + xerr[j]
			    Memr[bufr+ndx] = -delylen / delx
			    Memr[bufrerr+ndx] = (delylen*delylen/(delx*delx) *
				sigxsq + yerr[i] + yerr[j] + yrerr[i] +
				yrerr[j]) / (delx*delx)
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
		if (req_valr ("vmag", vmag, vmag, 1., 200.) != OK)
			next
		call gclear (gp)
		call plot_res2d (gp, x, y, xres, yres, w, npt, nx, ny, xf, yf, bxs, vmag)

	case 'I':
		call fatal (0, "INTERRUPT")

        case '?':
                call gpagefile (gp, KEYSFILE2, "maskalign cursor commands")

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

procedure	 plot_res2d (gp, x, y, xres, yres, w, npt, nx, ny, xf, yf, bxs, vmag)

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

#
# BOX_GRAPH: graph the box profile, find sky, show centers
#

int	procedure box_graph (x, y, zx, zy, tx, ty, nx, ny, title, xb, yb,
					xsz, ysz, xfwhm, yfwhm, cx, cy, zip)

real	x[nx], y[ny]		# position vectors
real	zx[nx], zy[ny]		# Intensity vectors for x,y cuts
real	tx[nx], ty[ny]		# Work      vectors for x,y cuts
int	nx, ny			# size of vectors/array
char	title[ARB]		# title of slit
real	xb, yb			# box centers in x, y
real	xsz, ysz		# box sizes in x, y
real	xfwhm, yfwhm		# x, y FWHM of stars
real	cx, cy			# returned x,y of star
int	zip			# "zip" mode (no interaction)

int	noff, i, ni
real	sky
real	gymin, gymax
real	g1xmin, g1xmax, g1ymin, g1ymax
real	g2xmin, g2xmax, g2ymin, g2ymax
real	grad

char	command[32]			# not sure if 32 is good
int	wcs, key
real	wx, wy
real	ndcx, ndcy
pointer	gp, gp1, gp2

int	req_valr()
real	saw_xcorr()
int	clgcur()
pointer	gopen()

begin
# Open the graphics stream
	gp = gopen ("stdgraph", NEW_FILE, STDGRAPH)
	gp1 = gopen ("stdgraph", APPEND, STDGRAPH)
	gp2 = gopen ("stdgraph", APPEND, STDGRAPH)
	call gsview (gp1, 0.07, 0.5, 0.2, 0.8)
	call gsview (gp2, 0.57, 1.0, 0.2, 0.8)

# Get Y-limits (same in both plots)
	call alimr (zx, nx, g1ymin, g1ymax)
	call alimr (zy, ny, g2ymin, g2ymax)
	g1ymin = min (g1ymin, g2ymin)
	g1ymax = max (g1ymax, g2ymax)
	gymin = g1ymin - 0.05 * (g1ymax - g1ymin)
	gymax = g1ymax + 0.05 * (g1ymax - g1ymin)
# Individual X-limits
	call alimr (x, nx, g1xmin, g1xmax)
	call alimr (y, ny, g2xmin, g2xmax)

# Plot the data
	call gclear (gp)
# First X:
	call gseti (gp1, G_FRAMECOLOR, 0)
	call gswind (gp1, g1xmin, g1xmax, gymin, gymax)
	call glabax (gp1, title, "x(pix)", "box profile")
	call gpline (gp1, x, zx, nx)
	call gmark (gp1, xb, g1ymin, GM_BOX, -xsz, 3.0)
	call gtext (gp1, xb, g1ymin, "BOX", "h=c;v=c;q=h;s=0.6")
	call gflush (gp1)

# now Y:
	call gseti (gp2, G_FRAMECOLOR, 0)
	call gswind (gp2, g2xmin, g2xmax, gymin, gymax)
	call glabax (gp2, title, "y(pix)", "")
	call gpline (gp2, y, zy, ny)
	call gmark (gp2, yb, g1ymin, GM_BOX, -ysz, 3.0)
	call gtext (gp2, yb, g1ymin, "BOX", "h=c;v=c;q=h;s=0.6")
	call gflush (gp2)

# Do an initial estimate of sky level (we _assume_ profile is centered!):
	i  = 0.5 * (nx - xsz) + 1.
	ni = 0.5 * (nx + xsz) + 1. - i + 1
	call asrtr (zx[i], tx, ni)
	i = max (ni * FRAC + 0.5, 1.)
	wx = tx[i]
	i  = 0.5 * (ny - ysz) + 1.
	ni = 0.5 * (ny + ysz) + 1. - i + 1
	call asrtr (zy[i], ty, ni)
	i = max (ni * FRAC + 0.5, 1.)
	wy = ty[i]
	sky = min (wx, wy)
# ... and copied from below ...
			call gseti (gp1, G_PLTYPE, GL_DASHED)
			call gamove (gp1, g1xmin, sky)
			call gadraw (gp1, g1xmax, sky)
			call gseti (gp1, G_PLTYPE, GL_SOLID)
			call gseti (gp2, G_PLTYPE, GL_DASHED)
			call gamove (gp2, g2xmin, sky)
			call gadraw (gp2, g2xmax, sky)
			call gseti (gp2, G_PLTYPE, GL_SOLID)
			call amaxkr (zx, sky, tx, nx)
			noff = xsz / 2. + 0.5
	                cx = saw_xcorr (x, tx, nx, xfwhm, xfwhm/2., noff, grad)
			call amaxkr (zy, sky, ty, ny)
			noff = ysz / 2. + 0.5
	                cy = saw_xcorr (y, ty, ny, yfwhm, yfwhm/2., noff, grad)
			call gamove (gp1, cx, gymin)
			call gadraw (gp1, cx, gymax)
			call gamove (gp2, cy, gymin)
			call gadraw (gp2, cy, gymax)

#	sky = INDEF
#	cx = 0.
#	cy = 0.

	if (zip == NO) {

	  while ( clgcur("coord", wx, wy, wcs, key, command, 32) != EOF ) {

	    if (key == 'q') {
		if (sky == INDEF)
			call printf ("Must set sky level (SPACE BAR)")
		else
			break
	    }

	    switch (key) {

		case ' ':
# remove old marks
			call gseti (gp1, G_PLTYPE, GL_CLEAR)
			call gamove (gp1, cx, gymin)
			call gadraw (gp1, cx, gymax)
			call gseti (gp1, G_PLTYPE, GL_SOLID)
			call gseti (gp2, G_PLTYPE, GL_CLEAR)
			call gamove (gp2, cy, gymin)
			call gadraw (gp2, cy, gymax)
			call gseti (gp2, G_PLTYPE, GL_SOLID)
			sky = wy
			call gseti (gp1, G_PLTYPE, GL_DASHED)
			call gamove (gp1, g1xmin, sky)
			call gadraw (gp1, g1xmax, sky)
			call gseti (gp1, G_PLTYPE, GL_SOLID)
			call gseti (gp2, G_PLTYPE, GL_DASHED)
			call gamove (gp2, g2xmin, sky)
			call gadraw (gp2, g2xmax, sky)
			call gseti (gp2, G_PLTYPE, GL_SOLID)
			call amaxkr (zx, sky, tx, nx)
			noff = xsz / 2. + 0.5
	                cx = saw_xcorr (x, tx, nx, xfwhm, xfwhm/2., noff, grad)
			call amaxkr (zy, sky, ty, ny)
			noff = ysz / 2. + 0.5
	                cy = saw_xcorr (y, ty, ny, yfwhm, yfwhm/2., noff, grad)
			call gamove (gp1, cx, gymin)
			call gadraw (gp1, cx, gymax)
			call gamove (gp2, cy, gymin)
			call gadraw (gp2, cy, gymax)

		case 'r':
			call gclear (gp)
# first X:
			call gswind (gp1, g1xmin, g1xmax, gymin, gymax)
			call glabax (gp1, title, "x", "box profile")
			call gpline (gp1, x, zx, nx)
			call gmark (gp1, xb, g1ymin, GM_BOX, -xsz, 3.0)
			call gtext (gp1, xb, g1ymin, "BOX", "h=c;v=c;q=h;s=0.6")
			if (sky != INDEF) {
				call gamove (gp1, cx, gymin)
				call gadraw (gp1, cx, gymax)
			}
			call gflush (gp1)
# then Y:
			call gswind (gp2, g2xmin, g2xmax, gymin, gymax)
			call glabax (gp2, title, "y", "")
			call gpline (gp2, y, zy, ny)
			call gmark (gp2, yb, g1ymin, GM_BOX, -ysz, 3.0)
			call gtext (gp2, yb, g1ymin, "BOX", "h=c;v=c;q=h;s=0.6")
			if (sky != INDEF) {
				call gamove (gp2, cy, gymin)
				call gadraw (gp2, cy, gymax)
			}
			call gflush (gp2)

		case 'f':
			call gctran (gp2, wx, wy, ndcx, ndcy, 1, 0)
			if (ndcx < 0.5) {
				call gseti (gp1, G_PLTYPE, GL_CLEAR)
				call gamove (gp1, cx, gymin)
				call gadraw (gp1, cx, gymax)
				call gseti (gp1, G_PLTYPE, GL_SOLID)
				call gctran (gp1, ndcx, ndcy, wx, wy, 0, 1)
				cx = wx
				call gamove (gp1, cx, gymin)
				call gadraw (gp1, cx, gymax)
			} else {
				call gseti (gp2, G_PLTYPE, GL_CLEAR)
				call gamove (gp2, cy, gymin)
				call gadraw (gp2, cy, gymax)
				call gseti (gp2, G_PLTYPE, GL_SOLID)
				call gctran (gp2, ndcx, ndcy, wx, wy, 0, 1)
				cy = wx
				call gamove (gp2, cy, gymin)
				call gadraw (gp2, cy, gymax)
			}

		case 'w':
			if (req_valr ("x_fwhm_px", xfwhm, xfwhm, 2., 20.) != OK)
				next
			if (req_valr ("y_fwhm_px", yfwhm, yfwhm, 2., 20.) != OK)
				next

		case 'i':
			call gclose (gp)
			return (ERR)
			break

		case 'z':
			zip = YES

		case '?':
		    call gpagefile (gp, KEYSFILE1, "boxfind cursor commands")

		case 'I':
			call fatal (0, "INTERRUPT")
	    }
	  }
	}

	call gclose (gp)
	return (OK)
end

int	procedure box_center (x, y, z, nx, ny, zx, zy, xsz, ysz, cx, cy)

real	x[nx], y[ny]		# position vectors
real	z[nx,ny]		# intensity array
int	nx, ny			# size of vectors/array
real	zx[nx], zy[ny]		# Work vectors for x,y cuts
real	xsz, ysz		# pixel sizes of boxes
real	cx, cy			# returned centers
real	grad			# new gradient, currently unused

int	i, j, i1, i2, j1, j2
int	nxrad, nyrad		# variable search radii

real	saw_xcorr()
begin

# Get the initial profile
	call amovkr (0., zx, nx)
	call amovkr (0., zy, ny)
	j1 = 0.5 * ny - 1
	j2 = j1 + 4
	do j = j1, j2 {
		call aaddr (zx, z[1,j], zx, nx)
	}
	i1 = 0.5 * nx - 1
	i2 = i1 + 4
	do j = 1, ny {
	    do i = i1, i2 {
		zy[j] = zy[j] + z[i,j]
	    }
	}
	call amulkr (zx, 1./real(j2-j1+1), zx, nx)
	call amulkr (zy, 1./real(i2-i1+1), zy, ny)

# calculate the maximum sensible search range:
	nxrad = 0.5 * (nx - xsz - 4. * X_HWID) - 0.5
	nyrad = 0.5 * (ny - ysz - 4. * Y_HWID) - 0.5

	cx = saw_xcorr (x, zx, nx, xsz, X_HWID, nxrad, grad)
	if (cx ==INDEF) {
		call eprintf ("Can't find box at %4.0f,%4.0f\n")
			call pargr (x[nx/2])
			call pargr (y[ny/2])
	}
	cy = saw_xcorr (y, zy, ny, ysz, Y_HWID, nyrad, grad)
	if (cy ==INDEF) {
		call eprintf ("Can't find box at %4.0f,%4.0f -- check positions on image and prescan value!\n")
			call pargr (x[nx/2])
			call pargr (y[ny/2])
	}

	if (cx == INDEF || cy == INDEF)
		return (ERR)

# Get the final profile:
	call amovkr (0., zx, nx)
	call amovkr (0., zy, ny)
	j1 = cy - y[1] + 1 - 0.5 * xsz + X_HWID
	j2 = cy - y[1] + 1 + 0.5 * xsz - X_HWID + 0.5
	do j = j1, j2 {
		call aaddr (zx, z[1,j], zx, nx)
	}
	i1 = cx - x[1] + 1 - 0.5 * xsz + X_HWID
	i2 = cx - x[1] + 1 + 0.5 * xsz - X_HWID + 0.5
	do j = 1, ny {
	    do i = i1, i2 {
		zy[j] = zy[j] + z[i,j]
	    }
	}
	call amulkr (zx, 1./real(j2-j1+1), zx, nx)
	call amulkr (zy, 1./real(i2-i1+1), zy, ny)

	cx = saw_xcorr (x, zx, nx, xsz, X_HWID, nxrad, grad)
	if (cx ==INDEF)
		call eprintf ("Got lost on recenter; check input positions")
	cy = saw_xcorr (y, zy, ny, ysz, Y_HWID, nyrad, grad)
	if (cy ==INDEF)
		call eprintf ("Got lost on recenter; check input positions")

	if (cx == INDEF || cy == INDEF)
		return (ERR)

# Print out centers for debugging
# call eprintf ("cx,cy final = %6.2f %6.2f\n")
# call pargr (cx)
# call pargr (cy)

	return (OK)
end

#
# SAW_XCORR: xcorr's a modified sawtooth with a vector to find centers
# This subroutine ALSO appears in MBOXFIND
#
# NOTE -- there can/will be an error if "sz-2*hwid" is larger than feature.
#
# Note -- as this stands, it is still not quite correct. If there are multiple
# crossings, the check on the slope SHOULD BE over the width of the feature;
# instead it is local only. This will probably produce the desired result in
# realistic cases.

real	procedure saw_xcorr (xgrid, z, nx, sz, hwid, noff, grad)

real	xgrid[nx], z[nx]		# position and intensity vectors
int	nx				# size of above
real	sz				# spacing of the peaks
real	hwid				# half-width of the peaks (pixels)
int	noff				# pixels to correlate across
real	grad				# returned gradient

int	nf, nc, xcen
int	ipos, ineg
int	i
real	xpeak, x, xdiff
real	maxderiv, cpos, cneg, fractpos
pointer	form, xcorr, deriv

real	adotr()

begin
# Create "form" vector; this is offset pos and neg. triangle functions
#          +
#         + +
#  +++++++   +++++++++   +++++++++.
#                     + +
#                      +
#
	nf = nx + 2 * noff
	nc = 1  + 2 * noff
	call calloc (form, nf, TY_REAL)
	call calloc (xcorr, nc, TY_REAL)
	call calloc (deriv, nc, TY_REAL)

	xcen = nf / 2 			# Center (to low ) of vectors
	xpeak = 0.5 * sz
	do i = 0, nf-1 {
#		x = i + 1 - (xcen + noff) + 1	# relative offset
		x = nf - (i+1) - xcen		# relative offset
		xdiff = abs (abs (x) - xpeak)
		if (x > 0.)
			Memr[form+i] =  max ((1. - xdiff / hwid), 0.)
		else
			Memr[form+i] = -max ((1. - xdiff / hwid), 0.)
	}

# ready to x-corr; this is actually backwards, but since "form" is symmetric we
# may cheat and assign the negative value:

	do i = 0, 2*noff
		Memr[xcorr+i] = adotr (z, Memr[form+2*noff-i], nx)
	do i = 0, 2*noff-1
		Memr[deriv+i] = 0.5 * (Memr[xcorr+i+1] - Memr[xcorr+i-1])
	Memr[deriv+2*noff] =  Memr[xcorr+2*noff] - Memr[xcorr+2*noff-1]
	Memr[deriv]        =  Memr[xcorr+1]      - Memr[xcorr]

	maxderiv = 0.
	cneg = INDEF
	do i = 0, 2*noff-1 {
		if (Memr[xcorr+i] <= 0. && Memr[xcorr+i+1] > 0.) {
			if (Memr[deriv+i] > maxderiv) {
				ineg = noff - i
				ipos = ineg - 1
				cneg = Memr[xcorr+i]
				cpos = Memr[xcorr+i+1]
				maxderiv = Memr[deriv+i]
			}
		}
	}
	if (cneg == INDEF) {
		x = INDEF
		return (x)
	}

	fractpos = cpos / (cpos - cneg)
	xcen = (nx+1) /2
	x = fractpos * xgrid[xcen-ineg] + (1. - fractpos) * xgrid[xcen-ipos]

	call mfree (deriv, TY_REAL)
	call mfree (xcorr, TY_REAL)
	call mfree (form, TY_REAL)

# TMP! (trial)
	grad = maxderiv

	return (x)	
end

############################# MOVE ELSEWHERE ######################

#
# DET_TO_ICS: convert "detector" coords into ICS coords
#

procedure	det_to_ics (xdet, ydet, ccd, xics, yics)

real	xdet, ydet				# x,y "detector" coords
real	ccd[NCCD,3]				# CCD geometry
real	xics, yics				# x,y ICS (returned)

int	n
real	xccd, yccd		# coords w/in indv ccd, normal orientation

int	det_chip ()
begin
	n = det_chip (xdet, ydet)

	xccd = xdet - CCDXPIX * (n-1)

	if (ydet > CCDYPIX) {
		yccd = ydet - CCDYPIX
	} else {
		yccd = ydet
	}

	call ccd_to_ics (xccd, yccd, ccd, n, xics, yics)
end

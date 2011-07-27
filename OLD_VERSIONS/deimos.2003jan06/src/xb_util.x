# Copied from ucsclris/maskalign.x and ucsclris/mboxfind.x
#
# XXX saw_xcorr() has grad added
# XXX TBD: review weighting -- should now be in X

# TBD:  review angle calc
# TBD:	move fitting into separate routine; call initially and with 'f'


include <gset.h>
include <math.h>
include	"futil.h"

# define	REL_HT	0.4		# relative height for crossings
# define	ID_CHSZ	9		# Character size of ID string
# define	LAB_SZ 80		# Char size of title/label strings
# define	X_HWID	1.5		# x half-width of triangle function
# define	Y_HWID	1.5		# y half-width of triangle function
# define	X_CRAD	13		# x centering radius
# define	Y_CRAD	13		# y centering radius
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

int	clgcur(), get_nearest()
real	promptr()
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
#### DEBUG
# do i = 1, npt {
# call eprintf ("DEBUG %6f %6f\n")
# call pargr (Memr[buferr+i-1])
# call pargr (yerr(i))
# }
	call get_lsqf0 (Memr[buft], yref, Memr[buferr], yrerr, w, npt, stats)
	by = OFFSET[stats]
	byerr = EOFFSET[stats]
	me1y = ME1[stats]
call eprintf ("DEBUG %6f %6f\n")
call pargr (byerr)
call pargr (me1y)

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
		w[i] = promptr( "Enter weight", SZ_LINE, w[i], 0.0, 1.0)
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
#			} else if (dely == 0.) {
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
		vmag = promptr( "Enter magnification factor", SZ_LINE, vmag, 0., INDEF)
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

bool procedure	box_graph (x, y, zx, zy, tx, ty, nx, ny, title, xb, yb,
						xsz, ysz, xfwhm, yfwhm, cx, cy)

real	x[nx], y[ny]		# position vectors
real	zx[nx], zy[ny]		# Intensity vectors for x,y cuts
real	tx[nx], ty[ny]		# Work      vectors for x,y cuts
int	nx, ny			# size of vectors/array
char	title[ARB]		# title of slit
real	xb, yb			# box centers in x, y
real	xsz, ysz		# box sizes in x, y
real	xfwhm, yfwhm		# x, y FWHM of stars
real	cx, cy			# returned x,y of star

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

int	clgcur()
real	promptr()
real	saw_xcorr()
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
			xfwhm = promptr( "Enter new x FWHM [pix]", SZ_LINE, xfwhm, 1., INDEF)
			yfwhm = promptr( "Enter new y FWHM [pix]", SZ_LINE, yfwhm, 1., INDEF)

		case 'i':
			call gclose (gp)
			return ((!OK))
			break

		case '?':
		    call gpagefile (gp, KEYSFILE1, "boxfind cursor commands")

		case 'I':
			call fatal (0, "INTERRUPT")
	    }
	}

	call gclose (gp)
	return (OK)
end

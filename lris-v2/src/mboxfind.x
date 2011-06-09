include	<imhdr.h>
include <gset.h>
include <math.h>

define		REL_HT	0.4		# relative height for crossings
define		ID_CHSZ	9		# Character size of ID string
define		LAB_SZ 80		# Char size of title/label strings
define		X_HWID	1.5		# x half-width of triangle function
define		Y_HWID	1.5		# y half-width of triangle function
define		X_CRAD	13		# x centering radius
define		Y_CRAD	13		# y centering radius
define		FRAC	0.333		# fractional lev. in sorted list for sky

define		PRE_COL 21		# Number of prescan columns per amp

define	KEYSFILE	"lris$lib/scr/boxfind.key"

#
# MBOXFIND: locate boxes and stars within them
#
# Since the boxes are relatively small, read in section and deal with the 
# 2-D image rather than cuts
#
# We currently have a kludge for robustness -- runs box-center twice, the
# second time after recentering image section on initial box center

procedure t_mboxfind()

char	image[SZ_FNAME]
char	input[SZ_FNAME]
char	output[SZ_FNAME]
pointer	im
int	inprepix			# prescan pix included in input x-coord.
int	nxbox, nybox			# Initial box size
real	xsz, ysz
real	xfwhm, yfwhm
real	xoff, yoff			# x,y shifts to apply to input box coord
pointer	fda, fdb

bool	bxstat
char	tchar, idstr[ID_CHSZ]
int	namp, prepix
int	npts, ndx, i, j
int	ncols, nlines
int	nx, ny
real 	xs, ys
real	xb, yb, xstar, ystar
int 	x1, y1, x2, y2
pointer	bufx, bufy, bufzx, bufzy, buftx, bufty
pointer	buf

real	rot, el, xflex, yflex		# parameters for flexure correction

bool	box_graph()
bool	imaccf()
int	clgeti(), imgeti()
int	fscan(), nscan()
real	clgetr(), imgetr()
pointer	immap(), open(), imgs2r()

begin
	call clgstr ("image", image, SZ_FNAME)
	nxbox = clgeti ("nxbox")
	nybox = clgeti ("nybox")
	xsz = clgetr ("xsz")
	ysz = clgetr ("ysz")
	xfwhm = clgetr ("xfwhm")
	yfwhm = clgetr ("yfwhm")
	xoff  = clgetr ("xoff")
	yoff  = clgetr ("yoff")
	inprepix = clgeti ("prescan")

	im = immap (image, READ_ONLY, 0)
	ncols = IM_LEN(im,1)
	nlines = IM_LEN(im,2)

	call clgstr ("input", input, SZ_FNAME)
	fda = open (input, READ_ONLY, TEXT_FILE)
	call clgstr ("output", output, SZ_FNAME)
	fdb = open (output, NEW_FILE, TEXT_FILE)

# Allocate arrays for marginal plots
	call malloc (bufx, nxbox, TY_REAL)
	call malloc (bufy, nybox, TY_REAL)
	call malloc (bufzx, nxbox, TY_REAL)
	call malloc (bufzy, nybox, TY_REAL)
	call malloc (buftx, nxbox, TY_REAL)
	call malloc (bufty, nybox, TY_REAL)

# Get relevant keywords ...
# ... for the number of prescan columns     # ?? FIX - need string
	if (imaccf (im, "NUMAMPS")) {
		namp = imgeti (im, "NUMAMPS")
	} else {
		call eprintf ("NUMAMPS missing; ONE AMP assumed!\n")
		namp = 1
	}

	if (imaccf (im, "PREPIX")) {
		prepix = namp * imgetr (im, "PREPIX")
	} else {
		call eprintf ("PREPIX missing; %2d assumed!\n")
			call pargi (PRE_COL)
		prepix = namp * PRE_COL
	}

# ... for the rotator angle and elevation for flexure mapping
	if (imaccf (im, "ROTPPOSN")) {
		rot = imgetr (im, "ROTPPOSN")
	} else {
		call eprintf ("ROTPPOSN missing; 90 (LRIS stow) assumed!\n")
		rot = 90.
	}

	if (imaccf (im, "EL")) {
		el = imgetr (im, "EL")
	} else {
		call eprintf ("EL missing; 0 (LRIS stow) assumed!\n")
		el = 0.
	}
	call flex_corr (rot, el, xflex, yflex)

# Print out summary of image:
	call printf ("  Prescan col.   flexure     offset     Total\n")
	call printf ("X:     %3d         %5.1f      %5.1f     %5.1f\n")
		call pargi (prepix)
		call pargr (xflex)
		call pargr (xoff)
		call pargr (prepix+xflex+xoff)
	call printf ("Y:     %3d         %5.1f      %5.1f     %5.1f\n\n")
		call pargi (0)
		call pargr (yflex)
		call pargr (yoff)
		call pargr (yflex+yoff)

	xoff = prepix+xflex+xoff
	yoff = yflex+yoff

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

# adjust x-values for differences in prescan pixels, both for offsets:
		xs = xs + xoff - inprepix
		ys = ys + yoff

		call gargwrd (idstr, ID_CHSZ)
		if (nscan() < 3)
			call strcpy ("(no ID)", idstr, ID_CHSZ)

# Now, we do a kludge -- run this twice, recentered on box the second time
	    xb = xs
	    yb = ys
	    do j = 1, 2 {
		x1 = xb - nxbox/2
		x2 = x1 + nxbox - 1
		y1 = yb - nybox/2
		y2 = y1 + nybox - 1

# checks on out-of-bounds
		x1 = max (1, x1)
		x2 = min (ncols, x2)
		y1 = max (1, y1)
		y2 = min (nlines, y2)

# Actual box length 
		nx = x2 - x1 + 1
		ny = y2 - y1 + 1

# Get the image section
		buf = imgs2r (im, x1, x2, y1, y2)

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
# Second time through loop with recentering on box
	    }

	    call printf ("Box center:  %6.2f %6.2f  (%4.1fx,%4.1fy removed) (del:%4.1f,%4.1f)\n")
		call pargr (xb-xoff)
		call pargr (yb-yoff)
		call pargr (xoff)
		call pargr (yoff)
		call pargr (xb-xs)
		call pargr (yb-ys)


# Now get star location
		bxstat = box_graph (Memr[bufx], Memr[bufy], Memr[bufzx],
			Memr[bufzy], Memr[buftx], Memr[bufty], nx, ny, idstr,
				xb, yb, xsz, ysz, xfwhm, yfwhm, xstar, ystar)

# Write out box/star coordinates
		if (bxstat == OK) {
		    call fprintf (fdb,
			"%7.2f %7.2f 0.5    %7.2f %7.2f  0.5   %7.2f  %7.2f\n")
			call pargr (xb-0.5*xsz-prepix)
			call pargr (xb+0.5*xsz-prepix)
			call pargr (yb-0.5*ysz)
			call pargr (yb+0.5*ysz)
			call pargr (xstar-prepix)
			call pargr (ystar)
		} else {
		    call fprintf (fdb, "# Skipped\n")
		}
		ndx = ndx + 1
	}
	npts = ndx

# Close up
	call imunmap (im)
	call mfree (bufty, TY_REAL)
	call mfree (buftx, TY_REAL)
	call mfree (bufzy, TY_REAL)
	call mfree (bufzx, TY_REAL)
	call mfree (bufy, TY_REAL)
	call mfree (bufx, TY_REAL)
	call close (fda)
	call close (fdb)
end

procedure	box_center (x, y, z, nx, ny, zx, zy, xsz, ysz, cx, cy)

real	x[nx], y[ny]		# position vectors
real	z[nx,ny]		# intensity array
int	nx, ny			# size of vectors/array
real	zx[nx], zy[ny]		# Work vectors for x,y cuts
real	xsz, ysz		# pixel sizes of boxes
real	cx, cy			# returned centers

int	i, j, i1, i2, j1, j2

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

	cx = saw_xcorr (x, zx, nx, xsz, X_HWID, X_CRAD)
	if (cx ==INDEF) {
		call eprintf ("Can't find box at %4.0f,%4.0f -- check positions on image and prescan value!\n")
			call pargr (x[nx/2])
			call pargr (y[ny/2])
		call fatal (0, "Can't find box-x -- check positions on image!")
	}
	cy = saw_xcorr (y, zy, ny, ysz, Y_HWID, Y_CRAD)
	if (cy ==INDEF) {
		call eprintf ("Can't find box at %4.0f,%4.0f -- check positions on image and prescan value!\n")
			call pargr (x[nx/2])
			call pargr (y[ny/2])
		call fatal (0, "Can't find box-y -- check positions on image!")
	}

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

	cx = saw_xcorr (x, zx, nx, xsz, X_HWID, X_CRAD)
	if (cx ==INDEF)
		call fatal (0, "Got lost on recenter; check input positions")
	cy = saw_xcorr (y, zy, ny, ysz, Y_HWID, Y_CRAD)
	if (cy ==INDEF)
		call fatal (0, "Got lost on recenter; check input positions")

# Print out centers for debugging
# call eprintf ("cx,cy final = %6.2f %6.2f\n")
# call pargr (cx)
# call pargr (cy)

end

#
# SAW_XCORR: xcorr's a modified sawtooth with a vector to find centers
#
# NOTE -- there can/will be an error if "sz-2*hwid" is larger than feature.
#
# Note -- as this stands, it is still not quite correct. If there are multiple
# crossings, the check on the slope SHOULD BE over the width of the feature;
# instead it is local only. This will probably produce the desired result in
# realistic cases.

real	procedure saw_xcorr (xgrid, z, nx, sz, hwid, noff)

real	xgrid[nx], z[nx]		# position and intensity vectors
int	nx				# size of above
real	sz				# spacing of the peaks
real	hwid				# half-width of the peaks
int	noff				# pixels to correlate across

int	nf, nc, xcen
int	ipos, ineg
int	i
real	xpeak, x, xdiff
real	maxderiv, cpos, cneg, fractpos
pointer	form, xcorr, deriv

real	vsum2()

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
		Memr[xcorr+i] = vsum2 (z, Memr[form+2*noff-i], nx)
	do i = 0, 2*noff-1
		Memr[deriv+i] = 0.5 * (Memr[xcorr+i+1] - Memr[xcorr+i-1])
	Memr[deriv+2*noff] =  Memr[xcorr+2*noff] - Memr[xcorr+2*noff-1]
	Memr[deriv]        =  Memr[xcorr+1]      - Memr[xcorr]

	maxderiv = 0.
	cneg = INDEF
	do i = 0, 2*noff-1 {
#call printf ("%3d  %6.1f %6.2f / %3d %3d\n")
#call pargi (i)
#call pargr (Memr[xcorr+i])
#call pargr (Memr[deriv+i])
#call pargi (ipos)
#call pargi (ineg)
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

	return (x)	
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

char	command[32]			# not sure if 32 is good
int	wcs, key
real	wx, wy
real	ndcx, ndcy
pointer	gp, gp1, gp2

int	clgcur()
real	clgetr()
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
	                cx = saw_xcorr (x, tx, nx, xfwhm, xfwhm/2., noff)
			call amaxkr (zy, sky, ty, ny)
			noff = ysz / 2. + 0.5
	                cy = saw_xcorr (y, ty, ny, yfwhm, yfwhm/2., noff)
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
	                cx = saw_xcorr (x, tx, nx, xfwhm, xfwhm/2., noff)
			call amaxkr (zy, sky, ty, ny)
			noff = ysz / 2. + 0.5
	                cy = saw_xcorr (y, ty, ny, yfwhm, yfwhm/2., noff)
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
			xfwhm = clgetr ("new_x_fwhm")
			yfwhm = clgetr ("new_y_fwhm")

		case 'i':
			call gclose (gp)
			return ((!OK))
			break

		case '?':
		    call gpagefile (gp, KEYSFILE, "boxfind cursor commands")

		case 'I':
			call fatal (0, "INTERRUPT")
	    }
	}

	call gclose (gp)
	return (OK)
end

#
# FLEX_CORR:  work out the flexure correction to the standand LRIS stow
# position (ROTPPOSN=90, EL=0).
# Currently defined values give 0.10 and 0.17 pixel rms in x and y.
#

define	AMP_X	-3.45		# pixel amplitude in x
define	AMPS_X	-0.4		# pixel amplitude in x (sinel, unmodulated)
define	OFF_X	 0.0		# pixel offset in x
define	PHI_X	 90.		# angle zeropoint in x
define	AMP_Y	-3.5		# pixel amplitude in y
define	AMPS_Y	 1.1		# pixel amplitude in y (sinel, unmodulated)
define	OFF_Y	 3.4		# pixel offset in y
define	PHI_Y	 0.		# angle zeropoint in y
define	RAMP0_Y	 1.2		# rotator felx. pixel amplitude in y
define	RAMP1_Y	 -0.4		# rotator felx. pixel amplitude in y
define	RPHI_Y	 5.		# rotator flex. angle zeropoint in y

procedure	flex_corr (rot, el, delx, dely)

real	rot				# rotator angle
real	el				# elevation
real	delx, dely			# returned offsets, in pixels

real	cosel, sinel			# cos, sin of elevation angle

begin
	cosel = cos (DEGTORAD(el))
	sinel = sin (DEGTORAD(el))

	delx =  AMP_X * cosel * sin (DEGTORAD (rot - PHI_X)) + OFF_X +
		AMPS_X * sinel

	dely =  AMP_Y * cosel * sin (DEGTORAD (rot - PHI_Y)) + OFF_Y +
		(RAMP0_Y + RAMP1_Y * cosel) * sin (2.*DEGTORAD (rot - RPHI_Y)) +
		AMPS_Y * sinel
end

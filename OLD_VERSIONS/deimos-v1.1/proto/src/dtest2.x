# NB -- **MUST** address some of the TMPs and FIXes!
# Makes MUCH more sense to use rtest routines for finding slits ...
#      .... although these may not identify a-boxes
# Need to work out anamorphic factor
# Need to put in adjustments for DWXL8
# Should put in weighting scheme.

include <math.h>
include <math/iminterp.h>
include "instrument.h"
include "align.h"		# used for finding boxes

define		NORD	5	# "order" of fit to derivatives
define	DEF_PROP_FILE	"deimos$lib/prop/redprop.txt"	# default properties

# Grating/focus params (always solved)
define	GX	1
define	GY	2
define	GZ	3
define	CF	4

define	NVAR	4

#
# QLSOLVE: Quick, limited solve ...
#

procedure qltest (fname, ccd, sys, map)

char	fname[SZ_FNAME]			# name of image
real	ccd[NCCD,3]			# CCD geometry
double	sys[NPARAM]			# system params
pointer	map[8]				# pointers to surf fits (1,2=amap;
					# 3,4=b; 5,6 unused, 7,8=brmap

int	nmeas				# No. measurements
pointer	bufxmm				# pointer to slit x
pointer	bufymm				# pointer to slit y 
pointer	bufwav				# pointer to wavelength, order
pointer	bufx, bufy			# pointers to x,y (pixel) meas'd values
pointer	bufxp, bufyp			# pointers to x,y (pixel) predicted
pointer	bufmat, bufpdx, bufpdy		# pointers to matrix array, partials
pointer	bufdel				# pointer to del-corrections

int	nv
int	iord
int	niter
int	ndx
real	rmsx, rmsy

#  these are needed for locating features
char	guiname[20]
int	i, k
int	nbox
pointer	bdat
real	xics, yics
real	xb, yb
int	nxbox, nybox
real	xsz, ysz
int	n
int	nw
int	stat, chip
double	a3[3,3]
real	scaling
real	w0
pointer	im0, mmap

real	wave[3]
double	holdo1, holdo2, holdo3

int	qxfm(), loc_box()
begin
	iord = 1
	nv = NVAR
	niter = 12	# TMP HARDCODE

# TMP!! TEST:
	nw = 3
	wave[1] = 0.55773
#	wave[2] = 0.83446
	wave[2] = 8430.17e-4
	wave[3] = 8885.83e-4
	xsz = 34.
	ysz = 28.
	nxbox = 3. * xsz
	nybox = 3. * ysz

# locate alignment boxes and set up wavelengths:
	call get_boxes (fname, bdat, nbox, guiname, xics, yics)

	if (nbox == 0)
		call fatal (0, "No Alignment Boxes Found!")

	nmeas = nbox * nw

# Allocate arrays
	call malloc (bufxmm, nmeas, TY_REAL)
	call malloc (bufymm, nmeas, TY_REAL)
	call malloc (bufwav, nmeas, TY_REAL)
	call malloc (bufx,   nmeas, TY_REAL)
	call malloc (bufy,   nmeas, TY_REAL)
	call malloc (bufxp,  nmeas, TY_REAL)
	call malloc (bufyp,  nmeas, TY_REAL)

	call malloc (bufmat, nv*(nv+1), TY_DOUBLE)
	call malloc (bufpdx, (nv+1), TY_REAL)
	call malloc (bufpdy, (nv+1), TY_REAL)
	call malloc (bufdel, nv, TY_REAL)

# Open mosaic image:
	call mos_init (fname, DEF_PROP_FILE, im0, mmap, n)

# Fill the location vectors
	ndx = 0
	do i = 1, nw {
		w0 = wave[i]
		do k = 0, nbox-1 {
			call gsetup (a3, sys)
			scaling = 1.
			stat  = qxfm (map, a3, sys, ccd, XMM(bdat,k), YMM(bdat,k), w0, scaling, xics, yics, XPX(bdat,k), YPX(bdat,k), chip, YES, YES)
		
			if (XPX(bdat,k) == INDEF || YPX(bdat,k) == INDEF) {
				call eprintf ("Box off chip -- skipping\n")
				next
			}

			xb = XPX(bdat,k)
			yb = YPX(bdat,k)
			stat = loc_box (xb, yb, xsz, ysz, mmap, chip, nxbox, nybox)

			if (stat == OK) {
				Memr[bufxmm+ndx] = XMM(bdat,k)
				Memr[bufymm+ndx] = YMM(bdat,k)
				Memr[bufwav+ndx] = w0
				Memr[bufx+ndx] = xb	# XPX(bdat,k)
				Memr[bufy+ndx] = yb	# YPX(bdat,k)
				ndx = ndx + 1
call eprintf ("%7.3f %7.3f %7.5f  %6.1f %6.1f\n")
call pargr (Memr[bufxmm+ndx-1])
call pargr (Memr[bufymm+ndx-1])
call pargr (Memr[bufwav+ndx-1])
call pargr (Memr[bufx+ndx-1])
call pargr (Memr[bufy+ndx-1])
			}
		}
	}
	nmeas = ndx

	call mos_free (im0, mmap)

# ready to solve ...

# FIX!
	X_OPT(sys) =    0.	# WRONG!!
	Y_OPT(sys) = -280.	# WRONG!!
	CAM_FOC(sys) = 382.0d0 * 1.0001
	call eprintf ("Input: (xyz) %8.5f %7.4f %7.4f\n")
		call pargd (RADTODEG (MU(sys)))
		call pargd (RADTODEG (GR_YERR(sys)))
		call pargd (RADTODEG (GR_ZERR(sys)))
	call eprintf ("Input: scaling = %7.5f\n")
		call pargd (CAM_FOC(sys)/382.0d0)

	holdo1 = MU(sys)
	holdo2 = GR_YERR(sys)
	holdo3 = GR_ZERR(sys)

# Now solve for image geometry:
	call qlsolve (Memr[bufxmm], Memr[bufymm], Memr[bufwav], iord,
	Memr[bufx], Memr[bufy], Memr[bufxp], Memr[bufyp], nmeas, niter,
	ccd, sys, map, rmsx, rmsy, Memd[bufmat], Memr[bufpdx], Memr[bufpdy],
	Memr[bufdel], nv, nv+1)



	if (rmsx > 2. && rmsy > 2.) {		# TMP! HARDCODES
		call fatal (0, "Solution did not converge!!")
	} else if (rmsy > 3.) {
		call eprintf ("Solution: (xyz) INDEF %7.4f %7.4f (prob OK)\n")
			call pargd (RADTODEG (GR_YERR(sys)))
			call pargd (RADTODEG (GR_ZERR(sys)))
		call eprintf ("Adjust dmu,droll,do3 by: (lots) %7.4f %7.4f\n")
			call pargd (RADTODEG (GR_YERR(sys) - holdo2))
			call pargd (RADTODEG (GR_ZERR(sys) - holdo3))
		call fatal (0, "Tilt solution did not converge!!")
	} else {
		call eprintf ("Solution: (xyz) %8.5f %7.4f %7.4f\n")
			call pargd (RADTODEG (MU(sys)))
			call pargd (RADTODEG (GR_YERR(sys)))
			call pargd (RADTODEG (GR_ZERR(sys)))
		call eprintf ("Solution: scaling = %7.5f\n")
			call pargd (CAM_FOC(sys)/382.0d0)

		call eprintf ("Adjust dmu,droll,do3 by: %7.4f %7.4f %7.4f\n")
			call pargd (RADTODEG (MU(sys) - holdo1))
			call pargd (RADTODEG (GR_YERR(sys) - holdo2))
			call pargd (RADTODEG (GR_ZERR(sys) - holdo3))
	}

# mfree ...

call eprintf ("REMINDER: NEED TO FIX HARDCODES, DETERMINE RANGE ...\n")
end




procedure qlsolve (xs, ys, w, iord, xobs, yobs, x, y, nmeas, niter, ccd, sys, map, rmsx, rmsy, mat, px, py, del, nv, nc)

real	xs[nmeas], ys[nmeas]		# SM coords of points
real	w[nmeas]			# wavelength
int	iord				# order of spectrum
real	xobs[nmeas], yobs[nmeas]	# observed pixel locations
real	x[nmeas], y[nmeas]		# returned best-fit positions
int	nmeas				# number of slits
int	niter				# number of iterations
real	ccd[NCCD,3]			# CCD geometry
double	sys[NPARAM]			# system params
pointer	map[ARB]			# vector of map pointers
real	rmsx, rmsy			# returned rms values
#
double	mat[nc,nv]
real	px[nc], py[nc]			# partials
real	del[nv]				# corrections
int	nv, nc				# nvar, nvar+1

double	a3[3,3]			 	# transforms
int	n				# CCD number
int	stat

pointer	asi					# pointer to interp. fit

int	i, j, k, nn
real	delang, delfoc
real	delta
double	sumx, sumy, sumxx, sumyy
real	wtx, wty
real	xccd, yccd, xics, yics

pointer	sfx, sfy

# XXX TMP -- decide
real	scaling

int	qxfm()
begin

# Set limits for delta's --> also used for partials: vinc = 2.*delpar / (NORD-1)
delang = 0.006	# 0.007
delfoc = 0.30	# 0.5


	call asiinit (asi, II_SPLINE3)


# CCD: define the geometry of the mosaic
	call ccd_geom (ccd, sys, YES)

	call eprintf ("\n rmsx rmsy\n")

	do nn = 1, niter {
		call amovkd (double(0.), mat, nv*nc)
		rmsx = 0.
		rmsy = 0.
		sumxx = 0.d0
		sumyy = 0.d0
		sumx = 0.
		sumy = 0.


	    do k = 1, nmeas {

# Evaluate the function and residuals:

# XXX qxfm as below in qnd
scaling = CAM_FOC(sys) / 382.0d0
		call gsetup (a3, sys)
		stat  = qxfm (map, a3, sys, ccd, xs[k], ys[k], w[k], scaling, x[k], y[k], xccd, yccd, n, YES, YES)
		

# Convert observed to ICS -- this is inefficient; should be done once PROVIDED
# the detector geometry is not changing.
		call get_ccdnum (xobs[k], yobs[k], n, xccd, yccd)
		call ccd_to_ics (xccd, yccd, ccd, n, xics, yics)
		px[nc] = xics - x[k]
		py[nc] = yics - y[k]

#		px[nc] = xobs[k] - x[k]
#		py[nc] = yobs[k] - y[k]
		sumx = sumx + px[nc]
		sumy = sumy + py[nc]
		sumxx = sumxx + px[nc] * px[nc]
		sumyy = sumyy + py[nc] * py[nc]

# Solve partials (note some of these are done analytically):

		px[CF] = (x[k] - X_OPT(sys)) / CAM_FOC(sys)
		py[CF] = (y[k] - Y_OPT(sys)) / CAM_FOC(sys)

		call qnd (MU(sys),      delang, xs[k], ys[k], w[k], asi,
		px[GX], py[GX], ccd, sys, a3, map, NO)

		call qnd (GR_YERR(sys), delang, xs[k], ys[k], w[k], asi,
		px[GY], py[GY], ccd, sys, a3, map, NO)

		call qnd (GR_ZERR(sys), delang, xs[k], ys[k], w[k], asi,
		px[GZ], py[GZ], ccd, sys, a3, map, NO)

# Accumulate data in matrix:
		wtx = 1.
		wty = 1.

		do j = 1, nv {
		    do i = 1, nc {
			mat[i,j] = mat[i,j] + wtx*(px[i] * px[j]) + wty*(py[i] * py[j])
		    }
		}
	    }
# Remove FIXED variables from fit
#	    do j = 1, nv {
#		if (fix[j] == YES) {
#			mat[j, j] = mat[j, j] + 1.
#			mat[nc,j] = mat[nc,j] + 0.
#		}
#	    }

# Solve matrix
	    call g2_elim (mat, nv)

# Get deltas
		delta = mat[nc,GX]
		del[GX] = 0.79 * max (min (delta, delang), -delang)

		delta = mat[nc,GY]
		del[GY] = 0.79 * max (min (delta, delang), -delang)

		delta = mat[nc,GZ]
		del[GZ] = 0.79 * max (min (delta, delang), -delang)

		delta = mat[nc,CF]
		del[CF] = 0.79 * max (min (delta, delfoc), -delfoc)

# TMP!  This could be cleaner, avoiding any calc if held fixed.
# Set fixed variable to zero delta
#	    do j = 1, nv {
#		if (fix[j] == YES)
#			del[j] = 0.
#	    }

	    rmsx = sqrt (sumxx / nmeas)
	    rmsy = sqrt (sumyy / nmeas)

	    call eprintf ( "%2d %5f %5f\n")
		call pargi (nn)
		call pargr (rmsx)
		call pargr (rmsy)

if (niter == 1)		# assume update is not desired!
	break

# Update parameters

		CAM_FOC(sys) = CAM_FOC(sys) + del[CF]
		MU(sys)      = MU(sys)      + del[GX]
		GR_YERR(sys) = GR_YERR(sys) + del[GY]
		if (ORDER(sys) == 0.)
			del[GZ] = 0.
		GR_ZERR(sys) = GR_ZERR(sys) + del[GZ]
	}

# Free work vectors
call eprintf ("avg: %f %f\n")
call pargr (sumx / nmeas)
call pargr (sumy / nmeas)
	call gsfree (sfy)
	call gsfree (sfx)
	call asifree (asi)

end

#
# QND: find numerical derivative by varying setup slightly. Works to ICS. QXFM
# version
#
#	call qnd (GR_ZERR(sys), delang, xs[k], ys[k], w[k], asi,
#		px[GZ], py[GZ], ccd, sys, a3, map, NO)

procedure qnd (varpar, delpar, xs, ys, w, asi, px, py, ccd, sys, a3, map, fixed)

double	varpar				# system parameter to vary 
real	delpar				# amount to vary parameter
real	xs, ys				# x,y on slit (mm)
real	w				# ref wavelength
pointer	asi				# pointer to fit structure
real	px, py				# partials in x, y
real	ccd[NCCD,3]			# CCD geometry
double	sys[NPARAM]			# system parameters for pt_xfm
double	a3[3,3]				# grating transform
pointer	map[ARB]			# vector of mappings
int	fixed				# Is parameter fixed?

int	n			# CCD number

int	i, ndx
int	noff
double	par0				# original value (must be restored!)
double	vinc
real	x[NORD], y[NORD]		# work arrays

real	deriv[2]

# NEW
int	stat
real	xpix, ypix
real	scaling
int	qxfm()

begin
# Check to see if anything must be done:
	if (fixed == YES) {
		px = 0.
		py = 0.
		return
	}

# store parameter
	par0 = varpar

	vinc = 2. * delpar / (NORD - 1)
	noff = NORD / 2

	ndx = 1
	do i = -noff, noff {
		varpar = par0 + vinc * i
		call gsetup (a3, sys)
# check types (xics vs x); or do we want pixel space??
# FIX (cleanup)
scaling = CAM_FOC(sys) / 382.0d0
		stat  = qxfm (map, a3, sys, ccd, xs, ys, w, scaling, x[ndx], y[ndx], xpix, ypix, n, YES, NO)
		
#		call pt_xfm (xs, ys, w, e1, a2, as, a3, a4, ccd, sys, x[ndx], y[ndx], n)
		ndx = ndx + 1
	}

	call asifit (asi, x, NORD)
	call asider (asi, real (noff+1), deriv, 2)
	px = deriv[2] / vinc
	
	call asifit (asi, y, NORD)
	call asider (asi, real (noff+1), deriv, 2)
	py = deriv[2] / vinc

# Restore parameter:
	varpar = par0
end


#
# LOC_BOX: locate a box/feature (somewhat inefficient in re-allocating
# vectors repeatedly, but it remains general this way).
#

int	procedure loc_box (xb, yb, xsz, ysz, mmap, chip, nxbox, nybox)

real	xb, yb			# init/fin box center
real	xsz, ysz		# px size of box
pointer	mmap			# pointer to image mapping
int	chip			# chip for feature
int	nxbox, nybox		# extraction box size

int	i, j
int	x1, x2, y1, y2
int	nx, ny
int	stat
real	zmax, zmin
pointer	buf
pointer bufbx, bufby, bufzx, bufzy

int	box_center()
pointer	chip_sect()
begin

# Allocate arrays for marginal plots
	    call malloc (bufbx, nxbox, TY_REAL)
	    call malloc (bufby, nybox, TY_REAL)
	    call malloc (bufzx, nxbox, TY_REAL)
	    call malloc (bufzy, nybox, TY_REAL)

# Loop through expected box positions to measure actual box/star positions

# Find the Box; we do a kludge -- run this twice, recentered on box the 2nd time
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


# Get the image section (limit checking and zero-fill included)
		    buf = chip_sect (x1, x2, y1, y2, mmap, chip, YES)

# TMP? Check for common error of amp dropout
		    call alimr (Memr[buf], nx*ny, zmax, zmin)
		    if (zmax == zmin) {
			call eprintf ("Subraster contains no useful data!\n")
			stat = ERR
			break
		    }
#		    if (zmax > SATURLEV)
#			call eprintf ("WARNING! Star may be saturated!!\007\n")

# Fill position vectors
		    do i = 0, nx-1 {
			Memr[bufbx+i] = i + x1
		    }
		    do i = 0, ny-1 {
			Memr[bufby+i] = i + y1
		    }
# Get the box position
		    stat = box_center (Memr[bufbx], Memr[bufby], Memr[buf],
			nx, ny, Memr[bufzx], Memr[bufzy], xsz, ysz, xb, yb)

		    if (stat != OK)
			break
		}

# free box vectors
	call mfree (bufzy, TY_REAL)
	call mfree (bufzx, TY_REAL)
	call mfree (bufby, TY_REAL)
	call mfree (bufbx, TY_REAL)

	return (stat)
end


#		    call fprintf (fd,
#			"%6.1f %6.1f  %s   %1s %5.2f # %8.3f %7.3f\n")
#			call pargr (XPX(bdat,i))
#			call pargr (YPX(bdat,i))
#			call pargstr (OBJNAME(bdat,i))
#			call pargc (PBAND(bdat,i))
#			call pargr (MAG(bdat,i))
#			call pargr (XMM(bdat,i))
#			call pargr (YMM(bdat,i))



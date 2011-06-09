define	MAX_RAD   4950.  #Max radius in pixels to accept pairs
define PRECOL 0
# define	PRECOL	96
# define	PRECOL	48


# TBD: resolve ICS vs MOSAIC coords (search for FIX)
#  Note -- this comes up for camera focus; currently hardcoded -- fix
# Feb03: modify to put in order for each


####
# Mod's: change so that multiple grating/angle info can be entered:
# 1. Set number  must be input per line
# 1a. Along with set number must be specified: grating, initial mu, ybar, order
# 2. solution: for each SET: solve for mu, yaw, tip, focus
# 4. For SYSTEM: solve coll (2 values), xy_opt, camerr, C2 C4, ccd_rot [cam_ang]

# MSOLVE: solve multi sets of pt/wavelength-pixel pairs to produce system params
#
# Modified 2001sep06 to work in ICS rather than individual pixel space

## 2001sep07: it appears that the solution is quite sensitive to numderiv 
# returns, which vary with fitting function, order and algorithm for stepping.
# Presumably this is also true of the step size.

## 2002jan22: There's a huge problem with the solutions for CCD_geom angles ...
# apparently this has arisen from change of CDD to ICS -- OF COURSE! there is
# no change in ICS as CCD is rotated ...

include <imhdr.h>
include <math.h>
include <math/iminterp.h>
include "instrument.h"

define		SZ_GRNM	16	# size in char of grating name in header
define		ID_CHSZ	16		# Character size of ID string

define		NORD	7	# "order" of fit to derivatives

# Grating/focus params (always solved)
define	GY	3
define	GZ	4
define	GX	2
define	CF	1
define	GL	5
define	GO	6

# Detector params
define	DX	1
define	DY	2
define	DZ	3

# System params (optionally fit)
define	LZ	1
define	LT	2
define	LP	3

define	TA	4
define	TP	5

define	CA	6
define	CP	7

define	MX	8
define	MY	9
define	MP	10

define	SZ	11	# ADDED
define	SY	12	# ADDED

define	NOPTVAR	12		# 7 optical system
define	NCCDVAR	24
define	NSYSVAR	36		# NOPTVAR + NCCDVAR


procedure t_msolve ()

char	input[SZ_FNAME]			# slit info file
char	output[SZ_FNAME]		# output extraction info
char	setinfo[SZ_FNAME]		# gratings table
char	dmap[SZ_FNAME]			# opt. distortion mapping
real	pa				# PA in degrees
char	newfile[SZ_FNAME]		# opt. rewrite of input w/ resid
# char	zmapname[SZ_FNAME]		# zeta-map name
char	grating[SZ_GRNM]		# name of grating
int	niter				# number of interations

#
pointer	fda, fdb, fdd			# file descriptor for in/output files
pointer	fds				# file descriptor for data set

real	ccd[NCCD,3]			# CCD geometry
double	sys[NPARAM]			# system params

int	nslit				# No. slits
pointer	bufxmm				# pointer to slit x
pointer	bufymm				# pointer to slit y 
pointer	bufwav, buford			# pointer to wavelength, order
pointer	bufx, bufy			# pointers to x,y (pixel) meas'd values
pointer	bufxad, bufyad			# pointers to x,y (pixel) adjustments
pointer	bufxp, bufyp			# pointers to x,y (pixel) predicted
pointer	bufset, bufpar			# pointers to dataset no., params
int	nset, nv
pointer	bufmat, bufpdx, bufpdy		# pointers to matrix array, partials
pointer	bufdel				# pointer to del-corrections
pointer	buffix				# pointer to fix/(vary) flag
real	grangle, roll, yaw, lines
int	order
int	npairs

char	tchar
int	ndx
int	i, j, m
int	off
int	iset

bool	rewrite
real	maxdev
real	dx, dy, dxmax, dymax, ds, dsmax
pointer	fdc

# real	relscale				# relative scale of 2 vectors
real	rmsx, rmsy

pointer	bufsq, bufsn

bool	strne()
int	access() 
int	fscan(), nscan()
int	clgeti()
real	clgetr()
pointer	open()

begin

# Read in parameters:
	call clgstr ("input", input, SZ_FNAME)
	fda = open (input, READ_ONLY, TEXT_FILE)

	call clgstr ("setinfo", setinfo, SZ_FNAME)
	fds = open (setinfo, READ_ONLY, TEXT_FILE)

# check for rewritten input file
	call clgstr ("new_input", newfile, SZ_FNAME)
	rewrite = strne (newfile, "")
	if (rewrite) {
	    if (access (newfile, 0, 0) == YES) {
		call eprintf ("newfile file %s exists!\n")
			call pargstr (newfile, SZ_FNAME)
		call fatal (0, "")
	    } else {
		fdc = open (newfile, NEW_FILE, TEXT_FILE)
	    }
	    maxdev = clgetr ("max_dev")
	}

#	call clgstr ("inst_config.sm2zeta", zmapname, SZ_FNAME)
#	if (access (zmapname, READ_ONLY, TEXT_FILE) == YES) {
#		fdz = open (zmapname, READ_ONLY, TEXT_FILE)
#	} else {
#		call fatal (0, "Cannot open zeta-map file! \n")
#	}

#	Y_BIN(sys) = 1.			# TMP !!

	COL_DST(sys) = clgetr ("trace.coll_zdst")
	COL_ERR(sys) = DEGTORAD(clgetr ("trace.coll_angle"))
	COL_PHI(sys) = DEGTORAD(clgetr ("trace.coll_phi"))

TNT_ANG(sys) = DEGTORAD(71.5 + clgetr ("trace.t2"))
TNT_PHI(sys) = DEGTORAD(90. + clgetr ("trace.p2"))

SL_ZERR(sys) = DEGTORAD(clgetr ("trace.os3"))
SL_YERR(sys) = DEGTORAD(clgetr ("trace.os2"))

	CAM_ANG(sys) = DEGTORAD(clgetr ("trace.cam_angle"))
	CAM_PHI(sys) = DEGTORAD(clgetr ("trace.cam_phi"))

	CAM_FOC(sys) = clgetr ("trace.cam_foc")

	X_OPT(sys) = clgetr ("trace.x_optaxis")
	Y_OPT(sys) = clgetr ("trace.y_optaxis")
	MOS_ROT(sys) = DEGTORAD(clgetr ("trace.mos_rotation"))

do i = 1, NCCD {
	CN_XERR(sys,i) = 0.
	CN_YERR(sys,i) = 0.
	CN_RERR(sys,i) = 0.
}

	niter = clgeti ("niter")

# Ready to deal with the data sets now:
# count the entries in the table table
	ndx = 0
	while (fscan (fds) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#') {
			next
		}
		ndx = ndx + 1
	}
	call seek (fds, BOF)

	call malloc (bufpar, ndx*6, TY_DOUBLE)

# count the entries in the table table
	ndx = 0
	nset = 0
	while (fscan (fds) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#') {
			next
		}
		call reset_scan()

		call gargi (iset)
		if (iset != nset+1)
			call eprintf ("Alignment error!!\n")

		call gargwrd (grating, SZ_GRNM)
		call gargr (grangle)
		call gargr (roll  )
		call gargr (yaw   )
		call gargr (lines )
		call gargi (order )
		call gargi (npairs)
		if (nscan() < 8)
			next
		if (npairs <= 0) {
			nset = nset + 1
			next
		}
	
		off = nset * 6 - 1
		Memd[bufpar+GX+off] = DEGTORAD(grangle)
		Memd[bufpar+GY+off] = DEGTORAD(roll)
		Memd[bufpar+GZ+off] = DEGTORAD(yaw)
		Memd[bufpar+GO+off] = order
		Memd[bufpar+GL+off] = lines * 1.e-3
		Memd[bufpar+CF+off] = CAM_FOC(sys)
		nset = nset + 1
		ndx = ndx + npairs
call eprintf ("Set %2d: %16s (%4.0f/mm, order=%02d) %6.3f %6.3f %6.3f %5d %4d\n")
call pargi (nset)
call pargstr (grating)
call pargr (lines)
call pargi (order)
call pargr (grangle)
call pargr (roll)
call pargr (yaw)
call pargi (npairs)
call pargi (ndx)
	}


# Allocate arrays
	call malloc (bufxmm, ndx, TY_DOUBLE)
	call malloc (bufymm, ndx, TY_DOUBLE)
	call malloc (bufwav, ndx, TY_DOUBLE)
	call malloc (buford, ndx, TY_INT)
	call malloc (bufx,   ndx, TY_REAL)
	call malloc (bufy,   ndx, TY_REAL)
	call malloc (bufxad, ndx, TY_REAL)
	call malloc (bufyad, ndx, TY_REAL)
	call malloc (bufxp,  ndx, TY_REAL)
	call malloc (bufyp,  ndx, TY_REAL)
	call malloc (bufset, ndx, TY_INT)	# which data "set"

# Get the features:
	ndx = 0
	while (fscan (fda) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#')
			next
		call reset_scan()

		call gargi (Memi[bufset+ndx])
		call gargd (Memd[bufxmm+ndx])
		call gargd (Memd[bufymm+ndx])
		call gargd (Memd[bufwav+ndx])
		call gargi (Memi[buford+ndx])
		call gargr (Memr[bufx  +ndx])
		call gargr (Memr[bufy  +ndx])

		if (nscan() < 7)
			next
		if (Memr[bufx+ndx] == INDEF || Memr[bufy+ndx] == INDEF)
			next
# TMP: test for radius (CRUDE):
		if (sqrt ((Memr[bufx+ndx]-4096.)**2+
				(Memr[bufy+ndx]-4096.)**2) > MAX_RAD) {
			call eprintf ("pair too extreme (set %d, w=%7f \n")
				call pargi (Memi[bufset+ndx])
				call pargd (Memd[bufwav+ndx])
			next
		}
		
		ndx = ndx + 1
	}
	call close (fdb)
	nslit = ndx

call eprintf ("CAVE!! TEST!! SCALING to x,ymm\n")
call amulkd (Memd[bufymm], (1.0d0 - 2.0d-5 * 17.), Memd[bufymm], nslit)
call amulkd (Memd[bufxmm], (1.0d0 - 2.0d-5 * 17.), Memd[bufxmm], nslit)

call eprintf ("Number of pairs: %d\n")
call pargi (nslit)



# Allocate arrays for solution:
	nv = NSYSVAR + 4 * nset
	call malloc (bufmat, nv*(nv+1), TY_DOUBLE)
	call malloc (bufpdx, (nv+1), TY_REAL)
	call malloc (bufpdy, (nv+1), TY_REAL)
	call malloc (bufdel, nv, TY_REAL)
	call malloc (buffix, nv, TY_INT)
	call amovki (NO, Memi[buffix], nv)

# TMP! FIX!
	call amulkd (Memd[bufwav], 1.0D-4, Memd[bufwav], nslit)

#########################
# NEW! Adjustments for distortion (and CCD motion?) -- PA dependent only
# However, we cheat somewhat, using the _observed_ location to give approx ICS
# position, which should be good enough for adjustments)

	call clgstr ("dmap", dmap, SZ_FNAME)
	if (strne (dmap, "")) {
		fdd = open (dmap, READ_ONLY, TEXT_FILE)
		pa = clgetr ("pa")

# ccd: define the geometry of the mosaic (this has not been done yet, so ...
	call ccd_geom (ccd, sys)

		call pa_adj (fdd, pa, Memr[bufx], Memr[bufy], Memr[bufxad],
						Memr[bufyad], nslit, ccd)
		call close (fdd)
	} else {
		call amovkr (0., Memr[bufxad], nslit)
		call amovkr (0., Memr[bufyad], nslit)
	}
#########################


# Now solve for image geometry:
	call multi_solve (Memd[bufxmm], Memd[bufymm], Memd[bufwav], Memi[buford],
	Memr[bufx], Memr[bufy], Memr[bufxad], Memr[bufyad], Memr[bufxp],
	Memr[bufyp], nslit, niter, ccd, sys, rmsx, rmsy, nset,
	Memd[bufpar], Memi[bufset], Memd[bufmat], Memr[bufpdx], Memr[bufpdy],
	Memr[bufdel], Memi[buffix], nv, nv+1)

# Calculate stats per set
# NB -- forced "set number" here too: FIX!?
	call calloc (bufsn, nset, TY_REAL)
	call calloc (bufsq, nset, TY_REAL)
	do i = 0, nslit-1 {
		j = Memi[bufset+i] - 1
		dx = Memr[bufxp+i] - Memr[bufx+i]
		dy = Memr[bufyp+i] - Memr[bufy+i]
		Memr[bufsq+j] = Memr[bufsq+j] + dx*dx + dy*dy
		Memr[bufsn+j] = Memr[bufsn+j] + 1.
	}
	do j = 0, nset-1 {
		call printf ("Set %2d: %6.2f px rms (n=%.0f)\n")
			call pargi (j+1)
			call pargr (sqrt (Memr[bufsq+j] / Memr[bufsn+j]))
			call pargr (Memr[bufsn+j])
	}
	call mfree (bufsn, TY_REAL)
	call mfree (bufsq, TY_REAL)

# Print out the input file w/ resids, or just the max resids:
	dxmax = 0.
	dymax = 0.
	dsmax = 0.
	if (rewrite) {
	    do i = 0, nslit-1 {
		dx = Memr[bufxp+i] - Memr[bufx+i]
		dy = Memr[bufyp+i] - Memr[bufy+i]
		ds = sqrt (dx*dx + dy*dy)
		dxmax = max (dxmax, abs (dx))
		dymax = max (dymax, abs (dy))
		dsmax = max (dsmax, ds)

		call fprintf (fdc,
		 "%2s %03d %8.3f %8.3f %9.2f %2d  %7.2f %7.2f  # %4.1f %4.1f\n")
		if (ds > maxdev)
			call pargstr ("##")
		else
			call pargstr ("  ")
		call pargi (Memr[bufset+i])
		call pargd (Memd[bufxmm+i])
		call pargd (Memd[bufymm+i])
		call pargd (Memd[bufwav+i]*1.D4)
		call pargi (Memi[buford+i])
		call pargr (Memr[bufx+i])
		call pargr (Memr[bufy+i])
		call pargr (dx)
		call pargr (dy)
	    }
	    call close (fdc)
	} else {
	    do i = 0, nslit-1 {
		dx = Memr[bufxp+i] - Memr[bufx+i]
		dy = Memr[bufyp+i] - Memr[bufy+i]
		ds = sqrt (dx*dx + dy*dy)
		dxmax = max (dxmax, abs (dx))
		dymax = max (dymax, abs (dy))
		dsmax = max (dsmax, ds)
	    }
	}

	call eprintf (" Max differences %5.2fx %5.2fy  %5.2ftot\n")
		call pargr (dxmax)
		call pargr (dymax)
		call pargr (dsmax)


# All OK, so open the output file
	call clgstr ("output", output, SZ_FNAME)
	if (access (output, 0, 0) == YES) {
		call eprintf ("output file %s exists!\n")
			call pargstr (output, SZ_FNAME)
		call fatal (0, "")
	}
	fdd = open (output, NEW_FILE, TEXT_FILE)


# Write output file:
#	call fprintf (fdd, "# Zeta Map:       %s\n")
#		call pargstr (zmapname)

	call fprintf (fdd, "# Info file: %s\n")
		call pargstr (setinfo)

	call fprintf (fdd, "# (Solved) parameters:\n")

	call fprintf (fdd, "#  COL_DST =%10.2f\t%s\n")
		call pargd (COL_DST(sys))
	    if (Memi[buffix-1+LZ] == NO)
		call pargstr ("(solved)")
	    else
		call pargstr ("(fixed) ")

	call fprintf (fdd, "#  COL_ERR =%10.5f\t%s\n")
		call pargd (RADTODEG (COL_ERR(sys)))
	    if (Memi[buffix-1+LT] == NO)
		call pargstr ("(solved)")
	    else
		call pargstr ("(fixed) ")

	call fprintf (fdd, "#  COL_PHI =%10.5f\t%s\n")
		call pargd (RADTODEG (COL_PHI(sys)))
	    if (Memi[buffix-1+LP] == NO)
		call pargstr ("(solved)")
	    else
		call pargstr ("(fixed) ")

	call fprintf (fdd, "#  TNT_ANG =%10.5f\t%s\n")
		call pargd (RADTODEG (TNT_ANG(sys)))
	    if (Memi[buffix-1+TA] == NO)
		call pargstr ("(solved)")
	    else
		call pargstr ("(fixed) ")

	call fprintf (fdd, "#  TNT_PHI =%10.5f\t%s\n")
		call pargd (RADTODEG (TNT_PHI(sys)))
	    if (Memi[buffix-1+TP] == NO)
		call pargstr ("(solved)")
	    else
		call pargstr ("(fixed) ")

	call fprintf (fdd, "#  CAM_ANG =%10.5f\t%s\n")
		call pargd (RADTODEG (CAM_ANG(sys)))
	    if (Memi[buffix-1+CA] == NO)
		call pargstr ("(solved)")
	    else
		call pargstr ("(fixed) ")

	call fprintf (fdd, "#  CAM_PHI =%10.5f\t%s\n")
		call pargd (RADTODEG (CAM_PHI(sys)))
	    if (Memi[buffix-1+CP] == NO)
		call pargstr ("(solved)")
	    else
		call pargstr ("(fixed) ")

	call fprintf (fdd, "#  X_OPT =%7.2f\t%s\n")
		call pargd (X_OPT(sys))
	    if (Memi[buffix-1+MX] == NO)
		call pargstr ("(solved)")
	    else
		call pargstr ("(fixed) ")

	call fprintf (fdd, "#  Y_OPT =%7.2f\t%s\n")
		call pargd (Y_OPT(sys))
	    if (Memi[buffix-1+MY] == NO)
		call pargstr ("(solved)")
	    else
		call pargstr ("(fixed) ")

	call fprintf (fdd, "#  MOS_ROT =%10.5f\t%s\n")
		call pargd (RADTODEG (MOS_ROT(sys)))
	    if (Memi[buffix-1+MP] == NO)
		call pargstr ("(solved)")
	    else
		call pargstr ("(fixed) ")

	call fprintf (fdd, "#  CCDGEOM %s:  xpix  ypix    zdeg\n")
	    if (Memi[buffix-1+MP] == NO)
		call pargstr ("(solved)")
	    else
		call pargstr ("(fixed) ")
	do m = 1, NCCD {
		call fprintf (fdd, "#  %2d: %7.2f  %7.2f  %7.4f\n")
			call pargi (m)
			call pargd (CN_XERR(sys,m))
			call pargd (CN_YERR(sys,m))
			call pargd (RADTODEG(CN_RERR(sys,m)))
	}

	do m = 1, nset {
		off = (m - 1) * 6 - 1
		MU(sys)      = Memd[bufpar+GX+off]
		GR_YERR(sys) = Memd[bufpar+GY+off]
		GR_ZERR(sys) = Memd[bufpar+GZ+off]
		CAM_FOC(sys) = Memd[bufpar+CF+off]

		call fprintf (fdd,
		"# %2d: MU=%9.5f;   ROLL=%8.5f;    YAW=%8.5f;    FOC=%8.3f\n")
			call pargi (m)
			call pargd (RADTODEG (MU(sys)))
			call pargd (RADTODEG (GR_YERR(sys)))
			call pargd (RADTODEG (GR_ZERR(sys)))
			call pargd (CAM_FOC(sys))
	}

	call fprintf (fdd, "#   Slider zerr= %7.4f; slider y-err= %7.4f \n")
		call pargd (RADTODEG (SL_ZERR(sys)))
		call pargd (RADTODEG (SL_YERR(sys)))

	call fprintf (fdd, "#   rmsx=%5.3f    rmsy=%5.3f   npair=%d\n")
		call pargr (rmsx)
		call pargr (rmsy)
		call pargi (nslit)


	call fprintf (fdd, "#\n# xmeas ymeas   xpred ypred    xmm  ymm  wave\n")

	do i = 0, nslit-1 {
		call fprintf (fdd, "%7.2f %7.2f    %7.2f %7.2f   #  %7.3f %7.3f %9.5f %2d  %2d\n")
		call pargr (Memr[bufx+i])
		call pargr (Memr[bufy+i])
		call pargr (Memr[bufxp+i])
		call pargr (Memr[bufyp+i])
		call pargd (Memd[bufxmm+i])
		call pargd (Memd[bufymm+i])
		call pargd (Memd[bufwav+i])
		call pargi (Memi[buford+i])
		call pargi (Memr[bufset+i])
	}
	call close (fdd)


# Release memory
	call mfree (buffix, TY_INT)
	call mfree (bufdel, TY_REAL)
	call mfree (bufpdx, TY_REAL)
	call mfree (bufpdy, TY_REAL)
	call mfree (bufmat, TY_DOUBLE)
	call mfree (bufset, TY_INT)
	call mfree (bufpar, TY_DOUBLE)
	call mfree (bufyp, TY_REAL)
	call mfree (bufxp, TY_REAL)
	call mfree (bufyad, TY_REAL)
	call mfree (bufxad, TY_REAL)
	call mfree (bufy, TY_REAL)
	call mfree (bufx, TY_REAL)
	call mfree (buford, TY_INT)
	call mfree (bufwav, TY_DOUBLE)
	call mfree (bufymm, TY_DOUBLE)
	call mfree (bufxmm, TY_DOUBLE)
end


#########################################################
## SUBROUTINES FOLLOW:
################################################

#
# MULTI_SOLVE: accumlate and perform fit
#

procedure multi_solve (xs, ys, w, iord, xobs, yobs, xadj, yadj, x, y, nslit, niter, ccd, sys, rmsx, rmsy, nset, spar, setndx, mat, px, py, del, fix, nv, nc)

double	xs[nslit], ys[nslit]		# SM coords of points
double	w[nslit]			# wavelength
int	iord[nslit]			# order of spectrum
real	xobs[nslit], yobs[nslit]	# observed pixel locations
real	xadj[nslit], yadj[nslit]	# PA-dependent distortion corrections
real	x[nslit], y[nslit]		# returned best-fit positions
int	nslit				# number of slits
int	niter				# number of iterations
real	ccd[NCCD,3]			# CCD geometry
double	sys[NPARAM]			# system params
real	rmsx, rmsy			# returned rms values
#
int	nset				# number of data sets
double	spar[6,nset]			# "set" data
int	setndx[nslit]
double	mat[nc,nv]
real	px[nc], py[nc]			# partials
real	del[nv]				# corrections
int	fix[nv]				# hold variable fixed?
int	nv, nc				# nvar, nvar+1

double	e1[3,3], a2[3,3], a3[3,3], a4[3,3]	# transforms
double	as[3,3]				# slider transform
int	n				# CCD number
int	stat

pointer	asi					# pointer to interp. fit

int	i, j, k, nn, m
int	ndx
real	delang, delphi, delpix, delfoc
real	delta
double	sumx, sumy, sumxx, sumyy
real	wtx, wty
real	xccd, yccd, xics, yics


pointer	sfx, sfy

bool	solve_geom		# solve mosaic geometry?

int	ident_ccd()
bool	clgetb()
int	btoi()
begin

# Set limits for delta's --> also used for partials: vinc = 2.*delpar / (NORD-1)
delphi = 0.04	# This is for collimator err phi
delang = 0.006	# 0.007
delfoc = 0.30	# 0.5
delpix = 30		# 50.


## START TO FIX THE VARIABLES!

fix[LZ] = btoi (!(clgetb ("vlz")))
fix[LT] = btoi (!(clgetb ("vlt")))
fix[LP] = btoi (!(clgetb ("vlp")))
fix[TA] = btoi (!(clgetb ("vta")))
fix[TP] = btoi (!(clgetb ("vtp")))
fix[CA] = btoi (!(clgetb ("vca")))
fix[CP] = btoi (!(clgetb ("vcp")))
fix[MX] = btoi (!(clgetb ("vmx")))
fix[MY] = btoi (!(clgetb ("vmy")))
fix[MP] = btoi (!(clgetb ("vmp")))
# NB: these should remain fixed until the roll/yaw for gratings held constant ..
fix[SZ] = YES
fix[SY] = YES

solve_geom = clgetb ("vdg")

if (solve_geom) {
	do m = NOPTVAR+1, NSYSVAR
		fix[m] = NO
# fix[MX] = YES
# fix[MY] = YES
# fix[MP] = YES
} else {
	do m = NOPTVAR+1, NSYSVAR
		fix[m] = YES
}


# If the order is 0, fix the YAW variable
	do m = 1, nset {
# call eprintf ("SPECIAL!  FIX YAW !! \n")
# 			ndx = (m-1) * 4 + NSYSVAR
#  			fix[GZ+ndx] = YES
		if (spar[GO,m] == 0.) {
			ndx = (m-1) * 4 + NSYSVAR
			fix[GZ+ndx] = YES
		}
	}


	call asiinit (asi, II_SPLINE3)


	call eprintf ("\n COL_ERR COL_PHI TN_ANG TN_PHI X_OPT  Y_OPT  ROT COL_DST rmsx rmsy\n")

	do nn = 1, niter {
	    call amovkd (double(0.), mat, nv*nc)
	    rmsx = 0.
	    rmsy = 0.
		sumxx = 0.d0
		sumyy = 0.d0
		sumx = 0.
		sumy = 0.

	    do k = 1, nslit {

# Assign values appropriate to set:
		m = setndx[k]
		MU(sys) = spar[GX,m]
		GR_YERR(sys) = spar[GY,m]
		GR_ZERR(sys) = spar[GZ,m]
		CAM_FOC(sys) = spar[CF,m]
		ORDER(sys) = iord[k]		# spar[GO,m]
		GRLINES(sys) = spar[GL,m]

# Evaluate the function and residuals:
		call setup (e1, a2, as, a3, a4, ccd, sys)
		call pt_xfm (xs[k], ys[k], w[k], e1, a2, as, a3, a4, ccd, sys, x[k], y[k], n)

# Add distortion predictions here:
		x[k] = x[k] + xadj[k]
		y[k] = y[k] + yadj[k]


# Convert observed to ICS -- this is inefficient; should be done once PROVIDED
# the detector geometry is not changing.
		call get_ccdnum (xobs[k], yobs[k], n, xccd, yccd)
		call ccd_to_ics (xccd, yccd, ccd, n, xics, yics)
		px[nc] = xics - x[k]
		py[nc] = yics - y[k]

# KLUDGE: this is a kludge that has to do with wild points -- bad!!
if (abs (px[nc]) > 90.) {
call eprintf ("DEBUG, px=%5.1f (line=%d, chip=%d; w=%5f, xyobs= %4.0f,%4.0f)\n")
call pargr (px[nc])
call pargi (k)
call pargi (n)
call pargd (w[k])
call pargr (xobs[k])
call pargr (yobs[k])
next
}

#		px[nc] = xobs[k] - x[k]
#		py[nc] = yobs[k] - y[k]
		sumx = sumx + px[nc]
		sumy = sumy + py[nc]
		sumxx = sumxx + px[nc] * px[nc]
		sumyy = sumyy + py[nc] * py[nc]

# Solve partials (note some of these are done analytically):
		call numderiv (COL_DST(sys), delang, xs[k], ys[k], w[k], asi,
			px[LZ], py[LZ], ccd, sys, e1, a2, as, a3, a4, fix[LZ])

		call numderiv (COL_ERR(sys), delang, xs[k], ys[k], w[k], asi,
			px[LT], py[LT], ccd, sys, e1, a2, as, a3, a4, fix[LT])

		call numderiv (COL_PHI(sys), delphi, xs[k], ys[k], w[k], asi,
			px[LP], py[LP], ccd, sys, e1, a2, as, a3, a4, fix[LP])

		call numderiv (TNT_ANG(sys), delang, xs[k], ys[k], w[k], asi,
			px[TA], py[TA], ccd, sys, e1, a2, as, a3, a4, fix[TA])

		call numderiv (TNT_PHI(sys), delang, xs[k], ys[k], w[k], asi,
			px[TP], py[TP], ccd, sys, e1, a2, as, a3, a4, fix[TP])

		call numderiv (CAM_ANG(sys), delang, xs[k], ys[k], w[k], asi,
			px[CA], py[CA], ccd, sys, e1, a2, as, a3, a4, fix[CA])

		call numderiv (CAM_PHI(sys), delang, xs[k], ys[k], w[k], asi,
			px[CP], py[CP], ccd, sys, e1, a2, as, a3, a4, fix[CP])

# These are technically not true, but close enough:
		if (fix[MX] == NO)		# X_OPT
			px[MX] = 1.
		else
			px[MX] = 0.
		py[MX] = 0.

		if (fix[MY] == NO)		# Y_OPT
			py[MY] = 1.
		else
			py[MY] = 0.
		px[MY] = 0.

		call numderiv (MOS_ROT(sys), delang, xs[k], ys[k], w[k], asi,
			px[MP], py[MP], ccd, sys, e1, a2, as, a3, a4, fix[MP])

# ADDED -- change to SL_xERR?
		call numderiv (SL_ZERR(sys), delang, xs[k], ys[k], w[k], asi,
			px[SZ], py[SZ], ccd, sys, e1, a2, as, a3, a4, fix[SZ])

		call numderiv (SL_YERR(sys), delang, xs[k], ys[k], w[k], asi,
			px[SY], py[SY], ccd, sys, e1, a2, as, a3, a4, fix[SY])


# Now for the CCD params (all zero except for chip n):
		call amovkr (0., px[NOPTVAR+1], NCCDVAR)
		call amovkr (0., py[NOPTVAR+1], NCCDVAR)
		if (solve_geom && n != 3) {
		    ndx = (n-1)*3 + NOPTVAR
		    px[DX+ndx] = -1.
		    py[DX+ndx] = 0.
		    px[DY+ndx] = 0.
		    py[DY+ndx] = -1.
		    call numder_ccd (CN_RERR(sys,n), delang, xs[k], ys[k], w[k],
	asi, px[DZ+ndx], py[DZ+ndx], ccd, sys, e1, a2, as, a3, a4, fix[DZ+ndx])
		}
			
			
# Below are set-specific; most are zero ...
		call amovkr (0., px[NSYSVAR+1], nset*4)
		call amovkr (0., py[NSYSVAR+1], nset*4)

# ... but some are not
		ndx = (m-1) * 4 + NSYSVAR

# FIX FIX!!!  TMP HARDCODES!
## FIX fix[CF] never defined...
fix[CF+ndx] = NO
		if (fix[CF+ndx] == NO) {	#CAM_FOC
			px[CF+ndx] = (x[k] - X_OPT(sys)) / CAM_FOC(sys)
			py[CF+ndx] = (y[k] - Y_OPT(sys)) / CAM_FOC(sys)
# 			px[CF+ndx] = (x[k] - 4144 - X_OPT(sys)) / CAM_FOC(sys)
# 			py[CF+ndx] = (y[k] - 4096 - Y_OPT(sys)) / CAM_FOC(sys)
		} else {
			px[CF+ndx] = 0.
			py[CF+ndx] = 0.
		}

		call numderiv (MU(sys), delang, xs[k], ys[k], w[k], asi,
		px[GX+ndx], py[GX+ndx], ccd, sys, e1, a2, as, a3, a4, fix[GX+ndx])

		call numderiv (GR_YERR(sys), delang, xs[k], ys[k], w[k], asi,
		px[GY+ndx], py[GY+ndx], ccd, sys, e1, a2, as, a3, a4, fix[GY+ndx])

		call numderiv (GR_ZERR(sys), delang, xs[k], ys[k], w[k], asi,
		px[GZ+ndx], py[GZ+ndx], ccd, sys, e1, a2, as, a3, a4, fix[GZ+ndx])

### FIX!  The following seems to be some weighting scheme
# Accumulate data in matrix:
		wtx = 1. #+ abs (x[k] - X_OPT(sys)) / 512.
		wty = 1. #+ exp (-1.* (x[k] - X_OPT(sys))**2 / 512.) + abs (y[k] - Y_OPT(sys)) / 1024.
		wtx = wtx #* 4.

		do j = 1, nv {
		    do i = 1, nc {
#call eprintf ("%d:px[%d]=%f, py[%d]=%f\n")
#call pargi (j)
#call pargi (i)
#call pargr (px[i])
#call pargi (i)
#call pargr (py[i])
			mat[i,j] = mat[i,j] + wtx*(px[i] * px[j]) + wty*(py[i] * py[j])
		    }
		}
	    }
# Remove FIXED variables from fit
	    do j = 1, nv {
		if (fix[j] == YES) {
			mat[j, j] = mat[j, j] + 1.
			mat[nc,j] = mat[nc,j] + 0.
		}
	    }
	    ndx = (3-1)*3 + NOPTVAR
	    j = DX+ndx
	    mat[j, j] = mat[j, j] + 1.
	    mat[nc,j] = mat[nc,j] + 0.
	    j = DY+ndx
	    mat[j, j] = mat[j, j] + 1.
	    mat[nc,j] = mat[nc,j] + 0.
	    j = DZ+ndx
	    mat[j, j] = mat[j, j] + 1.
	    mat[nc,j] = mat[nc,j] + 0.

# Solve matrix (getting _BIG_!)
	    call g2_elim (mat, nv)


## TMP -- 0.9/0.6 is cooling
	    delta = mat[nc,LZ]
	    del[LZ] = 0.79 * max (min (delta, delfoc), -delfoc)

	    delta = mat[nc,LT]
	    del[LT] = 0.79 * max (min (delta, delang), -delang)

	    delta = mat[nc,LP]
	    del[LP] = 0.79 * max (min (delta, delphi), -delphi)

	    delta = mat[nc,TA]
	    del[TA] = 0.79 * max (min (delta, delang), -delang)

	    delta = mat[nc,TP]
	    del[TP] = 0.79 * max (min (delta, delang), -delang)

	    delta = mat[nc,CA]
	    del[CA] = 0.79 * max (min (delta, delang), -delang)

	    delta = mat[nc,CP]
	    del[CP] = 0.79 * max (min (delta, delang), -delang)

	    delta = mat[nc,MX]
	    del[MX] = 0.79 * max (min (delta, delpix), -delpix)

	    delta = mat[nc,MY]
	    del[MY] = 0.79 * max (min (delta, delpix), -delpix)

	    delta = mat[nc,MP]
	    del[MP] = 0.79 * max (min (delta, delang), -delang)

# ADDED
	    delta = mat[nc,SZ]
	    del[SZ] = 0.79 * max (min (delta, delang), -delang)

	    delta = mat[nc,SY]
	    del[SY] = 0.79 * max (min (delta, delang), -delang)

# CCD specific
	    if (solve_geom) {
		do m = 1, NCCD {
		    ndx = (m-1)*3 + NOPTVAR
		    delta = mat[nc,DX+ndx]
		    del[DX+ndx] = 0.79 * max (min (delta, delpix), -delpix)
		    delta = mat[nc,DY+ndx]
		    del[DY+ndx] = 0.79 * max (min (delta, delpix), -delpix)
		    delta = mat[nc,DZ+ndx]
		    del[DZ+ndx] = 0.79 * max (min (delta, delang), -delang)
		}
	    } else {
		do m = 1, NCCD {
		    ndx = (n-1)*3 + NOPTVAR
		    del[DX+ndx] = 0.
		    del[DY+ndx] = 0.
		    del[DZ+ndx] = 0.
		}
	    } 

# Set-specific
	    do m = 1, nset {
		ndx = (m-1) * 4 + NSYSVAR

		delta = mat[nc,GX+ndx]
		del[GX+ndx] = 0.79 * max (min (delta, delang), -delang)

		delta = mat[nc,GY+ndx]
		del[GY+ndx] = 0.79 * max (min (delta, delang), -delang)

		delta = mat[nc,GZ+ndx]
		del[GZ+ndx] = 0.79 * max (min (delta, delang), -delang)

		delta = mat[nc,CF+ndx]
		del[CF+ndx] = 0.79 * max (min (delta, delfoc), -delfoc)
# call eprintf ("DEBUG!:  %d  %6f %6f %6f %6f\n")
# call pargi (ndx)
# call pargr (del[GX+ndx])
# call pargr (del[GY+ndx])
# call pargr (del[GZ+ndx])
# call pargr (del[CF+ndx])
	    }

# TMP!  This could be cleaner, avoiding any calc if held fixed.
# Set fixed variable to zero delta
	    do j = 1, nv {
		if (fix[j] == YES)
			del[j] = 0.
	    }

	    rmsx = sqrt (sumxx / nslit)
	    rmsy = sqrt (sumyy / nslit)

	    call eprintf ( "%2d %6.3f %6.3f %6.3f %6.3f %5.0f %5.0f %6.3f %6f %5f %5f\n")
		call pargi (nn)
		call pargd (RADTODEG (COL_ERR(sys)))
		call pargd (RADTODEG (COL_PHI(sys)))
		call pargd (RADTODEG (TNT_ANG(sys)))
		call pargd (RADTODEG (TNT_PHI(sys)))
		call pargd (X_OPT(sys))
		call pargd (Y_OPT(sys))
		call pargd (RADTODEG (MOS_ROT(sys)))
		call pargd (COL_DST(sys))
		call pargr (rmsx)
		call pargr (rmsy)

if (niter == 1)		# assume update is not desired!
	break

# Update parameters
	    COL_DST(sys) = COL_DST(sys) + del[LZ]
	    COL_ERR(sys) = COL_ERR(sys) + del[LT]
	    COL_PHI(sys) = COL_PHI(sys) + del[LP]
	    TNT_ANG(sys) = TNT_ANG(sys) + del[TA]
	    TNT_PHI(sys) = TNT_PHI(sys) + del[TP]
	    SL_ZERR(sys) = SL_ZERR(sys) + del[SZ]
	    SL_YERR(sys) = SL_YERR(sys) + del[SY]
	    CAM_ANG(sys) = CAM_ANG(sys) + del[CA]
	    CAM_PHI(sys) = CAM_PHI(sys) + del[CP]
	    X_OPT(sys)   = X_OPT(sys)   + del[MX]
	    Y_OPT(sys)   = Y_OPT(sys)   + del[MY]
	    MOS_ROT(sys) = MOS_ROT(sys) + del[MP]

# REVIEW: above sets del=0 for no solve; can remove this if()
	    if (solve_geom) {
		do m = 1, NCCD {
		    if (m == 3)				# Chip 3 is fiducial
			next
		    ndx = (m-1)*3 + NOPTVAR
# Must remove change to overall Mosaic; the sense should be correct
		    CN_XERR(sys,m) = CN_XERR(sys,m) + del[DX+ndx] #+ del[MX]
		    CN_YERR(sys,m) = CN_YERR(sys,m) + del[DY+ndx] #+ del[MY]
		    CN_RERR(sys,m) = CN_RERR(sys,m) + del[DZ+ndx] #+ del[MP]
		}
	    }

	    do m = 1, nset {
		ndx = (m-1) * 4 + NSYSVAR
		spar[CF,m] = spar[CF,m] + del[CF+ndx]
		spar[GX,m] = spar[GX,m] + del[GX+ndx]
		spar[GY,m] = spar[GY,m] + del[GY+ndx]
#		if (spar[GO,m] == 0.)
		if (spar[GO,m] == 0. || fix[GZ+ndx] == YES)
			del[GZ+ndx] = 0.
		spar[GZ,m] = spar[GZ,m] + del[GZ+ndx]
	    }
	}

# Free work vectors
call eprintf ("avg: %f %f\n")
call pargr (sumx / nslit)
call pargr (sumy / nslit)
	call gsfree (sfy)
	call gsfree (sfx)
	call asifree (asi)

### NB: because we have been working in ICS, must change predicted to mosaic
call eprintf ("Mosaic conversion at end of multi_solve...\n")
 	do k = 1, nslit {
 		n = ident_ccd (x[k], y[k])
 		call ics_to_ccd (x[k], y[k], ccd, n, x[k], y[k], stat)
		call mosim_coord(x[k], y[k], n)
# XXX REVIEW -- is this necessary??  Shouldn't these be flagged? Also, already+?
 		if (stat != ON_CHIP) {
 			x[k] = abs (x[k])
 			y[k] = abs (y[k])
#  			next
 		}
 	}
end




#
# G2_ELIM: procedure for gaussian elimination, n var's; double prec.
# (copy)

procedure g2_elim (a, n)

double	a[n+1,n]			# matrix to be solved
int	n				# number of variables

int	i, j, k
double	den, hold

begin
	do k = 1, n {
		den = a[k,k]
		if (den == 0.) {		# look for non-zero: switch
			do j = (k+1), n {
				if (a[k,k] != 0.) {
					do i = k, (n+1) {
						hold = a[i,j]
						a[i,j] = a[i,k]
						a[i,k] = hold
					}
				den = a[k,k]
				}
			}
			if (den == 0.)			# if still zero, skip
				next
		}
# call eprintf ("%d: %10.3g\n")
# call pargi(k)
# call pargd (den)
		do i = k, (n+1)
			a[i,k] = a[i,k] / den
		do j = 1, n
			if (j != k) {
				den = a[k,j]
				do i = k, (n+1)
					a[i,j] = a[i,j] - a[i,k] * den
			}
	}
end

#
# NUMDERIV: find numerical derivative by varying setup slightly. Works to ICS
#

procedure numderiv (varpar, delpar, xs, ys, w, asi, px, py, ccd, sys, e1, a2, as, a3, a4, fixed)

double	varpar				# system parameter to vary 
real	delpar				# amount to vary parameter
double	xs, ys				# x,y on slit (mm)
double	w				# ref wavelength
pointer	asi				# pointer to fit structure
real	px, py				# partials in x, y
real	ccd[NCCD,3]			# CCD geometry
double	sys[NPARAM]			# system parameters for pt_xfm
double	e1[3,3], a2[3,3], a3[3,3], a4[3,3]	# transforms
double	as[3,3]				# slider transform
int	fixed				# Is parameter fixed?

int	n			# CCD number

int	i, ndx
int	noff
double	par0				# original value (must be restored!)
double	vinc
real	x[NORD], y[NORD]		# work arrays

real	deriv[2]

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
		call setup (e1, a2, as, a3, a4, ccd, sys)
		call pt_xfm (xs, ys, w, e1, a2, as, a3, a4, ccd, sys, x[ndx], y[ndx], n)
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
# NUMDER_CCD: find numerical derivative for CCD params, since post-ICS
#

procedure numder_ccd (varpar, delpar, xs, ys, w, asi, px, py, ccd, sys, e1, a2, as, a3, a4, fixed)

double	varpar				# system parameter to vary 
real	delpar				# amount to vary parameter
double	xs, ys				# x,y on slit (mm)
double	w				# ref wavelength
pointer	asi				# pointer to fit structure
real	px, py				# partials in x, y
real	ccd[NCCD,3]			# CCD geometry
double	sys[NPARAM]			# system parameters for pt_xfm
double	e1[3,3], a2[3,3], a3[3,3], a4[3,3]	# transforms
double	as[3,3]				# slider transform
int	fixed				# Is parameter fixed?

int	n			# CCD number

int	i, ndx
int	stat
int	noff
double	par0				# original value (must be restored!)
double	vinc

real	xics, yics
real	x[NORD], y[NORD]		# work arrays

real	deriv[2]

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

# Solve for ICS location
	call setup (e1, a2, as, a3, a4, ccd, sys)
	call pt_xfm (xs, ys, w, e1, a2, as, a3, a4, ccd, sys, xics, yics, n)

# Now work with CCDs -- note that only ccd_geom must be called now ...
	ndx = 1
	do i = -noff, noff {
		varpar = par0 + vinc * i
		call ccd_geom (ccd, sys)
		call ics_to_ccd (xics, yics, ccd, n, x[ndx], y[ndx], stat)
		call mosim_coord(x[ndx], y[ndx], n)
		if (stat != ON_CHIP) {
			x[ndx] = abs (x[ndx])	## XXX ALREADY POSITIVE??
			y[ndx] = abs (y[ndx])
		}

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
# GET_CCDNUM: Hopefully tmp, designed to tell which chip & chip coords actual
# measured positions correspond to.  It is _loaded_ with HARDCODES!!
#

procedure	get_ccdnum (x, y, n, xccd, yccd)

real	x, y			# measured locations in image
int	n			# Chip number (returned)
real	xccd, yccd		# x,y on that CCD

begin
#	if (x < 49. || y < 1. || x > 8240. || y > 8192.) {
	if (x < 1. || y < 1. || x > 8192. || y > 8192.) {
		n = 0
call eprintf ("Whoops! `measurement' off image! (%f %f)\n")
			call  pargr (x)
			call  pargr (y)
		return
	}

	n = (x-49.)
	n = n / 2048 + 1

	xccd = x - 2048 * (n-1) - PRECOL

	if (y > 4096.) {
		n = n + 4
		yccd = y - 4096.
	} else {
		yccd = y
	}
end

#
# PA_ADJ: calculate PA-dependent adjustments for distortion 
#

procedure pa_adj (fd, pa, xobs, yobs, xadj, yadj, nslit, ccd)

pointer	fd				# file descriptor of mapping
real	pa				# PA in degrees
real	xobs[nslit], yobs[nslit]	# observed pixel locations
real	xadj[nslit], yadj[nslit]	# returned adjustment values
int	nslit				# number of slits
real	ccd[NCCD,3]			# CCD geometry

int	n
int	k
real	xdet, ydet			# correction for detector motion
real	xics, yics, xccd, yccd
real	xg, yg, xgdel, ygdel		# values in gravity system
real	theta, cost, sint
pointer	sfx, sfy			# pointers to x,y mappings

real	sum		# TMP

real	gseval()

begin
# Initialize the distortion mapping:
	call gs_ingest (fd, sfx, sfy)

# Set up the ICS-gravity xforms
# NB: check signs here, below
	theta = DEGTORAD(pa)
	cost = cos (theta)
	sint = sin (theta)

# OK, what about detector motion?  It works in opposite sense, has an ampl a,b
	xdet = -1.9 * cost
	ydet = -1.3 * sint
	
	do k = 1, nslit {
# convert observed into ICS:
		call get_ccdnum (xobs[k], yobs[k], n, xccd, yccd)
		call ccd_to_ics (xccd, yccd, ccd, n, xics, yics)
		xics = xics * 0.015	# convert to mm  TMP! HARDCODE
		yics = yics * 0.015	# convert to mm  TMP! HARDCODE

# ... convert into gravity system
		xg =  cost * xics + sint * yics
		yg = -sint * xics + cost * yics

# ... get the adjustments
		xgdel = gseval (sfx, xg, yg)
		ygdel = gseval (sfy, xg, yg)

# ... reconvert to ICS
		xadj[k] =  cost * xgdel - sint * ygdel + xdet
		yadj[k] =  sint * xgdel + cost * ygdel + ydet
	}

# TMP measure:  remove the average
	sum = 0.
	do k = 1, nslit
		sum = sum + xadj[k]
	sum = sum / nslit
	call aaddkr (xadj, -sum, xadj, nslit)
	sum = 0.
	do k = 1, nslit
		sum = sum + yadj[k]
	sum = sum / nslit
	call aaddkr (yadj, -sum, yadj, nslit)
call eprintf ("AVERAGES REMOVED!\n")

end


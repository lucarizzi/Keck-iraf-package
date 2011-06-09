# TDB immed: insert the NULL values appropriately
# OR, use the cfitsio routines
#
# TBD:  get various params into proper routines
# The malloc's in selector should be replaced by salloc's
# Note that slit-sep should be _added_ to all required lengths,
# The whole auto-selector issue needs to be revisited, with real lengths added.
# Resolution of slit lengths (extend, resolve conficts)
# Binding, slit editing
# Mapping slits to celestial coordinates  -- DONE!
# DEFINE the LEN issue

# Mapping slits to metal [mostly done] ****

# NB: it is probably best to get the central slit coords, etc, and then
# map onto the metal as we did before.  There are a few reasons for doing it
# this way, part of which is the curved-slit problem (assuming the object is
# more likely to be near the center than the edge of the slits); also, this
# is the direction things will likely go in the end; and 3rd it makes dealing
# with the PA=INDEF objects a bit easier, since they can be cut as vertical
# on the metal.

include	<math.h>
include	"deimos.h"
include	"dsimulator.h"

procedure	t_dsimulator()

char	objfile[SZ_FNAME]			# ID, RA, Dec, equ, prior,
char	output[SZ_FNAME]			# output file name
char	mdf[SZ_FNAME]				# mask design file name

pointer	fda, fdc				# in/out file descriptors

pointer	indat			# pointer to instr/tel data structure

int	ntarg
pointer	tdat			# pointer to data structure for targets

int	nslit
pointer	sdat			# pointer to data structure for slits

int	i	## TMP

bool	clgetb(), strne()
pointer	open()

begin
	call clgstr ("objfile", objfile, SZ_FNAME)
#	call clgstr ("output", output, SZ_FNAME)
	call clgstr ("mdf", mdf, SZ_FNAME)

# Read in telescope, default data:
	call data_init (indat)

# Open the input list of targets:
	fda = open (objfile, READ_ONLY, TEXT_FILE)

	fdc = open (mdf, NEW_FILE, TEXT_FILE)

# Read in the target data:
	call targ_init (fda, tdat, ntarg, indat)

# Calc the location of the telescope axis
	call fld2telax (indat)

# Calculate position wrt telescope axis:
	call tel_coords (tdat, ntarg, indat)


# Ready to display and be interactive
	call select_targ (tdat, ntarg, indat)

# Save output?

# Calculate the SLITS in celestial coordinates.
# It is here that slit lengths, edges must be resolved.  Does binding go here?
# Or is this a step that should actually be included above?
#

	call gen_slits (tdat, ntarg, sdat, nslit, indat)
#	call len_slits (tdat, ntarg, sdat, nslit, indat)
#####
call eprintf ("Slit lengthening LEFT OUT\007!")
#####


call eprintf (" calling sky_coord \n")
# Convert adopted slits into coordinated on the sky:
	call sky_coords (sdat, nslit, indat)



## For symmetry, we now recalculate the telescope coords and then mask coords.
## In the end, this operation may be performed at Keck immediately prior to
## mask fabrication.
#  Note that _refraction_ correction should probably go here, as a differential
#  correction, as it is tiny and specific to actual time and date of observation

call eprintf (" calling tel_coord \n")

# Calculate position wrt telescope axis:
	call tel_coords (sdat, nslit, indat)

 
 call eprintf (" calling mask_coords \n")
# Generate the mask coordinates for the slits.

 
 	call mask_coords (sdat, nslit)


 
#
call eprintf (" calling write_design ...\n")
#	call write_fitsio (mdf, indat, tdat, ntarg, sdat, nslit, objfile)
call eprintf ("DEBUG\n")
#	call write_design (indat, tdat, ntarg, sdat, nslit, fdc, objfile)
# call eprintf (" ASCII Surfcam Format ...\n")
 	call write_surfcam (sdat, nslit, fdc, "TEST_MASK1")

	call close (fdc)


# do i = 0, nslit-1 {
# call printf ("SLIT%d    %12.3h  %12.2h %7.2f  99 99  1 %6.2f \n")
# call pargi (INDEX(sdat,i))
# call pargd (RADTODEG(RA(sdat,i))/15.)
# call pargd (RADTODEG(DEC(sdat,i)))
# call pargr (2000.0)
# if (PA(sdat,i) != INDEF)
# 	call pargr (RADTODEG(PA(sdat,i)))
# call pargr (PA(sdat,i))
# }

#
call eprintf (" [recall: still must check proper x,y and tan mapping]\n")
#
end

	

#
# FLD2TELAX:  from field center and rotator PA, calc coords of telescope axis
#

procedure	fld2telax (indat)

pointer	indat

double	r, pa_fld

double	cosd, sind, cosa, sina, cost, sint, cosr, sinr

begin
# convert field center offset (arcsec) to radians
	r = DEGTORAD(sqrt (FLDCEN_X*FLDCEN_X + FLDCEN_Y*FLDCEN_Y) / 3600.)

# get PA of field center
	pa_fld = atan2 (FLDCEN_Y, FLDCEN_X)

	cosr = cos (r)
	sinr = sin (r)
	cosd = cos (DEC_FLD(indat))
	sind = sin (DEC_FLD(indat))

	cost = cos (PA_ROT(indat) + pa_fld)
	sint = sin (PA_ROT(indat) + pa_fld)

	sina = sinr * sint / cosd		# ASSUME not at dec=90
	cosa = sqrt (1. - sina*sina)

	RA_TEL(indat) = RA_FLD(indat) - asin (sina)
	DEC_TEL(indat) = asin ((sind*cosd*cosa - cosr*sinr*cost) /
				(cosr*cosd*cosa - sinr*sind*cost))

end

############################################################################
#

#
# SELECTOR: Does an auto selection of slits
# Should include an option for weighting to keep things toward the center.
# Note that y's sent to sel_rank are relative to starting y
# have it run from bottom to top
# -- need to work in segments to accomodate currently selected objects
# -- there was something else ...

procedure	selector (tdat, ntarg, nlist, minsep)

pointer	tdat
int	ntarg
int	nlist		# list to work on
real	minsep		# XXX min. separation -- probably should be in DEFDAT

int	nopt, npre		# number of options, prev. selected objects
int	i, ndx
int	ix			# starting index for search (saves time)
int	nselect				# Number of selected slits
real	xlow, xupp, xskip
pointer	bufx1, bufx2			# TMP? buffers for pre-sel. objs
pointer	bufn, bufx, bufp, bufsel	# TMP, for now

begin
	nopt = 0
	npre = 0
	do i = 0, ntarg-1 {
		if (SEL(tdat,i) == YES) {
			npre = npre + 1
		} else if (SAMPL(tdat,i) == nlist && STAT(tdat,i) == YES) {
			nopt = nopt + 1
		}
	}
	call malloc (bufx1, npre, TY_INT)
	call malloc (bufx2, npre, TY_REAL)
	call malloc (bufn, nopt, TY_INT)
	call malloc (bufx, nopt, TY_REAL)
	call malloc (bufp, nopt, TY_INT)
	call malloc (bufsel, nopt, TY_INT)

# Grep on previously selected objects and suitable options; fill vectors
	nopt = 0
	npre = 0
	do i = 0, ntarg-1 {
		if (SEL(tdat,i) == YES) {
			Memr[bufx1+npre] = XARCS(tdat,i) - LEN1(tdat,i)
			Memr[bufx2+npre] = XARCS(tdat,i) + LEN2(tdat,i)
			npre = npre + 1
		} else if (SAMPL(tdat,i) == nlist && STAT(tdat,i) == YES) {
			Memi[bufn+nopt] = INDEX(tdat,i)
			Memr[bufx+nopt] = XARCS(tdat,i)
			Memi[bufp+nopt] = PCODE(tdat,i)
			nopt = nopt + 1
		}
	}

# Sort the two lists
	call sel_sort (Memr[bufx1], Memr[bufx2], npre,
				Memi[bufn], Memr[bufx], Memi[bufp], nopt)

# The number of "gaps" to search is npre+1
	ndx = 0
	xlow = XLOW_LIM
	xskip = 0.
	nselect = 0
	if (nopt > 0) {
	    do i = 0, npre {
		if (i < npre) {
			xupp = Memr[bufx1+i]
			xskip = Memr[bufx2+i] - Memr[bufx1+i]
		} else {
			xupp = XUPP_LIM
		}
		if (xupp <= xlow)

			next

		call sel_rank (Memr[bufx], Memi[bufp], Memi[bufn],
		Memi[bufsel], nopt, ix, xlow, xupp, minsep, nselect)

		xlow = xupp + xskip
	    }
	}


#...select the mask slits
	if (nselect > 0) {
		do i = 0, nselect-1 {
			SEL(tdat,Memi[bufsel+i]) = YES
		}
	}


	call mfree (bufsel, TY_INT)
	call mfree (bufp, TY_INT)
	call mfree (bufx, TY_REAL)
	call mfree (bufn, TY_INT)
	call mfree (bufx2, TY_REAL)
	call mfree (bufx1, TY_INT)
end


# SEL_RANK: Select slits with priority ranking.
# The scheme is to find the next possible slit, and then to look up to one
# min-slit width away for higher-priority objects. The higher priorities are
# down-weighted depending on their distance.
#

procedure	sel_rank (x, pri, index, sel, npt, isel, x1, x2, minsep, nsel)

real	x[npt]				# Y (rel. to search start)
int	pri[npt]			# Priority
int	index[npt]			# Index of objects
int	sel[npt]			# selected objects
int	npt				# Number of objects
real	x1, x2				# xrange to fill
real	minsep				# minimum separation (arcsec)
int	nsel				# Number of selected objects
int	isel 				# starting index

int	i, j

real	xnext, xlook, xstop, xlast
real	prisel, prinorm

begin
# make sure we initialize to the right index ...
	if (nsel <= 0) {
		isel = 0
		nsel = 0
	}

	if (x2 - x1 < minsep)		# probably too restrictive
		return
					# There is no slit len in the routine!


# Start at half a slit length; stop inside half slit length
	xstop = min (x[npt], x2-0.5*minsep)
	xnext = x1 + 0.5 * minsep
	xlast = x1

# Loop through to end
	for (i = isel + 1; i <= npt; i = i + 1) {
		if (x[i] < xnext)
			next

		if (x[i] > xstop) {
			isel = i - 1
			break
		}

		isel = i
		prisel = pri[isel] / (x[isel] - xlast)
# Now look for higher priority to win out, over range (x_curr,xlook)
		xlook = min (x[isel]+minsep, xstop)
		if (isel < npt) {
			do j = isel+1, npt {
				if (x[j] >= xlook)
					break
				prinorm = pri[j] / (x[j] - xlast)
				if (prinorm > prisel) {
					isel = j
					prisel = prinorm
				}
			}
		}

		nsel = nsel + 1
		sel[nsel] = index[isel]
		xlast = x[isel]
		xnext = xlast + minsep
		i = isel			# Reset search start point
	}
end

#
# SEL_SORT: sort the selected/selection lists
#

procedure	sel_sort (px1, px2, npre, index, x, pri, nopt)

real	px1[npre], px2[npre]		# x-pre limits to sort
int	npre				# Number of prev. selections

int	index[nopt]			# Index of optional objects
real	x[nopt]				# x-opt list to sort
int	pri[nopt]			# Priority
int	nopt				# Number of optional objects

int	i, j
int	ihold, phold
real	xhold

begin

# Sort the preselected list in x (low-to-high)

	do i = 1, npre-1 {
	    do j = 1, npre-i {
		if (px1[j] > px1[j+1]) {
			xhold = px1[j+1]
			px1[j+1] = px1[j]
			px1[j] = xhold
			xhold = px2[j+1]
			px2[j+1] = px2[j]
			px2[j] = xhold
		}
	    }
	}

# Sort the list of optional objects in x (low-to-high)
	do i = 1, nopt-1 {
	    do j = 1, nopt-i {
		if (x[j] > x[j+1]) {
			xhold = x[j+1]
			x[j+1] = x[j]
			x[j] = xhold
			ihold = index[j+1]
			index[j+1] = index[j]
			index[j] = ihold
			phold = pri[j+1]
			pri[j+1] = pri[j]
			pri[j] = phold
		}
	    }
	}
end


#
# MASK_COORDS:  Convert (x,y) on sky to xmm,ymm on slitmask
# NB: Requires instrument params, INCLUDING THE APPROP. TEL FOC LENGTH
# Assumes that the XARCS,YARCS of the slit ends are the tan projections
#

procedure	mask_coords (sdat, nslit)

pointer	sdat
int	nslit

double	xfp, yfp			# x,y points in FP (tan projection)
double	xsm, ysm			# x,y points on the mask
double	pa

int	i
double	xoff, yoff			# offset, telaxis to origin of slitmask
real	sina, cosa

begin

# offset from telscope axis to slitmask origin, IN SLITMASK COORDS
	yoff = ZPT_YM * (1. - cos (DEGTORAD(M_ANGLE)))
	yoff = 0.	# XXX check!  Am not sure where the above comes from
	xoff = 0.

	do i = 0, nslit-1 {
# XXX For now, carry through the RELPA thing; in end, must be specified!
			if (RELPA(sdat,i) != INDEF) {
				cosa = cos (RELPA(sdat,i))
				sina = sin (RELPA(sdat,i))
			} else {
				cosa = 1.
				sina = 0.
			}
#			cosa = cos (RELPA(sdat,i))	# XXX 
#			sina = sin (RELPA(sdat,i))	# XXX 

			X1(sdat,i) = XARCS(sdat,i) - LEN1(sdat,i) * cosa
			Y1(sdat,i) = YARCS(sdat,i) - LEN1(sdat,i) * sina
			X2(sdat,i) = XARCS(sdat,i) + LEN2(sdat,i) * cosa
			Y2(sdat,i) = YARCS(sdat,i) + LEN2(sdat,i) * sina

# XXX cuidado!  I am not sure that the tan-projection of the rel PA is the
# same as the rel PA -- MUST CHECK! (This code comes from gen_slits)


# The focal plane coordinates are now simply a tan projection of (x,y) arcsec
# Need to verify that these are truly symmetric:
#		xfp = FL_TEL * tan (DEGTORAD(X1(sdat,i)/3600.))
#		yfp = FL_TEL * tan (DEGTORAD(Y1(sdat,i)/3600.)) / cos (DEGTORAD(X1(sdat,i)/3600.))


# X1,Y1 are now tan projections already!

		xfp = FL_TEL *  X1(sdat,i) / 206264.8D0
		yfp = FL_TEL * (Y1(sdat,i) - 0.5*SLWID(sdat,i)) / 206264.8D0
		pa = 0.
		call gnom_to_dproj (xfp, yfp, xfp, yfp)		# (allowed)
		call proj_to_mask (xfp, yfp, pa, xsm, ysm, pa)

		XMM1(sdat,i) = xsm + xoff
		YMM1(sdat,i) = ysm + yoff

		xfp = FL_TEL *  X2(sdat,i) / 206264.8D0
		yfp = FL_TEL * (Y2(sdat,i) - 0.5*SLWID(sdat,i)) / 206264.8D0
		pa = 0.
		call gnom_to_dproj (xfp, yfp, xfp, yfp)		# (allowed)
		call proj_to_mask (xfp, yfp, pa, xsm, ysm, pa)

		XMM2(sdat,i) = xsm + xoff
		YMM2(sdat,i) = ysm + yoff

		xfp = FL_TEL *  X2(sdat,i) / 206264.8D0
		yfp = FL_TEL * (Y2(sdat,i) + 0.5*SLWID(sdat,i)) / 206264.8D0
		pa = 0.
		call gnom_to_dproj (xfp, yfp, xfp, yfp)		# (allowed)
		call proj_to_mask (xfp, yfp, pa, xsm, ysm, pa)

		XMM3(sdat,i) = xsm + xoff
		YMM3(sdat,i) = ysm + yoff

		xfp = FL_TEL *  X1(sdat,i) / 206264.8D0
		yfp = FL_TEL * (Y1(sdat,i) + 0.5*SLWID(sdat,i)) / 206264.8D0
		pa = 0.
		call gnom_to_dproj (xfp, yfp, xfp, yfp)		# (allowed)
		call proj_to_mask (xfp, yfp, pa, xsm, ysm, pa)

		XMM4(sdat,i) = xsm + xoff
		YMM4(sdat,i) = ysm + yoff
	}

## Perhaps we want to force YMM4-YMM1 == YMM3-YMM2; for non-tilted slits,
## this should produce a cleaner edge; otherwise, jumps can occur.

end



#
# SKY_COORDS: Convert xarcs,yarcs in tel coords onto sky
#   Note that this routine, called infrequently, does not need to be efficient.
#

procedure	sky_coords (sdat, nslit, indat)

pointer	sdat
int	nslit
pointer	indat

int	i
double	r			# radius of object from tel-axis
double	phi			# PA on sky from tel-axis
double	sind, sina		# sin of dec, delta-RA
double	ra0, dec0, pa0		# RA, Dec and PA on axis

double	x, y

begin
	ra0  = RA_TEL(indat)
	dec0 = DEC_TEL(indat)
	pa0  = PA_ROT(indat)

	do i = 0, nslit-1 {
		x = 0.5 * (X1(sdat,i) + X2(sdat,i))
		y = 0.5 * (Y1(sdat,i) + Y2(sdat,i))

		r = sqrt (x*x + y*y)
		r = atan (r/206264.8D0)

		phi = atan2 (y, x) + pa0

		sind = sin (dec0) * cos (r) + cos (dec0) * sin (r) * cos (phi)

		sina = sin (r) * sin (phi) / sqrt (1. - sind*sind)

		DEC(sdat,i) = asin (sind)
		RA(sdat,i) = ra0 + asin (sina)
# PA(sdat,i) = already assigned 

# calc the centers and lengths of the slits

		XARCS(sdat,i) = 0.5 * (X1(sdat,i) + X2(sdat,i))
		YARCS(sdat,i) = 0.5 * (Y1(sdat,i) + Y2(sdat,i))

# XXX NB: by convention, slit length will be defined as TOTAL length
		x = X2(sdat,i) - X1(sdat,i)
		y = Y2(sdat,i) - Y1(sdat,i)
		LEN1(sdat,i) = 0.5 * sqrt (x*x + y*y)
		LEN2(sdat,i) = LEN1(sdat,i)

	}
end

#
# WRITE_DESIGN: generate the FITS table design
#

include	<time.h>
include	"ftab.h"

procedure	write_design (indat, tdat, ntarg, sdat, nslit, fd, objfile)

pointer	indat			# pointer to instr/telesc data struct

pointer	tdat			# pointer to targets data struct
int	ntarg

pointer	sdat			# pointer to slits data struct
int	nslit

pointer	fd			# output file descriptor
char	objfile[ARB]		# name of input file

int	i

char	ident[SZ_ID]
int	stat
int	nselect


char	fmtstr[SZ_LINE]				# format for writing table row
int	fmlen
char	date[SZ_TIME]					# time/date of file
pointer	kdat				# FITS table data structure
# int	inul
# real	rnul
# double	dnul
char	snul[1]
#####
char	creatask[60]
char	person[30]
char	unknown[3]

char	acode

char	fmt_mm[4], mm[2], deg[3]	# define to avoid filling declaration space

char	cnul[6]
int	inul
real	rnul
double	dnul

int	sscan()
begin
	call strcpy ("F9.3", fmt_mm, 4)
	call strcpy ("mm", mm, 2)
	call strcpy ("deg", deg, 3)

	fmlen = SZ_LINE

	call strcpy ("INDEF", cnul, 6)
	inul = -9999
	rnul = -9999.0		# Careful! will write it out longer if it
	dnul = -9999.0		# doesn't fit into the field

# SET UP FITS FILE
	call strcpy ("", snul, 1)
	call strcpy ("Phillips <phillips@ucolick.org>", person, 35)
	call strcpy ("DSIMULATOR -- June2000", creatask, 35)
	call strcpy ("???", unknown, 3)

	call ft_date (date, SZ_TIME)

# Write the primary HDU
	call ftab_whdu0 (fd)

call eprintf ("    table 1 ")
# Only include selected objects
	nselect = 0
	do i = 0, ntarg-1 {
		if (SEL(tdat,i) == YES)
			nselect = nselect + 1
	}

########  WRITE ObjectCat TABLE
	call ftab_init (kdat, 18, fmtstr)

	call ftcol_defi (kdat, "ObjectId",  "I6", snul, inul,  6, fmtstr, fmlen)
	call ftcol_defc (kdat, "OBJECT",   "A68", snul, cnul, 68, fmtstr, fmlen)
	call ftcol_defd (kdat, "RA_OBJ", "F12.8", deg, dnul, 12, fmtstr, fmlen)
	call ftcol_defd (kdat, "DEC_OBJ","F12.8", deg, dnul, 12, fmtstr, fmlen)
	call ftcol_defc (kdat, "RADESYS",   "A8", snul, cnul, 8, fmtstr, fmlen)
	call ftcol_defd (kdat, "EQUINOX", "F8.3",   "a", dnul,  8, fmtstr, fmlen)
	call ftcol_defd (kdat, "MJD-OBS", "F11.3",  "d", dnul, 11, fmtstr, fmlen)
	call ftcol_defr (kdat, "mag",      "F7.3", snul, rnul,  7, fmtstr, fmlen)
	call ftcol_defc (kdat, "pBand",      "A6", snul, cnul,  6, fmtstr, fmlen)
	call ftcol_defr (kdat, "RadVel", "F10.3",snul, rnul, 10, fmtstr, fmlen)
	call ftcol_defr (kdat, "MajAxis", "F9.2", "arcsec", rnul, 9, fmtstr, fmlen)
	call ftcol_defr (kdat, "MajAxPA", "F8.2",   deg, rnul, 8, fmtstr, fmlen)
	call ftcol_defr (kdat, "MinAxis", "F9.2", "arcsec", rnul, 9, fmtstr, fmlen)
	call ftcol_defr (kdat, "PM_RA",  "F9.4", "arcsec/a", rnul, 9, fmtstr, fmlen)
	call ftcol_defr (kdat, "PM_Dec", "F9.4", "arcsec/a", rnul, 9, fmtstr, fmlen)
	call ftcol_defr (kdat, "Parallax", "F7.4", "arcsec", rnul, 7, fmtstr, fmlen)
	call ftcol_defc (kdat, "ObjClass",  "A20", snul, cnul, 20, fmtstr, fmlen)
	call ftcol_defi (kdat, "CatFilePK",  "I6", snul, inul,  6, fmtstr, fmlen)

	call ftab_whead (fd, kdat, nselect)

# Boilerplate keywords 1
	call bpstd (fd, kdat, "ObjectCat", 1, 1, "ObjectId", snul)
	call bpadd (fd, kdat, 1, "CatFilePK", "CatFilePK", "CatFiles")

	call pkwstr (fd, kdat, "CATNAME", objfile, snul)

# Write out the OBJECT CAT TABLE
		
	do i = 0, ntarg-1 {
		if (SEL(tdat,i) == NO)
			next
		stat = sscan (DATLINE(tdat,i))
		call gargwrd (ident, SZ_ID)

		call sprintf (BUFROW(kdat), TWUSED(kdat), fmtstr)
			call pargi (INDEX(tdat,i))
			call pargstr (ident)
			call pargd   (RADTODEG(RA(tdat,i)))
			call pargd   (RADTODEG(DEC(tdat,i)))
			call pargstr ("")	# RADESYS
			call pargd   (STD_EQX(indat))	# XXX
			call pargr   (rnul)	#
			call pargr   (rnul)	# XXX mag
			call pargstr (cnul)	#pband
			call pargr   (rnul)	#
			call pargr   (rnul)	# maj ax
			if (PA(tdat,i) == INDEF)	# XXX should be specific
				call pargr (rnul)
			else
				call pargr (RADTODEG(PA(tdat,i)))
			call pargr   (rnul)	# min_ax
			call pargr   (rnul)	# PM_ra
			call pargr   (rnul)	# PM_dec
			call pargr   (rnul)	# Parallax
			if (PCODE(tdat,i) == CODE_GS)
				call pargstr ("Guide_Star")
			else if (PCODE(tdat,i) == CODE_AS)
				call pargstr ("Alignment_Star")
			else
				call pargstr ("Program_Target")
			call pargi   (1)
			
		call ftab_wrow (fd, kdat)
	}
	call ftab_free (fd, kdat)

call eprintf ("OK\n    table 2 ")
##### WRITE CatFiles TABLE
	call ftab_init (kdat, 2, fmtstr)

	call ftcol_defi (kdat, "CatFilePK",     "I6", snul, inul,   6, fmtstr, fmlen)
	call ftcol_defc (kdat, "CatFileName", "A255", snul, cnul, 255, fmtstr, fmlen)
#	call ftcol_defc (kdat, "CatFileName", "A68", snul, cnul, 68, fmtstr, fmlen)

	call ftab_whead (fd, kdat, 1)

# Boilerplate 2
	call bpstd (fd, kdat, "CatFiles", 1, 0, "CatFilePK", snul)

# Write out the OBJECT CAT TABLE
	call sprintf (BUFROW(kdat), TWUSED(kdat), fmtstr)
		call pargi (1)
		call pargstr (cnul)
			
	call ftab_wrow (fd, kdat)

	call ftab_free (fd, kdat)


call eprintf ("OK\n    table 3 ")
########  WRITE MaskDesign TABLE
	call ftab_init (kdat, 12, fmtstr)

	call ftcol_defi (kdat, "DesId",     "I11", snul, inul, 11, fmtstr, fmlen)
	call ftcol_defc (kdat, "DesName",   "A68", snul, cnul, 68, fmtstr, fmlen)
	call ftcol_defc (kdat, "DesAuth",   "A68", snul, cnul, 68, fmtstr, fmlen)
	call ftcol_defc (kdat, "DesCreat",  "A68", snul, cnul, 68, fmtstr, fmlen)
	call ftcol_defc (kdat, "DesDate",   "A19", snul, cnul, 19, fmtstr, fmlen)
	call ftcol_defd (kdat, "RA_PNT",  "F12.8", deg, dnul, 12, fmtstr, fmlen)
	call ftcol_defd (kdat, "DEC_PNT", "F12.8", deg, dnul, 12, fmtstr, fmlen)
	call ftcol_defc (kdat, "RADEPNT",    "A8", snul, cnul,  8, fmtstr, fmlen)
	call ftcol_defd (kdat, "EQUINPNT","F13.6",   "a", dnul, 13, fmtstr, fmlen)
	call ftcol_defd (kdat, "PA_PNT",  "F12.7", deg, dnul, 12, fmtstr, fmlen)
	call ftcol_defc (kdat, "DATE_PNT",  "A19", snul, cnul, 19, fmtstr, fmlen)
	call ftcol_defd (kdat, "LST_PNT",   "F8.3", deg, dnul,  8, fmtstr, fmlen)

	call ftab_whead (fd, kdat, 1)

# Boilerplate 3
	call bpstd (fd, kdat, "MaskDesign", 1, 0, "DesId", snul)

	call sprintf (BUFROW(kdat), TWUSED(kdat), fmtstr)
		call pargi (1)
		call pargstr ("output")	# XXX
		call pargstr (person)
		call pargstr (creatask)
		call pargstr (date)
		call pargd (RADTODEG(RA_TEL(indat)))
		call pargd (RADTODEG(DEC_TEL(indat)))
		call pargstr ("")
		call pargd (STD_EQX(indat))
		call pargd (RADTODEG(PA_ROT(indat)))
		call pargstr (date)		# XXX
		call pargd (RADTODEG(HA_TEL(indat)))
			
	call ftab_wrow (fd, kdat)

	call ftab_free (fd, kdat)


call eprintf ("OK\n    table 4 ")
########  WRITE DesiSlits TABLE

	call ftab_init (kdat, 9, fmtstr)

	call ftcol_defi (kdat, "dSlitId", "I11",  snul, inul, 11, fmtstr, fmlen)
	call ftcol_defi (kdat, "DesId",   "I11",  snul, inul, 11, fmtstr, fmlen)
	call ftcol_defr (kdat, "slitRA", "F12.8", deg, rnul, 12, fmtstr, fmlen)
	call ftcol_defd (kdat, "slitDec", "F12.8", deg, dnul, 12, fmtstr, fmlen)
	call ftcol_defc (kdat, "slitTyp",  "A1",  snul, snul,  1, fmtstr, fmlen)
	call ftcol_defr (kdat, "slitLen", "F11.3", "arcsec", rnul, 11, fmtstr, fmlen)
	call ftcol_defr (kdat, "slitLPA", "F8.3", deg, rnul,  8, fmtstr, fmlen)
	call ftcol_defr (kdat, "slitWid", "F11.3", "arcsec", rnul, 11, fmtstr, fmlen)
	call ftcol_defr (kdat, "slitWPA", "F8.3", deg, rnul,  8, fmtstr, fmlen)

	call ftab_whead (fd, kdat, nslit)

# Boilerplate 4
	call bpstd (fd, kdat, "DesiSlits", 1, 1, "dSlitID", "")
	call bpadd (fd, kdat, 1, "DesId", "DesId", "MaskDesign")

	do i = 0, nslit-1 {

		if (PCODE(sdat,i) == CODE_GS || PCODE(sdat,i) == CODE_RF)
			acode = 'G'
		else if (PCODE(sdat,i) == CODE_AS)
			acode = 'A'
		else
			acode = 'P'
		
		call sprintf (BUFROW(kdat), TWUSED(kdat), fmtstr)
			call pargi (INDEX(sdat,i))
			call pargi (1)
			call pargd (RADTODEG(RA(sdat,i)))
			call pargd (RADTODEG(DEC(sdat,i)))
			call pargc (acode)
			call pargr (LEN1(sdat,i)+LEN2(sdat,i))
			if (PA(sdat,i) == INDEF)	# XXX should be specific
				call pargr (0.)
			else
				call pargr (RADTODEG(PA(sdat,i)))
			call pargr (SLWID(sdat,i))
			call pargd (RADTODEG(PA_ROT(indat))+90.)	# XXX
			
		call ftab_wrow (fd, kdat)
	}

	call ftab_free (fd, kdat)


call eprintf ("OK\n    table 5 ")
########  WRITE SlitObjMap TABLE

	call ftab_init (kdat, 3, fmtstr)

	call ftcol_defi (kdat, "DesId",    "I11", snul, inul, 11, fmtstr, fmlen)
	call ftcol_defi (kdat, "ObjectId", "I11", snul, inul, 11, fmtstr, fmlen)
	call ftcol_defi (kdat, "dSlitId",  "I11", snul, inul, 11, fmtstr, fmlen)

	call ftab_whead (fd, kdat, nselect)

# Boilerplate 5
	call bpstd (fd, kdat, "SlitObjMap", 2, 3, "DesId", "ObjectId")
	call bpadd (fd, kdat, 1, "DesId", "DesId", "MaskDesign")
	call bpadd (fd, kdat, 2, "dSlitId", "dSlitId", "DesiSlits")
	call bpadd (fd, kdat, 3, "ObjectId", "ObjectId", "ObjectCat")


	do i = 0, ntarg-1 {
		if (SEL(tdat,i) == NO)
			next
	    call sprintf (BUFROW(kdat), TWUSED(kdat), fmtstr)
		call pargi (1)
		call pargi (i)
		call pargi (SLNDX(tdat,i))
			
	    call ftab_wrow (fd, kdat)
	}

	call ftab_free (fd, kdat)


call eprintf ("OK\n    table 6 ")
########  WRITE MaskBlu TABLE
	call ftab_init (kdat, 16, fmtstr)

	call ftcol_defi (kdat, "BluId",     "I11", snul, inul, 11, fmtstr, fmlen)
	call ftcol_defi (kdat, "DesId",     "I11", snul, inul, 11, fmtstr, fmlen)
	call ftcol_defc (kdat, "BluName",   "A68", snul, cnul, 68, fmtstr, fmlen)
	call ftcol_defc (kdat, "BluObsvr",  "A68", snul, cnul, 68, fmtstr, fmlen)
	call ftcol_defc (kdat, "BluCreat",  "A68", snul, cnul, 68, fmtstr, fmlen)
	call ftcol_defc (kdat, "BluDate",   "A19", snul, cnul, 19, fmtstr, fmlen)
	call ftcol_defd (kdat, "LST_Use",  "F8.3",deg, dnul,  8, fmtstr, fmlen)
	call ftcol_defc (kdat, "Date_Use",  "A19", snul, cnul, 19, fmtstr, fmlen)
	call ftcol_defc (kdat, "TELESCOP",  "A68", snul, cnul, 68, fmtstr, fmlen)
	call ftcol_defc (kdat, "RefrAlg",   "A68", snul, cnul, 68, fmtstr, fmlen)
	call ftcol_defd (kdat, "AtmTempC", "F5.1","degC", dnul, 5, fmtstr, fmlen)
	call ftcol_defd (kdat, "AtmPres",  "F6.1","mmHg", dnul, 6, fmtstr, fmlen)
	call ftcol_defd (kdat, "AtmHumid", "F5.3", "K/m", dnul, 5, fmtstr, fmlen)
	call ftcol_defd (kdat, "AtmTTLap", "F7.5", snul, dnul,  7, fmtstr, fmlen)
	call ftcol_defd (kdat, "RefWave",  "F7.2", "nm", dnul,  7, fmtstr, fmlen)
	call ftcol_defc (kdat, "DistMeth",  "A68", snul, cnul, 68, fmtstr, fmlen)

	call ftab_whead (fd, kdat, 1)

# Boilerplate 6
	call bpstd (fd, kdat, "MaskBlu", 1, 1, "BluId", snul)
	call bpadd (fd, kdat, 1, "DesId", "DesId", "MaskDesign")

	call sprintf (BUFROW(kdat), TWUSED(kdat), fmtstr)
		call pargi (1)
		call pargi (1)
		call pargstr ("output")		# XXX
		call pargstr (person)
		call pargstr (creatask)
		call pargstr (date)
		call pargd   (RADTODEG(HA_TEL(indat)+RA_TEL(indat)))
		call pargstr (date)		# Date_use XXX
		call pargstr ("Keck II")
		call pargstr (cnul)		# RefrAlg
		call pargd   (TEMP(indat))
		call pargd   (PRES(indat))
		call pargd   (0.)		# TMP: no clue  AtmHumid XXX
		call pargd   (0.)		# TMP: no clue  AtmTTLap XXX
		call pargd   (WAVER(indat))
		call pargstr (cnul)		# DistAlg
			
	call ftab_wrow (fd, kdat)

	call ftab_free (fd, kdat)


call eprintf ("OK\n    table 7 ")
##### write the BluSlits table
	call ftab_init (kdat, 11, fmtstr)

	call ftcol_def (kdat, "bSlitId", "I11", snul, inul, 11, fmtstr, fmlen)
	call ftcol_def (kdat, "BluId",   "I11", snul, inul, 11, fmtstr, fmlen)
	call ftcol_def (kdat, "dSlitId", "I11", snul, inul, 11, fmtstr, fmlen)
	call ftcol_def (kdat, "slitX1", fmt_mm, mm, dnul, 9, fmtstr, fmlen)
	call ftcol_def (kdat, "slitY1", fmt_mm, mm, dnul, 9, fmtstr, fmlen)
	call ftcol_def (kdat, "slitX2", fmt_mm, mm, dnul, 9, fmtstr, fmlen)
	call ftcol_def (kdat, "slitY2", fmt_mm, mm, dnul, 9, fmtstr, fmlen)
	call ftcol_def (kdat, "slitX3", fmt_mm, mm, dnul, 9, fmtstr, fmlen)
	call ftcol_def (kdat, "slitY3", fmt_mm, mm, dnul, 9, fmtstr, fmlen)
	call ftcol_def (kdat, "slitX4", fmt_mm, mm, dnul, 9, fmtstr, fmlen)
	call ftcol_def (kdat, "slitY4", fmt_mm, mm, dnul, 9, fmtstr, fmlen)

	call ftab_whead (fd, kdat, nslit)

# Boilerplate 7
	call bpstd (fd, kdat, "BluSlits", 1, 2, "bSlitId", "")
	call bpadd (fd, kdat, 1, "BluId", "BluId", "MaskBlu")
	call bpadd (fd, kdat, 2, "dSlitId", "dSlitId", "DesiSlits")

	do i = 0, nslit-1 {
		call sprintf (BUFROW(kdat), TWUSED(kdat), fmtstr)
			call pargi (i)			# XXX
			call pargi (1)			# XXX
			call pargi (i)			# XXX
			call pargd (XMM1(sdat,i))
			call pargd (YMM1(sdat,i))
			call pargd (XMM2(sdat,i))
			call pargd (YMM2(sdat,i))
			call pargd (XMM3(sdat,i))
			call pargd (YMM3(sdat,i))
			call pargd (XMM4(sdat,i))
			call pargd (YMM4(sdat,i))
		call ftab_wrow (fd, kdat)
	}

	call ftab_free (fd, kdat)


call eprintf ("OK\n    table 8 ")
### append the RDBmap:
	call rdb_map (fd)


call eprintf ("OK\n")

end

#
# RDB_MAP: Append pre-existing RDBmap
#

procedure	rdb_map (fd)

pointer	fd			# output file descriptor

pointer	fdm			# file descriptor for RDBmap

pointer	open()
begin
	fdm = open (RDB_MAP_FILE, READ_ONLY, TEXT_FILE)

	call fcopyo (fdm, fd)

	call close (fdm)
end

#
# SLIT_SORT: sort the slits in x; keep track of the indices only
#

procedure	slit_sort (x, index, n)

real	x[n]				# x-pre limits to sort
int	index[n]			# Index of slit
int	n				# Number of prev. selections


int	i, j
int	ihold
real	xhold

begin

# Sort the slit list in x (low-to-high)

	do i = 1, n-1 {
	    do j = 1, n-i {
		if (x[j] > x[j+1]) {
			xhold = x[j+1]
			x[j+1] = x[j]
			x[j] = xhold
			ihold = index[j+1]
			index[j+1] = index[j]
			index[j] = ihold
		}
	    }
	}
end

#
# LEN_SLITS: adjust slit lengths to fit -- perhaps should be integral part of
# gen_slits
#

procedure	len_slits (tdat, ntarg, sdat, nslit, indat)

pointer	tdat			# pointer to target struct
int	ntarg

pointer	sdat			# pointer to slits struct
int	nslit

pointer	indat			# pointer to instrument struct

int	i
pointer	sp
pointer	bufx			# pointer to x-pos
pointer	bufi			# pointer to index vector

int	pc1, pc2
int	ndx1, ndx2
real	xlow, xupp, xcen
real	del1, del2, tana

begin
	call smark (sp)
	call salloc (bufx, nslit, TY_REAL)
	call salloc (bufi, nslit, TY_INT)

	call amovr (XARCS(sdat,0), Memr[bufx], nslit)
	call amovi (INDEX(sdat,0), Memi[bufi], nslit)

	call slit_sort (Memr[bufx], Memi[bufi], nslit)

## XXX must add extension to mask edge in here.

	do i = 0, nslit-2 {
		ndx1 = Memi[bufi+i]
		ndx2 = Memi[bufi+i+1]
		pc1 = PCODE(sdat,ndx1)
		pc2 = PCODE(sdat,ndx2)

# If both are alignment boxes, just go on ...
		if (pc1 == CODE_AS && pc2 == CODE_AS)		# no problem
			next

# We will need to recalculate something ...

		xlow = X2(sdat,ndx1) + SLIT_GAP(indat)
		xupp = X1(sdat,ndx2) - SLIT_GAP(indat)
		xcen = 0.5 * (xlow + xupp)

		if (pc1 == CODE_AS) {
			del1 = 0.
			del2 = X1(sdat,ndx2) - xlow
		} else if (pc2 == CODE_AS) {
			del1 =  xupp - X2(sdat,ndx1)
			del2 = 0.
		} else {
			del1 = xcen - 0.5*SLIT_GAP(indat) - X2(sdat,ndx1)
			del2 = X1(sdat,ndx2) - (xcen + 0.5*SLIT_GAP(indat))
		}

		X2(sdat,ndx1) = X2(sdat,ndx1) + del1
		if (del1 != 0. && RELPA(sdat,ndx1) != INDEF) {
			tana = tan (RELPA(sdat,ndx1))
			Y2(sdat,ndx1) = Y2(sdat,ndx1) + del1 * tana
		}

		X1(sdat,ndx2) = X1(sdat,ndx2) - del2
		if (del2 != 0. && RELPA(sdat,ndx2) != INDEF) {
			tana = tan (RELPA(sdat,ndx2))
			Y1(sdat,ndx2) = Y1(sdat,ndx2) - del2 * tana
		}

# The centers and lengths have now changed -- but defer to sky_coords()

# XXX report status here
#		if (del1 < 0) {			# conflict
#		} else if (xlow < xupp) {		# lengthen
#		}
		
	}

# XXX Check against objects, etc, and reflect status
		
	call sfree (sp)
end

# NB: problem with CVACCUM -- weights MUST be a variable as it is CHANGED in
# place if weighting != WTS_USER

# MEMORY PROBLEM -- probably cvfit values out-of-bounds on everlapping slits
# NB: PROBLEM WITH CVACCUM -- writes wherever in memory if x<xmin or x>xmax

# We want to change the sky-spectrum build-up onto a _slit-by-slit_ basis.
# This should be (a) faster; (b) much better at handling overlapping slitlets
# and (c) should handle overlaps _properly_.

# QLook tests

include	<imhdr.h>
include	<imio.h>
include	<math.h>
include	<math/curfit.h>
include "instrument.h"
include "redux.h"
include "fitsio.h"

define	NORDSKY	5000
define	NORDSLT	4
define	NSLTPTS	13		# number of points to fit in slit-mm mapping

define	NWAV	6	# XXX TMP! number of wavelengths for plotting
define	NWSTEP	40	# Number of steps to fit slit edge curves

procedure	t_rtest ()

char	fname[SZ_FNAME]
char	ximage[SZ_FNAME]
char	oimage[SZ_FNAME]
int	ichip1, ichip2			# First, last chip paris for processing
bool	regions_only			# stop at generating regions file?

char	guiname[SZ_FNAME]
pointer	im, imo, imx
int	i
int	nslits
pointer	sdat				# pointer to slit info from header
pointer	cdat				# pointer to slit info per chip
# pointer	fd			# output file descriptor



pointer	fda, fdb
pointer	fdc
pointer	map[8]				# pointers to surf fits (1,2=amap;3,4=b)
					# 5,6 unused, 7,8=brmap

double	sys[NPARAM]			# system parameters
real	ccd[NCCD,3]			# CCD geometry
double	a3[3,3]				# grating transform
real	xics, yics			# pixel values in ICS
real	scaling
int	stat
int	chip

int	m
int	nchp[4]				# array for number of slits on chip
int	ns				# number of slits per chip
real	w0

int	j, nx, ny
int	bluechip, redchip

# int	fits_ize()
int	qxfm(), qxfm_init()

bool	clgetb()
int	clgeti()
real	clgetr()
pointer	immap()
begin
	call clgstr ("image", fname, SZ_FNAME)
	call sprintf (ximage, SZ_FNAME, "%s[%d]")
		call pargstr (fname); call pargi (0)
	im = immap (ximage, READ_ONLY, 0)

# get slider, grating angle, etc and maps
	call get_system_info (im, sys, fda, fdb, fdc)

	ichip1 = clgeti ("chip1")
	ichip2 = clgeti ("chip2")
	regions_only = clgetb ("regions_only")

	call get_slit_info (fname, sdat, nslits, guiname)

# Not all slits are in correct orientation; fix
	call rect_slit_info (sdat, nslits)

# remap slitnumbers:
	do i = 0, nslits-1 {
		SLITNO(sdat,i) = i
	}

# Now get qmodel values (most others are hidden in get_system_info):
	MU(sys) = MU(sys) + DEGTORAD(clgetr ("dmu"))
	GR_YERR(sys) = GR_YERR(sys) + DEGTORAD(clgetr ("droll3"))
	GR_ZERR(sys) = GR_ZERR(sys) + DEGTORAD(clgetr ("do3"))
	scaling = 1. + clgetr ("qmodel.scale_adj")

# Initialize the mappings
	stat = qxfm_init (fda, fdb, map, a3, sys, ccd)
	call gs_ingest (fdc, map[7], map[8])		# reverse beta mapping
	call close (fda)
	call close (fdb)
	call close (fdc)

# CHECK!: USED?
	w0 = 1.e-4 * clgetr ("wave")

#
# TEST QLSOLVE:
	if (clgetb ("tweak"))
		call qltest (fname, ccd, sys, map)
#	call fatal (0, "INTENDED STOP: **QLTEST**")

#
# Write out the "regions" file mapping
#
	call mk_regions_file (map, a3, sys, ccd, sdat, nslits, scaling)
	if (regions_only)
		call fatal (0, "Intended Stop")


# Count chip assignments, based on where object lands
	call amovki (0, nchp, 4)
	do i = 0, nslits-1 {
# calculate the mapping
		stat  = qxfm (map, a3, sys, ccd, XMM(sdat,i), YMM(sdat,i), w0,
		scaling, xics, yics, XPX(sdat,i), YPX(sdat,i), chip, YES, NO)

#		stat1 = qxfm (map, a3, sys, ccd, X1MM(sdat,i), Y1MM(sdat,i), w0,
#		scaling, xics, yics, X1PX(sdat,i), Y1PX(sdat,i), chip, YES, NO)

#		stat2 = qxfm (map, a3, sys, ccd, X2MM(sdat,i), Y2MM(sdat,i), w0,
#		scaling, xics, yics, X2PX(sdat,i), Y2PX(sdat,i), chip, YES, NO)

# Assign chip numbers to the slits; keep count
		if (chip < 1 || chip > 8) {
			call eprintf ("Trouble: chip=%d\n")
				call pargi (chip)
			OBJCHIP[sdat,i] = 0
		} else {
			m = chip
			if (chip > 4)
				m = m - 4
			OBJCHIP[sdat,i] = m
			nchp[m] = nchp[m] + 1
		}
	}

call eprintf ("Number of slits per chip: %d %d %d %d\n")
call pargi (nchp[1])
call pargi (nchp[2])
call pargi (nchp[3])
call pargi (nchp[4])

	imo = immap ("aout.fits[0/8]", NEW_COPY, im)
	IM_PIXTYPE(imo) = TY_REAL
	call imdelf (imo, "BZERO")
	call imdelf (imo, "BSCALE")
	call imunmap (imo)

# FOR EACH CHIP PAIR ...	# XXX probably need a switch for turning off
# certain parts ...
	do chip = ichip1, ichip2 {
	    bluechip = chip
	    redchip = chip + 4

	    call sprintf (ximage, SZ_FNAME, "%s[%d]")
		call pargstr (fname); call pargi (bluechip)
	    call sprintf (oimage, SZ_FNAME, "%s[%d/8]")
		call pargstr ("aout.fits"); call pargi ((chip-1)*2+1)
call eprintf ("%s --> %s\n"); call pargstr (ximage); call pargstr (oimage)

	    imx = immap (ximage, READ_ONLY, 0)
	    imo = immap (oimage, NEW_COPY, imx)
	    IM_PIXTYPE(imo) = TY_REAL
	    call imdelf (imo, "BZERO")
	    call imdelf (imo, "BSCALE")
	    nx = IM_LEN(imo,1)
	    ny = IM_LEN(imo,2)

	    ns = nchp[chip]

	    call chp_alloc (cdat, ns)

	    j = 0
	    do i = 0, nslits-1 {
		if (OBJCHIP(sdat,i) != bluechip)
			next
call eprintf ("DEBUG: assigning %d slit %d\n"); call pargi (j); call pargi (SLITNO(sdat,i))
		SLITNO(cdat,j) = SLITNO(sdat,i)
		j = j + 1
	    }
# Need slits to be in order:
	    call xsort_slits (sdat, cdat, ns)

# Map edges of slits on image
		call map_edges (chip, imx, sdat, cdat, ns, map, a3, sys, ccd, scaling)
# Fit each slit on mask, fit the input angles for speed
		call fit_slit_inp_ang (sdat, cdat, ns, map, a3, sys, ccd, scaling)

		call do_something (im, imx, imo, chip, sdat, cdat, ns, map, a3, sys, ccd, scaling)
	
call eprintf ("DEBUG: unmapping imo ...\n")
	    call imunmap (imo)
call eprintf ("DEBUG: unmapping imx ...\n")
	    call imunmap (imx)

# Do the red chip now:
	    redchip = chip + 4
	    call sprintf (ximage, SZ_FNAME, "%s[%d]")
		call pargstr (fname); call pargi (redchip)	# XXX
	    call sprintf (oimage, SZ_FNAME, "%s[%d/8]")
		call pargstr ("aout.fits"); call pargi (chip*2)	# XXX
call eprintf ("%s --> %s\n"); call pargstr (ximage); call pargstr (oimage)

	    imx = immap (ximage, READ_ONLY, 0)
	    imo = immap (oimage, NEW_COPY, imx)
	    IM_PIXTYPE(imo) = TY_REAL
	    call imdelf (imo, "BZERO")
	    call imdelf (imo, "BSCALE")
	    nx = IM_LEN(imo,1)
	    ny = IM_LEN(imo,2)

		call map_edges (redchip, imx, sdat, cdat, ns, map, a3, sys, ccd, scaling)

		call do_something (im, imx, imo, redchip, sdat, cdat, ns, map, a3, sys, ccd, scaling)
	
call eprintf ("DEBUG: unmapping imo ...\n")
	    call imunmap (imo)
call eprintf ("DEBUG: unmapping imx ...\n")
	    call imunmap (imx)

# End of red chip section


call eprintf ("DEBUG: freeing the fit structures ...\n")
	    do i = 0, ns-1 {
#		call cvfree (CVGAM[cdat,i])
#		call cvfree (CVALP[cdat,i])
		call cvfree (CVS[cdat,i])
		call cvfree (CVO[cdat,i])
		call cvfree (CV2[cdat,i])
		call cvfree (CV1[cdat,i])
	    }
call eprintf ("DEBUG: freeing the slit structures ...\n")
	    call chp_free (cdat, ns)
call eprintf ("DEBUG: done!\n")
	}
	
	
	call imunmap (im)

end


#
# GET_SYSTEM_INFO: get info about grating/slider/etc and maps
#

define	SZ_KWVAL	72

procedure	get_system_info (im, sys, fda, fdb, fdc)

pointer	im				# image header pointer
double	sys[NPARAM]			# system parameters
pointer	fda, fdb, fdc			# pointers to mappings

char	amap[SZ_FNAME], bmap[SZ_FNAME]	# input mappings
char	brmap[SZ_FNAME]			# input reverse mappings
pointer	fd

char	gname[SZ_KWVAL]			# grating name
char	tname[SZ_KWVAL]			# test name
int	slider				# slider number
real	gtilt				# grating tilt
int	tslider	
int	got_info
real	a0, a1, o2, o3, gmm

bool	streq()
real	imgetr()
int	imgeti(), imaccf()
int	fscan(), nscan()
pointer	open()

begin
	if (imaccf (im, "GRATEPOS") == YES) {
		slider = imgeti (im, "GRATEPOS")
		if (slider < 2 || slider > 4)
			call fatal (0, "Illegal GRATEPOS value")
	} else {
		call fatal (0, "NO GRATING KEYWORD")
	}

	if (slider == 3) {
	    if (imaccf (im, "G3TLTVAL") == YES)
		gtilt = imgetr (im, "G3TLTVAL")
	    else
		call fatal (0, "NO G3TLTVAL")
	} else if (slider == 4) {
	    if (imaccf (im, "G4TLTVAL") == YES)
		gtilt = imgetr (im, "G4TLTVAL")
	    else
		call fatal (0, "NO G4TLTVAL")
	} else {
		gtilt = 0.
	}

	if (imaccf (im, "GRATENAM") == YES) {
		call imgstr (im, "GRATENAM", gname, SZ_KWVAL)
	} else {
		call fatal (0, "NO GRATENAM KEYWORD")
	}

# Compare to table; get gmm, roll, o3, etc
	got_info = NO
	fd = open ("deimos$mappings/grating_table.txt", READ_ONLY, TEXT_FILE)
	while (fscan (fd) != EOF) {
		call gargwrd (tname, SZ_KWVAL)
		if (streq (tname, "") || streq (tname, "#"))
			next
		call gargi (tslider)
		if (streq (gname, tname) && slider == tslider) {
			call gargr (gmm)
			call gargr (o2)
			call gargr (o3)
			call gargr (a1)
			call gargr (a0)
			if (nscan() < 7)
				next
			got_info = YES
		}
	}
	if (got_info == NO)
		call fatal (0, "No table info for grating")

	MU(sys) = DEGTORAD(a1 * gtilt + a0)
	GR_YERR(sys) = DEGTORAD(o2)
	GR_ZERR(sys) = DEGTORAD(o3)
	GRLINES(sys) = 1.e-3 * gmm

	if (MU(sys) < -18.d0 || slider == 2)
		ORDER(sys) = 0
	else
		ORDER(sys) = 1

# Get maps
	if (slider == 2) {
		call clgstr ("deimos.amf2map", amap, SZ_FNAME)
		call clgstr ("deimos.bmf2map", bmap, SZ_FNAME)
		call clgstr ("deimos.bmr2map", brmap, SZ_FNAME)
	} else if (slider == 3) {
		call clgstr ("deimos.amf3map", amap, SZ_FNAME)
		call clgstr ("deimos.bmf3map", bmap, SZ_FNAME)
		call clgstr ("deimos.bmr3map", brmap, SZ_FNAME)
	} else if (slider == 4) {
		call clgstr ("deimos.amf4map", amap, SZ_FNAME)
		call clgstr ("deimos.bmf4map", bmap, SZ_FNAME)
		call clgstr ("deimos.bmr4map", brmap, SZ_FNAME)
	}

	fda = open (amap, READ_ONLY, TEXT_FILE)
	fdb = open (bmap, READ_ONLY, TEXT_FILE)
	fdc = open (brmap, READ_ONLY, TEXT_FILE)
	
	call eprintf ("Got: slider=%d grating=%s gmm=%7.2f mu=%7.3f yerr=%6.3f zerr=%6.3f\n")
		call pargi (slider)
		call pargstr (gname)
		call pargr (gmm)
		call pargd (RADTODEG(MU(sys)))
		call pargd (RADTODEG(GR_YERR(sys)))
		call pargd (RADTODEG(GR_ZERR(sys)))

	call eprintf ("Got amap =  '%s'\n")
		call pargstr (amap)
	call eprintf ("Got bmap =  '%s'\n")
		call pargstr (bmap)
	call eprintf ("Got brmap = '%s'\n")
		call pargstr (brmap)

end

procedure	mk_regions_file (map, a3, sys, ccd, sdat, nslits, scaling)

pointer	map[8]				# pointers to surf fits (1,2=amap;3,4=b)
					# 5,6 unused, 7,8=brmap
double	sys[NPARAM]			# system parameters
real	ccd[NCCD,3]			# CCD geometry
double	a3[3,3]				# grating transform
pointer	sdat
int	nslits
real	scaling
pointer	fdg

int	i, j
real	xobj, xs1, xs2
real	yobj, ys1, ys2
real	xics, yics			# pixel values in ICS
int	stat, stat1, stat2
int	chip
real	xfudge
real	w0
real	wave[NWAV]		# XXX TMP!

int	qxfm()

real	clgetr()
pointer	open()
begin
	fdg = open ("aout.txt", APPEND, TEXT_FILE)

	xfudge = clgetr ("xreg_off")
# XXX TMP!
wave[1] = 5577.34e-4
wave[2] = 6948.94e-4
wave[3] = 7316.28e-4
wave[4] = 8430.17e-4
wave[5] = 8885.83e-4
wave[6] = 7600.00e-4

#	call fprintf (fdg, "global wcs=wcsp\n")
	do j = 1, NWAV {
	    call fprintf (fdg, "global color=%s\n")
	    switch (j) {
		case 1:
			call pargstr ("cyan")
		case 2:
			call pargstr ("green")
		case 3:
			call pargstr ("yellow")
		case 4:
			call pargstr ("red")
		case 5:
			call pargstr ("magenta")
		case 6:
			call pargstr ("blue")
	    }
	    w0 = wave[j]

	    do i = 0, nslits-1 {
# calculate the mapping
		stat  = qxfm (map, a3, sys, ccd, XMM(sdat,i), YMM(sdat,i), w0,
		scaling, xics, yics, xobj, yobj, chip, YES, YES)

		stat1 = qxfm (map, a3, sys, ccd, X1MM(sdat,i), Y1MM(sdat,i), w0,
		scaling, xics, yics, xs1, ys1, chip, NO, YES)

		stat2 = qxfm (map, a3, sys, ccd, X2MM(sdat,i), Y2MM(sdat,i), w0,
		scaling, xics, yics, xs2, ys2, chip, NO, YES)

		if (yobj > 0. && yobj < 8193.) {
			call fprintf (fdg, "line %7.1f %7.1f  %7.1f %7.1f # line=0 0 text={%s}\n")
				call pargr (xs1+xfudge)
				call pargr (ys1)
				call pargr (xs2+xfudge)
				call pargr (ys2)
				call pargstr (OBJNAME(sdat,i))

			call fprintf (fdg, "cross point %7.1f %7.1f\n")
				call pargr (xobj+xfudge)
				call pargr (yobj)
		}
	    }
	}
	call close (fdg)
end

#
# GET_SLIT_INFO: gets info on alignment stars from Mask Design FITS tables.
# Generic (I think);  the only hardcode is that the XMM,YMM must be in order
# following XMM1.

procedure	get_slit_info (fname, sdat, nbox, guiname)

char	fname[ARB]		# FITS file name w/o extn
pointer	sdat			# box data struct
int	nbox			# number of boxes
char	guiname[ARB]		# GI Name for slitmask

char	xname[SZ_FNAME]		# FITS name with extention
char	xline[80]		# XXX hardcode
%	character*80 f77nam
%	character*80 f77lin
%	character*80 f77nul

int	lu
int	stat
int	nr

int	i, j, k
int	n, ndx
int	nslit
int	nobj
int	ic1, ic2, ic3, ic4
real	xa, xb, ya, yb
real	xf

bool	isnul			# Null value(s) present?
int	inul			# dummy null value
real	rnul			# dummy null value

int	posn_hdu(), get_col_ndx()
#bool	strne()
begin
	rnul = INDEF
	call f77pak ("???", f77nul, 80)

	stat = 0		# apparently needed for FITSIO

	call f77pak (fname, f77nam, SZ_FNAME)

	call ftgiou (lu, stat)
	call ftnopn (lu, f77nam, READONLY, stat)

# Get DesiSlits table (3) for slitTyp=A (col 5)

	if (posn_hdu (lu, "DesiSlits") != OK)
		call fatal (0, "")

	if (get_col_ndx (lu, "slitTyp", ic1) != OK)
		call fatal (0, "column must exist!")
	if (get_col_ndx (lu, "dSlitId", ic2) != OK)
		call fatal (0, "column must exist!")

	call ftgnrw (lu, nr, stat)
	nslit = nr

	call rdx_alloc (sdat, nslit)

	ndx = 0
	do i = 0, nslit-1 {
		j = i + 1
		call ftgcvs (lu, ic1, j, 1, 1, f77nul, f77lin, isnul, stat)
		call f77upk (f77lin, xline, 1)
#		if (strne (xline, "A") && strne (xline, "P"))
#			next

		call ftgcvj (lu, ic2, j, 1, 1, inul, SLITNO(sdat,ndx), isnul, stat)

# XXX ... and box size (col 6,8)?

# increment counter
		ndx = ndx + 1
	}
	nbox = ndx


# Get BluSlits table; search dSlitId, get average of xmm, ymm
# NB: the order in the table is important here; partially checked

	if (posn_hdu (lu, "BluSlits") != OK)
		call fatal (0, "")

	if (get_col_ndx (lu, "dSlitId", ic1) != OK)
		call fatal (0, "column must exist!")
	if (get_col_ndx (lu, "slitX1", ic2) != OK)
		call fatal (0, "column must exist!")
	if (get_col_ndx (lu, "slitY4", ic4) != OK)
		call fatal (0, "column must exist!")

# check format:
	if (ic4 != ic2+7 || nr != nslit) {
		call eprintf ("%s, nr=%d stat=%d\n")
		call pargstr (xname)
		call pargi (nr)
		call pargi (stat)
		call fatal (0, "Unexpected format in BluSlits table")
	}

	do i = 0, nbox-1 {
		ndx = SLITNO(sdat,i)
		do k = 0, nslit-1 {
			j = k + 1
			call ftgcvj (lu, ic1, j, 1, 1, inul, n, isnul, stat)
			if (n != ndx)
				next

			call ftgcve (lu, ic2,   j, 1, 1, rnul, xa, isnul, stat)
			call ftgcve (lu, ic2+6, j, 1, 1, rnul, xb, isnul, stat)
			X2MM(sdat,i) = 0.5 * (xa + xb)

			call ftgcve (lu, ic2+1, j, 1, 1, rnul, ya, isnul, stat)
			call ftgcve (lu, ic2+7, j, 1, 1, rnul, yb, isnul, stat)
			Y2MM(sdat,i) = 0.5 * (ya + yb)

			call ftgcve (lu, ic2+2, j, 1, 1, rnul, xa, isnul, stat)
			call ftgcve (lu, ic2+4, j, 1, 1, rnul, xb, isnul, stat)
			X1MM(sdat,i) = 0.5 * (xa + xb)

			call ftgcve (lu, ic2+3, j, 1, 1, rnul, ya, isnul, stat)
			call ftgcve (lu, ic2+5, j, 1, 1, rnul, yb, isnul, stat)
			Y1MM(sdat,i) = 0.5 * (ya + yb)
		}
	}
			


# Get SlitObjMap table; search on dSlitId to get ObjectId

	if (posn_hdu (lu, "SlitObjMap") != OK)
		call fatal (0, "")

	if (get_col_ndx (lu, "dSlitId", ic1) != OK)
		call fatal (0, "column must exist!")
	if (get_col_ndx (lu, "ObjectId", ic2) != OK)
		call fatal (0, "column must exist!")
	if (get_col_ndx (lu, "TopDist", ic3) != OK)
		call fatal (0, "column must exist!")
	if (get_col_ndx (lu, "BotDist", ic4) != OK)
		call fatal (0, "column must exist!")

	call ftgnrw (lu, nr, stat)
	nobj = nr
	do i = 0, nbox-1 {
		ndx = SLITNO(sdat,i)
		do k = 0, nobj-1 {
			j = k + 1
			call ftgcvj (lu, ic1, j, 1, 1, inul, n, isnul, stat)
			if (n != ndx)
				next

			call ftgcvj (lu, ic2, j, 1, 1, inul, OBJNO(sdat,i), isnul, stat)

# work out the object location, based on fractions along slit
# NB THIS NEEDs CHECKING, PLUS ERROR IN SENSE ...   XXX
			call ftgcve (lu, ic3, j, 1, 1, 0., xa, isnul, stat)
			call ftgcve (lu, ic4, j, 1, 1, 0., xb, isnul, stat)
			xf = xa / (xa + xb)
			XMM(sdat,i) = X2MM(sdat,i) * xf + X1MM(sdat,i) * (1.-xf)
			YMM(sdat,i) = Y2MM(sdat,i) * xf + Y1MM(sdat,i) * (1.-xf)
			break
		}
	}



# Get ObjectCat table (1); search ObjectId (Col 1), get OBJECT (col 2), etc

	if (posn_hdu (lu, "ObjectCat") != OK)
		call fatal (0, "")

### WORK HERE
	if (get_col_ndx (lu, "ObjectId", ic1) != OK)
		call fatal (0, "column must exist!")
	if (get_col_ndx (lu, "OBJECT", ic2) != OK)
		call fatal (0, "column must exist!")
	if (get_col_ndx (lu, "mag", ic3) != OK)
		call fatal (0, "column must exist!")
	if (get_col_ndx (lu, "pBand", ic4) != OK)
		call fatal (0, "column must exist!")

	do i = 0, nbox-1 {
		ndx = OBJNO(sdat,i)
		do k = 0, nobj-1 {
			j = k + 1
			call ftgcvj (lu, ic1, j, 1, 1, inul, n, isnul, stat)
			if (n != ndx)
				next

			call ftgcvs (lu, ic2, j, 1, 1, f77nul, f77lin, isnul, stat)
			call f77upk (f77lin, OBJNAME(sdat,i), SZ_ID-1)
			call ftgcve (lu, ic3, j, 1, 1, -99., MAG(sdat,i), isnul, stat)
			call ftgcvs (lu, ic4, j, 1, 1, f77nul, f77lin, isnul, stat)
			call f77upk (f77lin, xline, 1)
			PBAND(sdat,i) = xline[1]
			if (xline[2] != EOS)
				call eprintf ("Warning; pband truncated\n")
			break
		}
	}


# Get MaskBlu table for guiname

	if (posn_hdu (lu, "MaskBlu") != OK)
		call fatal (0, "")

	if (get_col_ndx (lu, "guiname", ic1) != OK)
		call fatal (0, "column must exist!")

	call ftgcvs (lu, ic1, 1, 1, 1, f77nul, f77lin, isnul, stat)
	if (isnul)
		call strcpy ("", guiname, SZ_GUINAM)
	else
		call f77upk (f77lin, guiname, SZ_GUINAM)


	call ftclos (lu, stat)


end

#
# RECT-SLIT_INFO: Rectify inconsistencies in slit into -- eg x1<x2
# Eventually, top and bottom will go here, too ...
#

procedure	rect_slit_info (sdat, nslits)

pointer	sdat
int	nslits

int	i
real	hold
begin
	do i = 0, nslits-1 {
		if (X2MM(sdat,i) < X1MM(sdat,i)) {
			hold = X2MM(sdat,i)
			X2MM(sdat,i) = X1MM(sdat,i)
			X1MM(sdat,i) = hold
			hold = Y2MM(sdat,i)
			Y2MM(sdat,i) = Y1MM(sdat,i)
			Y1MM(sdat,i) = hold
		}
	}
end

#
# RDX_ALLOC: allocate arrays for Redux (stolen from targ_init)
#
procedure	rdx_alloc (sdat, nbox)

pointer	sdat
int	nbox

begin
# Allocate the vectors
	call malloc (sdat, RDATLEN, TY_STRUCT)
	call malloc (PTSLTNO(sdat), nbox, TY_INT)
	call malloc (PTOBJNO(sdat), nbox, TY_INT)
	call malloc (PTMAG(sdat), nbox, TY_REAL)
	call malloc (PTPBAND(sdat), nbox, TY_CHAR)	## XXX should be general
	call malloc (PTXMM(sdat), nbox, TY_REAL)
	call malloc (PTYMM(sdat), nbox, TY_REAL)
	call malloc (PTX1MM(sdat), nbox, TY_REAL)
	call malloc (PTY1MM(sdat), nbox, TY_REAL)
	call malloc (PTX2MM(sdat), nbox, TY_REAL)
	call malloc (PTY2MM(sdat), nbox, TY_REAL)
	call malloc (PTXPX(sdat), nbox, TY_REAL)
	call malloc (PTYPX(sdat), nbox, TY_REAL)
	call malloc (PTX1PX(sdat), nbox, TY_REAL)
	call malloc (PTY1PX(sdat), nbox, TY_REAL)
	call malloc (PTX2PX(sdat), nbox, TY_REAL)
	call malloc (PTY2PX(sdat), nbox, TY_REAL)
	call malloc (PTCHIP(sdat), nbox, TY_INT)

# Allocate the string name array
	call malloc (PTNAME(sdat), nbox*SZ_ID, TY_CHAR)

end

#
# CHP_ALLOC: allocate arrays for CHIP Redux (stolen from targ_init)
#
procedure	chp_alloc (cdat, nslit)

pointer	cdat
int	nslit

begin
# Allocate the vectors
	call malloc (cdat, CDATLEN, TY_STRUCT)
	call malloc (PTSLTNO(cdat), nslit, TY_INT)
	call calloc (PTOVERL(cdat), nslit, TY_INT)	# Note zeroing
	call malloc (PTCV1(cdat), nslit, TY_STRUCT)
	call malloc (PTCV2(cdat), nslit, TY_STRUCT)
	call malloc (PTCVO(cdat), nslit, TY_STRUCT)
	call malloc (PTCVS(cdat), nslit, TY_STRUCT)
	call malloc (PTCVA(cdat), nslit, TY_STRUCT)
	call malloc (PTCVG(cdat), nslit, TY_STRUCT)

	return

entry	chp_free (cdat, nslit)

	call mfree (PTCVG(cdat), TY_STRUCT)
	call mfree (PTCVA(cdat), TY_STRUCT)
	call mfree (PTCVS(cdat), TY_STRUCT)
	call mfree (PTCVO(cdat), TY_STRUCT)
	call mfree (PTCV2(cdat), TY_STRUCT)
	call mfree (PTCV1(cdat), TY_STRUCT)
	call mfree (PTOVERL(cdat), TY_INT)
	call mfree (PTSLTNO(cdat), TY_INT)

end


# GET_TBL_NDX: get index of a particular table
# Currently not used

procedure	get_tbl_ndx (fname, desname, open_it, lu, ndx)

char	fname[ARB]		# name of image
char	desname[ARB]		# desired table name
int	open_it			# open the file?
int	lu			# returned (open) unit number for FITSIO
int	ndx			# returned 0-indexed number of HDU

int	stat

%	character*80 f77nam
%	character*80 extnam

int	access()
begin

	ndx = -1
	stat = 0	# need for FITSIO, which doesn't clean this out

	if (access (fname, 0, 0) == NO) {
		call eprintf ("File does not exist (%-s)\n")
			call pargstr (fname)
		return
	}

	call f77pak (fname, f77nam, SZ_FNAME)
	call f77pak (desname, extnam, SZ_FNAME)
	
	if (open_it == YES) {
		call ftgiou (lu, stat)
		call ftnopn (lu, f77nam, READONLY, stat)
	}

	call ftmnhd (lu, ANY_HDU, extnam, 0, stat)

	if (stat == BAD_HDU_NUM) {
		call eprintf ("Extension does not exist (%-s)\n")
			call pargstr (desname)
		return
	}

	call ftghdn (lu, ndx)
	ndx = ndx - 1
end

#
# XSORT SLITS: arrange slit numbers for increasing x-order
#

procedure	xsort_slits (sdat, cdat, ns)

pointer	sdat				# pointer to slits struct
pointer	cdat				# pointer to chip_slits struct
int	ns

int	ndx
int	i
pointer	sp, bufx

begin
	call smark (sp)
	call salloc (bufx, ns, TY_REAL)

	do i = 0, ns-1 {
		ndx = SLITNO (cdat,i)
		Memr[bufx+i] = X1MM(sdat,ndx)		# could sort on ...?
	}

	call ndxsortr (Memr[bufx], SLITNO(cdat,0), ns)


	call sfree (sp)

end

#
# NDXSORTR: Sort a list of real values (low to high) and carry index
#

procedure	ndxsortr (r, ndx, n)

real	r[n]				# Real vector to sort
int	ndx[n]				# index vector
int	n				# number of points

int	i, j
real	h
int	ih

begin
# Sort the list (low-to-high)
	do i = 1, n-1 {
	    do j = 1, n-i {
		if (r[j] > r[j+1]) {
			h = r[j+1]
			r[j+1] = r[j]
			r[j] = h

			ih = ndx[j+1]
			ndx[j+1] = ndx[j]
			ndx[j] = ih
		}
	    }
	}
end


procedure	do_something (im0, im1, im2, chip, sdat, cdat, ns, map, a3, sys, ccd, scaling)

pointer	im0, im1, im2		# input[0], input, output image pointers
int	chip
pointer	sdat
pointer	cdat
int	ns

pointer	map[8]				# pointers to surf fits (1,2=amap;3,4=b)
double	sys[NPARAM]			# system parameters
real	ccd[NCCD,3]			# CCD geometry
double	a3[3,3]				# grating transform
real	scaling

int	i, j
int	nx, ny
int	nx2
int	precol, postcol

int	i1, i2, j1, j2, npix
int	off
int	stat
real	bias
pointer	bufline
pointer	bufchip, bufwave, bufsky

char	datasec[80]

int	get_datasec()
real	vsum1()

# XXX totally tmp
int	n1, n2, clgeti()
real	cveval()

int	imgeti()
real	imgetr()
pointer	imgl2r(), impl2r(), imgs2r()

begin
# XXX TMP:
	n1 = clgeti ("n1")
	n2 = clgeti ("n2")

	nx = IM_LEN(im1,1)
	ny = IM_LEN(im1,2)
	precol  = imgeti (im0, "PRECOL")
	postcol = imgeti (im0, "POSTPIX")

	IM_LEN(im2,1) = nx - precol - postcol
	nx2 = IM_LEN(im2,1)

# This is very cheesy -- we want the full-up solution of mos_init()

	if (get_datasec (im1, i1, i2, j1, j2) != OK)
		call fatal (0, "DATASEC failed")
	call sprintf (datasec, 80, "[%d:%d,%d:%d]")
		call pargi (i1-precol)
		call pargi (i2-precol)
		call pargi (j1)
		call pargi (j2)
	call impstr (im2, "DATASEC", datasec)
	call imputr (im2, "CRPIX1P", imgetr (im1, "CRPIX1P")-precol)
	call imputr (im2, "CRPIX1",  imgetr (im1, "CRPIX1")-precol)
	

# Cheap bias subtract XXX
	i1 = nx - postcol + 1
	i2 = nx
	npix = (i2 - i1 + 1) * ny
	bias = vsum1 (Memr[imgs2r (im1, i1, i2, 1, ny)], npix) / npix
	bias = -bias

# Also pretty cheesy: Allocate space for an array for the chip:

	call malloc (bufchip, nx2*ny, TY_REAL)
	call calloc (bufwave, nx2*ny, TY_REAL)
	call calloc (bufsky,  nx2*ny, TY_REAL)

	if (chip > 4) {
	    do j = 1, ny {
		bufline = imgl2r(im1, j)
		off = (4096-j) * nx2 + nx2 - 1
		do i = 0, nx2-1 {
			Memr[bufchip+off-i] = Memr[bufline+precol+i] + bias
		}

#if (mod (j, 400) == 0) {
#call eprintf ("ACCUMSKY %d\n"); call pargi (j)
#}
	    }
	} else {
	    do j = 1, ny {
		bufline = imgl2r(im1, j)
		off = (j-1) * nx2
		call aaddkr (Memr[bufline+precol], bias, Memr[bufchip+off], nx2)

#if (mod (j, 400) == 0) {
#call eprintf ("ACCUMSKY %d\n"); call pargi (j)
#}
	    }
	}

# solve the sky:
	do i = 0, ns-1 {
		call fit_slit_sky (Memr[bufchip], Memr[bufwave], Memr[bufsky],
				nx2, ny, map, a3, sys, ccd, scaling,
				chip, sdat, cdat, i, ns, stat, n1, n2)
		if (stat != OK) {
			call eprintf ("Sky fit stat=%d (x= %5f %5f)\n")
				call pargi (stat)
				call pargr (cveval (CV1(cdat,i), 2048.))
				call pargr (cveval (CV2(cdat,i), 2048.))
		}
	}

# remove sky and write out image:
call eprintf ("DEBUG: writing out ... ns=%d\n"); call pargi (ns)
	if (chip > 4) {
	    do j = 1, ny {
		bufline = impl2r (im2, j)
		off = (4096-j) * nx2 + nx2 - 1
		do i = 0, nx2-1 {
		    Memr[bufline+i] = Memr[bufchip+off-i] - Memr[bufsky+off-i]
		}
	    }
	} else {
	    do j = 1, ny {
		off = (j-1) * nx2
		call asubr (Memr[bufchip+off], Memr[bufsky+off], Memr[impl2r (im2, j)], nx2)
	    }
	}

call eprintf ("DEBUG: Finished chip %d\n")
call pargi (chip)
	call mfree (bufsky, TY_REAL)
	call mfree (bufwave, TY_REAL)
	call mfree (bufchip, TY_REAL)
#call eprintf ("DEBUG: memory free\n")

	

end

#
# MAP_EDGES: map slit edges onto pixels; check for overlap and resolve these
# cases.
#

procedure	map_edges (chip, im, sdat, cdat, ns, map, a3, sys, ccd, scaling)

int	chip
pointer	im
pointer	sdat
pointer	cdat
int	ns

pointer	map[8]				# pointers to surf fits (1,2=amap;3,4=b)
double	sys[NPARAM]			# system parameters
real	ccd[NCCD,3]			# CCD geometry
double	a3[3,3]				# grating transform
real	scaling

int	i, j, ndx
int	nx, ny
real	wt

int	stat
real	xpix, ypix
real	x1, y1
real	x2, y2
real	y
real	wmin, wmax, winc, w1, w2
real	w
real	xics, yics			# pixel values in ICS
int	nordedg				# order of polynomial for edge fitting
real	yemn, yemx			# min, max value for edge fitting

int	qxfm()
real	qwave()

real	cveval()
begin
	nordedg = 4			# XXX OK??
	yemn = -400.			# XXX Hardcode
	yemx = 4400.			# XXX Hardcode
	do i = 0, ns-1 {
		call cvinit (CV1[cdat,i], LEGENDRE, nordedg, yemn, yemx)
		call cvinit (CV2[cdat,i], LEGENDRE, nordedg, yemn, yemx)
		call cvinit (CVO[cdat,i], LEGENDRE, nordedg, yemn, yemx)
	}

	wt = 1.
	wmin = 0.6000			# XXX WHOA!
	wmax = 0.9600			# XXX WHOA!
#	wmin = 0.2400			# XXX WHOA!
#	wmax = 0.6400			# XXX WHOA!
	winc = 0.0100
call eprintf ("REMINDER -- need WAVE RANGE CALC'd\n")
	
	nx = IM_LEN(im,1)
	ny = IM_LEN(im,2)

	do i = 0, ns-1 {
		ndx = SLITNO(cdat,i)
# Work out wavelength range (mask_to_pre, ccd_to_post); some iteration needed?
##### WORK
		xpix = XPX(sdat,ndx)
		ypix = 0.5
		call ccd_to_ics (xpix, ypix, ccd, chip, xics, yics)
		w1 = qwave (map, a3, sys, ccd, XMM(sdat,ndx), YMM(sdat,ndx),
							scaling, xics, yics)
		ypix = 4096.5
		call ccd_to_ics (xpix, ypix, ccd, chip, xics, yics)
		w2 = qwave (map, a3, sys, ccd, XMM(sdat,ndx), YMM(sdat,ndx),
							scaling, xics, yics)

		wmin = w1 - 0.1 * (w2 - w1)
		wmax = w2 + 0.1 * (w2 - w1)
		winc = (wmax-wmin) / NWSTEP
call eprintf ("DEBUG: slit edges fit on range %6f--%6f, steps %6f\n")
call pargr (wmin)
call pargr (wmax)
call pargr (winc)
##### WORK

# Work out fit of x(y)
		do j = 0, NWSTEP {
			w = winc * j + wmin
			stat  = qxfm (map, a3, sys, ccd, X1MM(sdat,ndx), Y1MM(sdat,ndx), w, scaling, xics, yics, x1, y1, chip, NO, NO)

# If well off chip, keep going
			if (y1 < yemn || y1 > yemx)	# XXX
				next

			call cvaccum (CV1(cdat,i), y1, x1, wt, WTS_USER)

			stat  = qxfm (map, a3, sys, ccd, X2MM(sdat,ndx), Y2MM(sdat,ndx), w, scaling, xics, yics, x2, y2, chip, NO, NO)

			if (y2 < yemn || y2 > yemx)	# XXX
				next

			call cvaccum (CV2(cdat,i), y2, x2, wt, WTS_USER)
		}

		call cvsolve (CV1(cdat,i), stat)
		if (stat != OK) {
			call eprintf ("CVSOLVE1: stat=%d\n")
				call pargi (stat)
		}
		call cvsolve (CV2(cdat,i), stat)
		if (stat != OK) {
			call eprintf ("CVSOLVE2: stat=%d\n")
				call pargi (stat)
		}

		if (i > 0) {
		    do j = 1, 4096, 128  {
			y = j 
			x1 = cveval (CV1(cdat,i), y)
			x2 = cveval (CV2(cdat,i-1), y)
			call cvaccum (CVO(cdat,i), y, 0.5*(x2+x1), wt, WTS_USER)
		    }
		    call cvsolve (CVO(cdat,i), stat)
		    if (stat != OK) {
			call eprintf ("CVSOLVEO: stat=%d\n")
				call pargi (stat)
		    }
		}

#do k = 1, 4096, 128 {
#call printf ("tile %d; image; line %7.1f %7.1f  %7.1f %7.1f # line=0 0\n")
#ypix = k
#xpix = cveval (CV1(cdat,i), ypix)
#call pargi (chip)
#call pargr (xpix)
#call pargr (ypix)
#ypix = k + 128
#xpix = cveval (CV1(cdat,i), ypix)
#call pargr (xpix)
#call pargr (ypix)
#}
#do k = 1, 4096, 128 {
#call printf ("tile %d; image; line %7.1f %7.1f  %7.1f %7.1f # line=0 0\n")
#ypix = k
#xpix = cveval (CV2(cdat,i), ypix)
#call pargi (chip)
##call pargr (xpix)
#call pargr (ypix)
#ypix = k + 128
#xpix = cveval (CV2(cdat,i), ypix)
#call pargr (xpix)
#call pargr (ypix)
#}


## Now get the wavelength extremes for each slit and init the sky fits
		ypix = 0.5
		xpix = cveval (CV1(cdat,i), ypix)
		call ccd_to_ics (xpix, ypix, ccd, chip, xics, yics)
		w1 = qwave (map, a3, sys, ccd, X1MM(sdat,ndx), Y1MM(sdat,ndx),
							scaling, xics, yics)
		xpix = cveval (CV2(cdat,i), ypix)
		call ccd_to_ics (xpix, ypix, ccd, chip, xics, yics)
		w2 = qwave (map, a3, sys, ccd, X2MM(sdat,ndx), Y2MM(sdat,ndx),
							scaling, xics, yics)
		wmin = min (w1, w2)	# XXX OK since we've flipped upper chips

		ypix = 4096.5
		xpix = cveval (CV1(cdat,i), ypix)
		call ccd_to_ics (xpix, ypix, ccd, chip, xics, yics)
		w1 = qwave (map, a3, sys, ccd, X1MM(sdat,ndx), Y1MM(sdat,ndx), scaling, xics, yics)
		xpix = cveval (CV2(cdat,i), ypix)
		call ccd_to_ics (xpix, ypix, ccd, chip, xics, yics)
		w2 = qwave (map, a3, sys, ccd, X2MM(sdat,ndx), Y2MM(sdat,ndx), scaling, xics, yics)
		wmax = max (w1, w2)	# XXX OK since we've flipped upper chips

		call cvinit (CVS[cdat,i], SPLINE3, NORDSKY, wmin, wmax)
call eprintf ("Slit %d: lambda range= %5f %5f\n")
call pargi (i)
call pargr(wmin)
call pargr(wmax)

	}
end


bool	procedure	test_stat (i, mask)

int	i
int	mask

begin
	if (and (i, mask) == 0)
		return (true)
	else
		return (false)
end

#
# QWAVE: get a wavelength based on maps, slit and pixels coords
#

real	procedure qwave (map, a3, sys, ccd, xmm, ymm, scaling, xics, yics)

pointer	map[8]			# pointer to maps
double	a3[3,3]			# grating transform
double	sys[NPARAM]		# system parameters
real	ccd[NCCD,3]		# CCD geometry
real	xmm, ymm		# x,y slitmask coords
real	scaling			# scaling adjustment
real	xics, yics		# pixel values in ICS
real	wave			# wavelength (um) returned

double	r[3]
double	alpha, beta, gamma
real	tanx, tany		# should be double; check mapping evals
real	x, y

real	gseval()

begin

# Get mapping and convert to r[3]
	tanx = gseval (map[1], xmm, ymm)
	tany = gseval (map[2], xmm, ymm)
# call eprintf ("tanx,y: %5f %5f\n")
# call pargr (tanx)
# call pargr (tany)

	r[3] = -1.d0 / sqrt (1. + tanx*tanx + tany*tany)
	r[1] = r[3] * tanx
	r[2] = r[3] * tany

# xform into grating system
	call gen_xfm (r, a3, YES)

# convert to alpha,gamma
	alpha = -atan2 (-r[2], -r[3])
	gamma = atan2 (r[1], sqrt (r[3]*r[3]+r[2]*r[2]))

# (we can check the gamma as error)

# Convert pixels to beta, gamma:

	x = xics / scaling
	y = yics / scaling

	tanx = gseval (map[7], x, y)
	tany = gseval (map[8], x, y)

	r[3] = -1.d0 / sqrt (1. + tanx*tanx + tany*tany)
	r[1] = r[3] * tanx
	r[2] = r[3] * tany

# xform into grating system
	call gen_xfm (r, a3, YES)

# convert to beta,gamma
	beta = -atan2 (-r[2], -r[3])
	gamma = atan2 (r[1], sqrt (r[3]*r[3]+r[2]*r[2]))

# XXX at this point, gamma is negative wrt gamma from above.  Since cos(gamma)
# is insensitive, it doesn't matter, but it warrants later resolution!!

# Apply the grating equation

	wave = cos (gamma) / (ORDER(sys)*GRLINES(sys)) * (sin (beta) + sin (alpha))
	return (wave)
end



procedure	accum_sky (line, wave, nx, chip, y, sdat, cdat, ns, map, a3, sys, ccd, scaling) 

real	line[nx], wave[nx]
int	nx
int	chip
real	y
pointer	sdat
pointer	cdat
int	ns

pointer	map[8]				# pointers to surf fits (1,2=amap;3,4=b)
double	sys[NPARAM]			# system parameters
real	ccd[NCCD,3]			# CCD geometry
double	a3[3,3]				# grating transform
real	scaling

int	i, j, i1, i2
int	ndx
real	x1, x2
real	xs1, xs2, ys1, ys2
real	xslit, yslit, xics, yics, f, x

real	qwave()

real	cveval()
begin
# zero out wavelength line
	call amovkr (0., wave, nx)

# go through loop of slits
	do i = 0, ns-1 {
		ndx = SLITNO(cdat,i)
		x1 = cveval (CV1(cdat,i), y)
		x2 = cveval (CV2(cdat,i), y)
		xs1 = X1MM(sdat,ndx)
		xs2 = X2MM(sdat,ndx)
		ys1 = Y1MM(sdat,ndx)
		ys2 = Y2MM(sdat,ndx)

# calculate wavelength
		i1 = x1 + 0.5
		i2 = x2 + 0.5
		i1 = max (i1, 1)
		i2 = min (i2, 2048)
		do j = i1, i2 {
			x = j
			f = (x2 - x) / (x2 - x1)	# XXX technically wrong
							# In addition, see (*)
			xslit =	f * xs1 + (1.-f) * xs2
			yslit = f * ys1 + (1.-f) * ys2
			call ccd_to_ics (x, y, ccd, chip, xics, yics)
			wave[j] = qwave (map, a3, sys, ccd, xslit, yslit, scaling, xics, yics)
		}

	}
#
# (*) There's a problem that slit edge has dual purpose, in defining where we
# work AND in defining where the slit edge falls.  It's normally not a problem,
# but in the case of slitlet overlap it becomes a mess.
#

end

procedure	sub_sky (line, wave, nx, y, cdat, ns)

real	line[nx], wave[nx]
int	nx
real	y
pointer	cdat
int	ns

int	i, k, i1, i2
pointer	cv

real	cveval()
begin
# go through loop of slits
	do k = 0, ns-1 {
		i1 = cveval (CV1(cdat,k), y) + 0.5
		i2 = cveval (CV2(cdat,k), y) + 0.5
		i1 = max (i1, 1)
		i2 = min (i2, 2048)
		cv = CVS(cdat,k)
		do i = i1, i2 {
			line[i] = line[i] - cveval (cv, wave[i])
#			line[i] = cveval (cv, wave[i])
		}
	}
end


#
# FIT_SLIT_SKY: For a given slit, calc wave, collect sky values, fit sky and
# then add to the sky array.  We end with 3 arrays: data, wave, sky
#
# FUTURE UPGRADE: for speed, fit slit tanx, tany as functions of f. The
# evaluation of these functions will be _much_ faster than the high-order
# surface fit evaluation!  This alone should save significant time.


procedure	fit_slit_sky (data, wave, sky, nx, ny, map, a3, sys, ccd, scaling, chip, sdat, cdat, is, ns, stat, n1, n2)

real	data[nx,ny]
real	wave[nx,ny]
real	sky[nx,ny]
int	nx, ny
#
pointer	map[8]				# pointers to surf fits (1,2=amap;3,4=b)
double	sys[NPARAM]			# system parameters
real	ccd[NCCD,3]			# CCD geometry
double	a3[3,3]				# grating transform
real	scaling
#
int	chip
pointer	sdat
pointer	cdat
int	is				# number of THIS slit
int	ns				# number of slits
int	stat
int	n1,n2			# XXX TMP!

int	ndx, nxl, npt
int	i, j, n
int	i1, i2
int	ia1[4096], ia2[4096]		# range of pixels for "all" slit
int	ib1[4096], ib2[4096]		# range of pixels assigned slit
int	ic1[4096], ic2[4096]		# range of pixels for "clean" slit
int	lenbuf
int	off
int	over_low, over_upp, overlap	# overlapping slit flags
real	x1, x2
#real	xs1, xs2, ys1, ys2
real	xa, xb, xc
real	xslit, yslit, xics, yics, f, x
real	alp, gam			# evaluated alpha, gamma
real	y
real	cutoff		# TMP?
real	rnsq, g, nsig, res		# all to do with weighting scheme
pointer	cvalo, cvahi			# curfit edges "all" slit
pointer	cvblo, cvbhi			# curfit edges "all" slit
pointer	cvclo, cvchi			# curfit edges "clean" slit
pointer	cvf
pointer	bufdat, bufwav, bufwts, buffit	# sky fitting vectors
pointer	buflam				# This is array of ALL waves for slit

real	qqwave()

real	cveval()
real	cvstatr()

begin
call eprintf ("Fitting sky for slit %d\n"); call pargi (is)

# Problem: Wave must be calculated over the _entire_ slit, but will only be
# saved in wave array over the _assigned_ region.  There is still a mess here.
# Obviously, we must save a larger array of wavelengths, only part of which
# is stored in the final array and of which only a smaller fraction gets
# stored in the sky-fitting vector

# There are three regions defined:
# (a) the full slit coverage
# (b) the assigned-slit region; overlaps split in 2
# (c) the clean-sky region; devoid of overlap & slit edges
# Obviously, we want to collect sky for fitting from "c".  We must carry
# "a" for building up the sky vector; "b" is saved in the output array for
# future extractions

# Solution: create/store local array of wavelengths; what gets used then
# depends on (ia1,ia2); (ib1,ib2); (ic1,ic2)

	cvalo = CV1(cdat,is) 
	cvahi = CV2(cdat,is) 
	cvblo = cvalo
	cvbhi = cvahi
	cvclo = cvalo
	cvchi = cvahi
	cvf  = CVS(cdat,is)

# Check for overlap:
	y = 2048.		# XXX TEST
	overlap = NO
	over_low = NO
	over_upp = NO
# lower edge ...
	if (is > 0) {
		if (cveval (CV2(cdat,is-1), y) > cveval (cvalo, y)) {
			over_low = YES
			overlap  = YES
call eprintf ("SLIT OVERLAP (low) DETECTED (%6f > %6f)\n")
call pargr (cveval (CV2(cdat,is-1), y)); call pargr (cveval (cvalo, y))
		}
	}
# upper edge ...
	if (is < ns-1) {
		if (cveval (CV1(cdat,is+1), y) < cveval (cvahi, y)) {
			over_upp = YES
			overlap  = YES
call eprintf ("SLIT OVERLAP (upp) DETECTED (%6f < %6f)\n")
call pargr (cveval (CV1(cdat,is+1), y)); call pargr (cveval (cvahi, y))
		}
	}

	if (over_low == YES) {
		cvblo = CVO(cdat,is)
		cvclo = CV2(cdat,is-1)
	}
	if (over_upp == YES) {
		cvbhi = CVO(cdat,is+1)
		cvchi = CV1(cdat,is+1)
	}

# NEW: assign limits as vectors
	if (overlap == NO) {
		do j = 1, ny {
			y = j
			ia1[j] = cveval (cvalo, y) + 0.5
			ia2[j] = cveval (cvahi, y)
			ic1[j] = cveval (cvclo, y) + 5.5
			ic2[j] = cveval (cvchi, y) - 5.
		}
		call amaxki (ia1,    1, ia1, ny)
		call aminki (ia2, 2048, ia2, ny)
		call amaxki (ic1,    1, ic1, ny)
		call aminki (ic2, 2048, ic2, ny)
		call amovi (ia1, ib1, ny)
		call amovi (ia2, ib2, ny)
	} else {
		do j = 1, ny {
			y = j
			xa = cveval (cvalo, y)
			xb = cveval (cvblo, y)
			xc = cveval (cvclo, y)
			xc = max (xc, xa + 5.)
			xc = max (xc, xb)
			ia1[j] = xa + 0.5
			ib1[j] = xb + 0.5
			ic1[j] = xc + 0.5

			xa = cveval (cvahi, y)
			xb = cveval (cvbhi, y)
			xc = cveval (cvchi, y)
			xc = min (xc, xa - 5.)
			xc = min (xc, xb)
			ia2[j] = xa + 0.5
			ib2[j] = xb		# Note lack of +0.5
			ic2[j] = xc + 0.5	
		}
		call amaxki (ia1,    1, ia1, ny)
		call aminki (ia2, 2048, ia2, ny)
		call amaxki (ib1,    1, ib1, ny)
		call aminki (ib2, 2048, ib2, ny)
		call amaxki (ic1,    1, ic1, ny)
		call aminki (ic2, 2048, ic2, ny)

# The following kills the separate "clean" region
		if (over_upp == YES)
			call amaxi (ic2, ib2, ic2, ny)
		if (over_low == YES)
			call amini (ic1, ib1, ic1, ny)
	}



# Calc size and alloc array
	y = 1.				# XXX hardcode
	x1 = cveval (cvalo, y)
	x2 = cveval (cvahi, y)
	i1 = (x2 - x1) + 0.5
	y = 4096.			# XXX hardcode
	x1 = cveval (cvalo, y)
	x2 = cveval (cvahi, y)
	i2 = (x2 - x1) + 0.5
	nxl = max (i1, i2)
	lenbuf = nxl * ny
	call malloc (bufdat, lenbuf, TY_REAL)
	call malloc (bufwav, lenbuf, TY_REAL)
	call malloc (bufwts, lenbuf, TY_REAL)
	call malloc (buffit, lenbuf, TY_REAL)
	call calloc (buflam, lenbuf, TY_REAL)	 # Different!


# Assign the wavelengths
	ndx = SLITNO(cdat,is)
#	xs1 = X1MM(sdat,ndx)
#	xs2 = X2MM(sdat,ndx)
#	ys1 = Y1MM(sdat,ndx)
#	ys2 = Y2MM(sdat,ndx)

# calculate wavelengths and fill fitting vectors
	ndx = 0				# offset for fitting vectors
#	do j = 1, ny {
	do j = n1, n2 {
	    y = j
	    x1 = cveval (cvalo, y)
	    x2 = cveval (cvahi, y)
	    off = nxl * (j-1)
	    do i = ia1[j], ia2[j] {
		x = i
		f = (x2 - x) / (x2 - x1)	# XXX technically wrong
						# In addition, see (*)
#		xslit =	f * xs1 + (1.-f) * xs2
#		yslit = f * ys1 + (1.-f) * ys2
		call ccd_to_ics (x, y, ccd, chip, xics, yics)
		alp = cveval (CVALP(cdat, is), f)
		gam = cveval (CVALP(cdat, is), f)
		Memr[buflam+off] = qqwave (map, a3, sys, ccd, scaling, alp, gam, xics, yics)
		off = off + 1
	    }



# FOR NOW, do this:
# Move wavelengths into array
	    i1 = ib1[j]
	    off = nxl * (j-1) + (i1 - ia1[j])
	    n = ib2[j]-i1+1
	    call amovr (Memr[buflam+off], wave[i1,j], n)

# Move "clean" regions into fit vectors
	    i1 = ic1[j]
	    n = ic2[j] - i1 + 1
	    call amovr (data[i1,j], Memr[bufdat+ndx], n)
	    call amovr (wave[i1,j], Memr[bufwav+ndx], n)
	    call amovkr (1., Memr[bufwts+ndx], n)
	    ndx = ndx + n

if (mod (j, 800) == 0) {
call eprintf ("QWAVE %d (%d--%d, %d--%d, %d--%d)\n"); call pargi (j)
call pargi (ia1[j])
call pargi (ia2[j])
call pargi (ib1[j])
call pargi (ib2[j])
call pargi (ic1[j])
call pargi (ic2[j])
}

	}

# Check for overwrite problems ...
	off = off - 1
	if (off > lenbuf) {
		call eprintf ("FIT_SLIT_SKY: off (%d) exceeds buf.len (%d)\n")
			call pargi (off)
			call pargi (lenbuf)
		call fatal (0, "BUFFERS exceeded!")
	}
	npt = ndx
	if (npt > lenbuf) {
		call eprintf ("FIT_SLIT_SKY: npt (%d) exceeds buf.len (%d)\n")
			call pargi (npt)
			call pargi (lenbuf)
		call fatal (0, "BUFFERS exceeded!")
	}


# Check for fit input problems ...
	if (npt > 100) {
		call alimr (Memr[bufwav], npt, x1, x2)
call eprintf ("DEBUG: range= %8f %8f\n")
call pargr (x1)
call pargr (x2)
		if (x1 < cvstatr (cvf, CVXMIN) || x2 > cvstatr (cvf, CVXMAX)) {
			call eprintf ("PROBLEM!! wavelegnths (%8f--%8f) exceed range (%8f--%8f)!\n")
				call pargr (x1)
				call pargr (x2)
				call pargr (cvstatr (cvf, CVXMIN))
				call pargr (cvstatr (cvf, CVXMAX))

			call eprintf ("Trying to reconstruct ...\n")
			do i = 0, npt-1 {
			  if (Memr[bufwav+i] == x1) {
			  call eprintf ("Offending wavelen. at ndx %d of %d\n")
				call pargi (i+1)
				call pargi (npt)
			  }
			}
				call fatal (0, "exiting")
		}
	} else {
		call eprintf ("NO GOOD SKY TO FIT -- SKIPPING!\n")
		return
	}


#call eprintf ("DEBUG: fitting %d points\n"); call pargi (npt)
	call cvfit (cvf, Memr[bufwav], Memr[bufdat], Memr[bufwts], npt, WTS_USER, stat)


##################### Iterative section ########################
# XXX set some nominal values; should be flexible!!
	g = 1.3
	rnsq = (4. / g) ** 2
	nsig = 3.
	f = 1. / (nsig*nsig)

# The following seems a necessary step -- get rid of massive CRs so convergence
# can be quicker.

	cutoff = 300.
	call cvvector (cvf, Memr[bufwav], Memr[buffit], npt)
	call asubr (Memr[bufdat], Memr[buffit], Memr[buffit], npt)
	do i = 0, npt-1 {
		if (Memr[buffit+i] > cutoff)	# XXX Just remove uppers
			Memr[bufwts+i] = 0.
	}
	call cvfit (cvf, Memr[bufwav], Memr[bufdat], Memr[bufwts], npt, WTS_USER, stat)
	cutoff =  80.
	call cvvector (cvf, Memr[bufwav], Memr[buffit], npt)
	call asubr (Memr[bufdat], Memr[buffit], Memr[buffit], npt)
	do i = 0, npt-1 {
		if (Memr[buffit+i] > cutoff)	# XXX Just remove uppers
			Memr[bufwts+i] = 0.
	}
	call cvfit (cvf, Memr[bufwav], Memr[bufdat], Memr[bufwts], npt, WTS_USER, stat)


	do j = 1, 2 {			# XXX TMP! Should be dynamic limits
		call cvvector (cvf, Memr[bufwav], Memr[buffit], npt)

		do i = 0, npt-1 {
			res = Memr[bufdat+i] - Memr[buffit+i]
			x = f * res*res / (rnsq + Memr[buffit+i]/g)
			Memr[bufwts+i] = exp (-x*x)
		}
		call cvfit (cvf, Memr[bufwav], Memr[bufdat], Memr[bufwts], npt, WTS_USER, stat)
	}

#call alimr (Memr[buffit], npt, x1, x2)
#call eprintf ("DEBUG: fit  range= %8f %8f\n")
#call pargr (x1)
#call pargr (x2)
#call alimr (Memr[bufdat], npt, x1, x2)
#call eprintf ("DEBUG: data range= %8f %8f\n")
#call pargr (x1)
#call pargr (x2)

	call mfree (bufwts, TY_REAL)
	call mfree (bufwav, TY_REAL)
	call mfree (bufdat, TY_REAL)

# Write sky
#	do j = 1, ny {
	do j = n1, n2 {
		off = nxl * (j-1)
		i1 = ib1[j]
		n = ib2[j] - i1 + 1
		call cvvector (cvf, Memr[buflam+off], Memr[buffit], n)
		call aaddr (sky[i1,j], Memr[buffit], sky[i1,j], n)
	}

	call mfree (buflam, TY_REAL)
	call mfree (buffit, TY_REAL)
end

# FIT_SLIT_INP_ANG:
# The input angles for each slit (as a function of fraction along slit) will
# change slowly and smoothly.  Therefore, rather than calculate them with
# the very high-order fits used in the mappings, we refit each slit with
# a low-order fit in order speed up that part of the optical model.
#
procedure	fit_slit_inp_ang (sdat, cdat, ns, map, a3, sys, ccd, scaling)

pointer	sdat
pointer	cdat
int	ns

pointer	map[8]				# pointers to surf fits (1,2=amap;3,4=b)
double	sys[NPARAM]			# system parameters
real	ccd[NCCD,3]			# CCD geometry
double	a3[3,3]				# grating transform
real	scaling

# From qwave()
double	r[3]
double	alpha, gamma
real	tanx, tany		# should be double; check mapping evals
real	gseval()

real	f[NSLTPTS]			# vector of fractions
real	wt[NSLTPTS]			# vector of weights
real	alp[NSLTPTS]			# vector of alphas
real	gam[NSLTPTS]			# vector of gammas
int	i, j, ndx
int	ier
real	fmin, fmax, finc
real	xs1, xs2, ys1, ys2
real	xmm, ymm

begin
# To begin, we define a vector of fractions. We extend it beyond [0.,1.] for
# stability:

	finc = 1. / (NSLTPTS - 3.)
	fmin = -finc
	fmax = 1. + finc
	do i = 1, NSLTPTS {
		f[i] = fmin + (i-1) * finc
		wt[i] = 1.
	}

	
	do i = 0, ns-1 {
		call cvinit (CVALP(cdat,i), LEGENDRE, NORDSLT, fmin, fmax)
		call cvinit (CVGAM(cdat,i), LEGENDRE, NORDSLT, fmin, fmax)

# Now, fill up the vectors with data
		ndx = SLITNO(cdat,i)
		xs1 = X1MM(sdat,ndx)
		xs2 = X2MM(sdat,ndx)
		ys1 = Y1MM(sdat,ndx)
		ys2 = Y2MM(sdat,ndx)
		do j = 1, NSLTPTS {
			xmm = f[j] * xs1 + (1.-f[j]) * xs2
			ymm = f[j] * ys1 + (1.-f[j]) * ys2
# get mapping and convert to r[3]
			tanx = gseval (map[1], xmm, ymm)
			tany = gseval (map[2], xmm, ymm)

			r[3] = -1.d0 / sqrt (1. + tanx*tanx + tany*tany)
			r[1] = r[3] * tanx
			r[2] = r[3] * tany

# xform into grating system
			call gen_xfm (r, a3, YES)

# convert to alpha,gamma
			alpha = -atan2 (-r[2], -r[3])
			gamma = atan2 (r[1], sqrt (r[3]*r[3]+r[2]*r[2]))

			alp[j] = alpha
			gam[j] = gamma
		}

		call cvfit (CVALP(cdat,i), f, alp, wt, NSLTPTS, WTS_USER, ier)
		if (ier != OK) {
			call eprintf ("FIT_SLIT_INP_ANG: alpha fit err %d\n")
				call pargi (ier)
		}

		call cvfit (CVGAM(cdat,i), f, alp, wt, NSLTPTS, WTS_USER, ier)
		if (ier != OK) {
			call eprintf ("FIT_SLIT_INP_ANG: gamma fit err %d\n")
				call pargi (ier)
		}
	}
end

#
# QQWAVE: get a wavelength based on maps, input angles and pixels coords
#

real	procedure qqwave (map, a3, sys, ccd, scaling, alphin, gammin, xics, yics)


pointer	map[8]			# pointer to maps
double	a3[3,3]			# grating transform
double	sys[NPARAM]		# system parameters
real	ccd[NCCD,3]		# CCD geometry
real	scaling			# scaling adjustment
real	alphin
real	gammin
real	xics, yics		# pixel values in ICS
real	wave			# wavelength (um) returned

double	r[3]
double	alpha, beta, gamma
real	tanx, tany		# should be double; check mapping evals
real	x, y

real	gseval()

begin
	alpha = alphin
	gamma = gammin

# (we can check the gamma as error)

# Convert pixels to beta, gamma:

	x = xics / scaling
	y = yics / scaling

	tanx = gseval (map[7], x, y)
	tany = gseval (map[8], x, y)

	r[3] = -1.d0 / sqrt (1. + tanx*tanx + tany*tany)
	r[1] = r[3] * tanx
	r[2] = r[3] * tany

# xform into grating system
	call gen_xfm (r, a3, YES)

# convert to beta,gamma
	beta = -atan2 (-r[2], -r[3])
	gamma = atan2 (r[1], sqrt (r[3]*r[3]+r[2]*r[2]))

# XXX at this point, gamma is negative wrt gamma from above.  Since cos(gamma)
# is insensitive, it doesn't matter, but it warrants later resolution!!

# Apply the grating equation

	wave = cos (gamma) / (ORDER(sys)*GRLINES(sys)) * (sin (beta) + sin (alpha))
	return (wave)
end

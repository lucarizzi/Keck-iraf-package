# FAKEMAP: Map RA, DEC into X,Y on slitmask -- LRIS version with expanded
# limits for DEIMOS
#
#  7Mar97:  modify to work with the Keck distortion map
# 21Jul97:  revised to add slit_sep and TV guider coords (as in sun version)
# 24sep97:  add slit-ccd mapping; added update to PA for precession
# 26mar98: Allow PCODE=0 objects as reference -- carry through, but no slits
# 06apr98: Passes the updated delta PA.
#
# This version of mapmask calls routines modified for DEIMOS.  Therefore,
# the conversion to slitmask coords does _NOT_ make use of the "touch-up"
# fit, but is done analytically.
#

include	<math.h>
include	<error.h>
include "keck.h"

# LRIS.H -- specific parameters for LRIS
#

define	PPLDIST	19947.D0	# Distance to Telescope Exit Pupil (mm)
define	M_RCURV	2192.1D0	# Mask radius of curvature (mm)
define	M_ANGLE	8.06D0		# Mask angle (tilt) in degrees

define	MM_PIX	0.15767		# mm per pixel
define	ASEC_PIX 0.21		# arcsec per pixel (unbinned)

define	ZPT_YM	197.		# Axial Dist to SMCS y=0 (mm) (305-108)
# MASK_HT0 is calculated as ZPT_YM * tan (M_ANGLE) - MASK_ZINTCPT
define	MASK_HT0 7.0		# Height (mm) above datum at SMCS y=0 (CHECK)
#
define	MASK_X0	88.0  		# Half of Mask X size (215.9) (mm)
define	MASK_Y0	360.0		# Half of Mask Y size (335.6) (mm)
define	Y_BAR	360.0		# y of bar in mm in SMCS
define	HW_BAR	0.1		# half-width of bar in mm
define	FP_XCEN	305.		# Mask x-center (mm) in FP coord
define	FP_YCEN	0.		# Mask y-center (mm) in FP coord
define	MIN_MY	3.5		# Minimum mask y-value
define	MAX_MY	716.5		# Maximum mask y-value

define	CCDXCEN	1024.5		# CCD x-center
define	CCDYCEN	3000.		# CCD y-center
define	PIX_SZ	0.024		# pixel sz (mm)

define		SZ_INLN	128	# Maximum length of input line
define		SZ_ID	32	# Maximum length of ID string
define		CODE_RF	0	# Priority code for reference (addn) objects
define		CODE_GS	-1	# Priority code for guide star
define		CODE_AS	-2	# Priority code for alignment star



define	TV_XOFF		231.8 	# edge of TV CCD (mm, FPCS)
define	TV_YOFF		31.5	# edge of TV CCD (mm, FPCS)
define	TV_MMPIX	0.1835	# pixel scale: mm (FP) / pix (TV)

procedure	t_fakemap()

char	objfile[SZ_FNAME]			# ID, prior, mag, RA, Dec
char	output[SZ_FNAME]			# output file name
real	pa_mask
double	ha0					# Hour angle (field center)
double	pres, wave0, temp
double	epoch					# Epoch for PA precession
real	sep_slit				# minimum separation of slits
real	def_ybot, def_ytop, box_rad		# default (min) sizes
bool	verbose					# verbose output?
pointer	fda, fdb				# in/out file descriptors

int	pcode				# formerly "priority"
double	equinox				# copied value for coord. equinox
real	magn, pa, ybot, ytop
real	delpa				# for precessing PA

# char	sortfile[SZ_FNAME]		# temporary sort file
# pointer	fdx

char	tchar, id[SZ_ID]
int	npts, ndx, i, j
int	nslits				# number of actual slits
real	y1, y2, ytot, ylow, yupp, yb, yt	# used for verbose output
double	guide_ra, guide_dec		# saved coords of guide star
pointer	bufid, bufndx, bufmag
pointer	bufx, bufy, bufy1, bufy2, bufpc, bufpa, bufyb, bufyt
pointer	bufxccd, bufyccd		# ccd vectors

double	cosa, sina
double	pa_tel, sinpat, cospat		# PA of tel center, sin, cos of same
double	ra0, dec0, ra, dec		# field center, object
double	ha				# HA of object
double	rlat				# latitude in radians
double	xaxis				# tot. offset (arcsec) to tel.
double	ra_ref, dec_ref			# refracted coords
double	delra_tel, deldec_tel		# apparent offsets (rad) to tel.
double	ra_tel, dec_tel			# apparent RA, dec (rad) of tel.

double	xfp, yfp			# "ideal" coords in focal plane system
double	x0, y0				# calculated x,y of field center
double	delx, dely			# calculated delta x,y from field center
#double	x, y, ccdx, ccdy
real	x, y, ccdx, ccdy
double	xd, yd				# double for ccd-calc
double	dpa				# double for pa adjust
double	dzero				# double zero for dec kludge
bool	guide_star			# is there a guide star?
bool	specific_pa			# is there a specific PA?
bool	adjust_pa			# adjust the PA?

int	xtv, ytv			# x,y of guide star in TV pixels
real	xgs, ygs			# x,y of guide star in FPCS (mm)

bool	clgetb(), strne()
real	clgetr()
int	map_sort()
int	fscan(), nscan()
pointer	open()

begin
	call clgstr ("objfile", objfile, SZ_FNAME)
	call clgstr ("output", output, SZ_FNAME)
	ha0 = clgetr ("ha0")
	temp = clgetr ("temp")
	pres = clgetr ("pressure")
	wave0 = clgetr ("lambda_cen") * 1.e-4		# in microns
	sep_slit = clgetr ("sep_slit") * MM_ARCS	# sep. bet. slits (mm)
	def_ybot = clgetr ("lower_min")			# min. length in arcsec
	def_ytop = clgetr ("upper_min")			# min. length in arcsec
	box_rad = 0.5 * clgetr ("box_sz")		# box 1/2-length in arcs
	verbose = clgetb ("verbose")
	adjust_pa = clgetb ("adjust_pa")
	epoch = clgetr ("epoch")

        fda = open (objfile, READ_ONLY, TEXT_FILE)

# Count the entries
	ndx = 0
	while (fscan (fda) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#' || nscan() == 0) {
			next
		}
		call reset_scan()
		ndx = ndx + 1
	}
	npts = ndx
	call seek (fda, BOF)

# Allocate arrays
	call malloc (bufx, npts, TY_REAL)
	call malloc (bufy, npts, TY_REAL)
	call malloc (bufy1, npts, TY_REAL)
	call malloc (bufy2, npts, TY_REAL)
	call malloc (bufpc, npts, TY_INT)
	call malloc (bufpa, npts, TY_REAL)
	call malloc (bufyb, npts, TY_REAL)
	call malloc (bufyt, npts, TY_REAL)
	call malloc (bufmag, npts, TY_REAL)
	call malloc (bufxccd, npts, TY_REAL)
	call malloc (bufyccd, npts, TY_REAL)
	call malloc (bufid, npts*SZ_ID, TY_CHAR)
	call amovkc (EOS, Memc[bufid], npts*SZ_ID)
	call malloc (bufndx, npts, TY_INT)

# Now get Field Center and work out its lris coords
	while (fscan(fda) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#' || nscan() == 0) {
			next
		}
		call reset_scan()
		call gargwrd (id, SZ_ID)
		call gargi (pcode)
		call gargd (equinox)
		call gargd (ra0)
		call gargd (dec0)
		call gargr (pa_mask)
		if (nscan() < 6 || strne (id, "CENTER"))
		    call fatal (0, "CENTER not first line or poorly formatted")

# Convert to radians
		ha0  = DEGTORAD (15. * ha0)
		ra0  = DEGTORAD (15. * ra0)
		dec0 = DEGTORAD (dec0)
		break
	}

	rlat = DEGTORAD(OBS_LAT)

# Work out apparent tel. center based on refracted coords of field center
	call rad_refract (ra0, dec0, ha0, rlat, pres, temp, wave0,
								ra_ref, dec_ref)
	cosa = cos (DEGTORAD(-pa_mask))
	sina = sin (DEGTORAD(-pa_mask))

# NB NB: the following is only approx.!  The telescope center needs to be worked
# out more accurately!!
	xaxis = FP_XCEN / MM_ARCS
	delra_tel = DEGTORAD(xaxis/3600.) * cosa / cos (dec_ref)
	deldec_tel = DEGTORAD(xaxis/3600.) * sina

	ra_tel = ra_ref - delra_tel
	dec_tel = dec_ref - deldec_tel
call eprintf ("%12.8f  %12.8f\n")
call pargd (delra_tel*206265.)
call pargd (-deldec_tel*206265.)

# Here's the correct version:
	xaxis = DEGTORAD (FP_XCEN / MM_ARCS / 3600.)	# angular offset to tel.
	dec_tel = asin (cos (xaxis) * sin (dec_ref) -
					sin (xaxis) * cos (dec_ref) * sina)
	delra_tel = asin (sin (xaxis) * cosa / cos (dec_tel))
	ra_tel = ra_ref - delra_tel
call eprintf ("%12.8f  %12.8f\n")
call pargd (delra_tel*206265.)
call pargd ((dec_tel - dec_ref)*206265.)
# ... and corrected PA:
	sinpat = cos (delra_tel) * -sina - sin (delra_tel) * cosa * sin (dec_ref)
	cospat = cos (dec_ref) / cos (dec_tel) * cosa
	pa_tel = RADTODEG (atan2 (sinpat, cospat))
call eprintf ("PA(tel) = %12.8f \n")
call pargd (pa_tel)

#	cosa = cos (DEGTORAD(-pa_tel))
#	sina = sin (DEGTORAD(-pa_tel))

# Don't need cospat, sinpat
# Actually, we want to preserve the PA as pa_tel. This becomes more complicated:
	delra_tel = asin (sin (xaxis) * cosa / cos (dec_ref))
	ra_tel = ra_ref - delra_tel
	dec_tel = atan ((sin (xaxis) * cos (xaxis) * -sina + sin (dec_ref) * cos (dec_ref) * cos (delra_tel)) / (cos (dec_ref) ** 2 - sin (xaxis) **2))

# Here's the calculated field center
	call fp_coord (ra_ref, dec_ref, ra_tel, dec_tel, cosa, sina, xfp, yfp)
call eprintf ("%12.8f %12.8f \n")
call pargd (xfp)
call pargd (yfp)

	call gnom_to_dproj (xfp, yfp, xfp, yfp)		# (allowed: radial corr)
call eprintf ("%12.8f %12.8f \n")
call pargd (xfp)
call pargd (yfp)
# TMP! for tracing (see below)
# call printf ("%9.4f %9.4f \n")
# call pargd (xfp)
# call pargd (yfp)

	dzero = 0.					# kludge for dec double
	call lris_coord (xfp, yfp, dzero, dzero, x0, y0, dzero, dpa)

#
# Loop through the rest of the file, keeping (x,y,pr,pa,y1,y2) and input lines
	ndx = 0
	guide_star = false
	while (fscan(fda) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#' || nscan() == 0) {
			next
		}
		call reset_scan()
		call gargwrd (id, SZ_ID)
		call gargi (pcode)
		call gargr (magn)
		call gargd (ra)
		call gargd (dec)
		call gargr (pa)
		call gargr (ybot)
		call gargr (ytop)
		if (nscan() < 5) {
		    call eprintf ("Poorly-formatted data line -- skipped\n")
		    next
		}

# If guide (pickoff mirror) star, save RA and Dec
		if (pcode == CODE_GS) {
			if (guide_star) {
				call eprintf ("Additional Guide Star ignored\n")
				next
			}
			guide_star = true
			guide_ra  = ra
			guide_dec = dec
		}

# Additional slit info to be carried along: individual PA ...
		if (nscan() < 6 || pa == INDEF) {
			specific_pa = false
			pa = pa_mask
		} else {
			specific_pa = true
		}
		pa = mod ((pa - pa_mask + 360.), 180.)
		if (pa > 90.)
			pa = pa - 180.

# ... and lengths
		if (nscan() < 7 || ybot == INDEF)
			ybot = def_ybot
		if (nscan() < 8 || ytop == INDEF)
			ytop = def_ytop
# Check codes and force the right sizes:
		if (pcode == CODE_GS) {
			ybot = 0.5 / MM_ARCS		# give it 1.0mm
			ytop = 0.5 / MM_ARCS
		} else if (pcode == CODE_AS) {
			ybot = box_rad
			ytop = box_rad
		} else if (pcode == CODE_RF) {
			ybot = 0.			# for safety
			ytop = 0.
		}

# Ready to calculate the mask coords. First convert RA,Dec to radians
		ra  = DEGTORAD (15. * ra)
		dec = DEGTORAD (dec)

# refract the coordinates; work out focal plane, lris coords:
		ha = ha0 + (ra0 - ra)
		call rad_refract (ra, dec, ha, rlat, pres, temp, wave0,
								ra_ref, dec_ref)

		call fp_coord (ra_ref, dec_ref, ra_tel, dec_tel, cosa, sina,
								xfp, yfp)

		call gnom_to_dproj (xfp, yfp, xfp, yfp)		# (allowed)
# TMP!! This is temporary, to set up input for tracing:
# call printf ("%8.3f %8.3f\n")
# call pargd (xfp)
# call pargd (yfp)

		dpa = pa
		call lris_coord (xfp, yfp, x0, y0, delx, dely, dpa, dpa)
		if (specific_pa && adjust_pa)
			pa = dpa

		x = delx + MASK_X0			# not used; TMP?
		y = dely + MASK_Y0			# not used; TMP?

# Check to make sure it's on mask
		if ((y < MIN_MY || y > MAX_MY) && pcode != CODE_GS) {
			call eprintf ("Object %s off mask -- skipped\n")
				call pargstr (id, SZ_ID)
			next
		}

		Memi[bufndx+ndx] = ndx
		Memi[bufpc+ndx] = pcode
		Memr[bufx+ndx] = delx
		Memr[bufy+ndx] = dely
		Memr[bufpa+ndx] = pa
		Memr[bufyb+ndx] = ybot * MM_ARCS
		Memr[bufyt+ndx] = ytop * MM_ARCS
		call strcpy (id, Memc[bufid+ndx*SZ_ID], SZ_ID)
		Memr[bufmag+ndx] = magn
		ndx = ndx + 1
	}
	npts = ndx
	call close (fda)

## OK, we have the (x,y) values. Rewrite input sorted in y..................

	nslits = map_sort (Memr[bufy], npts, Memr[bufx], Memr[bufyb],
			Memr[bufyt], Memi[bufndx], Memi[bufpc], guide_star)

# Guide star (if present) now last in list; if present, save x,y FPCS
	if (Memi[bufpc+npts-1] == CODE_GS) {
		xgs = FP_XCEN - Memr[bufx+npts-1] * cos (DEGTORAD (M_ANGLE))
		ygs = FP_YCEN - Memr[bufy+npts-1]
	}


## Now ready to work out the slit lengths.............................
# Recall, these are based on delta values so we can treat the bar sensibly
# Only objects which fall on mask included in list, and the guide star -- if
# present -- is at end.

# First, assign nominal limits to slits:

	call asubr (Memr[bufy], Memr[bufyt], Memr[bufy2], npts)
	call aaddr (Memr[bufy], Memr[bufyb], Memr[bufy1], npts)

# ... then make the edges butt
	call slit_adjust (Memr[bufy], Memr[bufyb], Memr[bufyt], Memr[bufy1],
				Memr[bufy2], Memi[bufpc], nslits, sep_slit)


# If verbose, print out differences between request and actual lengths:
	if (verbose) {
		call printf ("\nLen('')  Len(mm)   bot.,top    (requested)     del1,del2   Code  Ident\n")

		do i = 0, nslits-1 {
			ndx = Memi[bufndx+i]
			y = Memr[bufy+i]
			y1 = Memr[bufy1+i]
			y2 = Memr[bufy2+i]
			yb = Memr[bufyb+i]
			yt = Memr[bufyt+i]
			ytot = y1 - y2
			ylow = y1 - y  + 4.9e-4		# last term for roundoff
			yupp = y  - y2 + 4.9e-4

			call printf (
	"%6.2f %8.2f %6.2f,%5.2f  (%5.2f,%5.2f)  %5.2f,%5.2f %1s %4d  %s\n")
				call pargr (ytot / MM_ARCS)
				call pargr (ytot)
				call pargr (ylow)
				call pargr (yupp)
				call pargr (yb)
				call pargr (yt)
				call pargr (ylow - yb)
				call pargr (yupp - yt)
# ... work out status code
			    if (yb > ylow || yt > yupp) {
				if ((ylow + yupp) >= (yt + yb))
					call pargc ("^")	# Size OK
				else
					call pargc ("*")	# Too short
			    } else {
				call pargc (" ")		# All OK
			    }
			    call pargi (Memi[bufpc+i])
			    call pargstr (Memc[bufid+ndx*SZ_ID], SZ_ID)
		}
	}

# ... Now update the coords to the center of the mask
	call aaddkr (Memr[bufx],  MASK_X0, Memr[bufx],  npts)
	call aaddkr (Memr[bufy],  MASK_Y0, Memr[bufy],  npts)
	call aaddkr (Memr[bufy1], MASK_Y0, Memr[bufy1], npts)
	call aaddkr (Memr[bufy2], MASK_Y0, Memr[bufy2], npts)

# Find appropriate slit limits for CODE_RF objects:
	do i = nslits, npts-1 {
		if (Memi[bufpc+i] == CODE_RF) {
			y = Memr[bufy+i]
			do j = 0, nslits-1 {
				if (Memr[bufy1+j] > y && Memr[bufy2+j] < y) {
					Memr[bufy1+i] = Memr[bufy1+j]
					Memr[bufy2+i] = Memr[bufy2+j]
					break
				}
			}
		}
	}


# Now work out predicted CCD coordinates
	call ccd_map (Memr[bufx], Memr[bufy], Memr[bufxccd], Memr[bufyccd], npts)

# One final kludge: precession will cause a change in PA -- until precession
# is worked into the whole package, calculate the updated PA for now.
# Precession approximation from Lang, using m,n for 1975.0
	delpa = RADTODEG (atan ((epoch - equinox) * 9.7157e-5 * sin (ra0) / cos (dec0)))

### Finally ready to assemble output file
# Open and write header info
        fdb = open (output, NEW_FILE, TEXT_FILE)

	call fprintf (fdb,
	"###########\n##\n##  NL.MAPMASK:   Input = %s\n##\n###########\n##\n")
		call pargstr (objfile, SZ_FNAME)

	call fprintf (fdb,
	"##     Lambda (A)     T (C)      P (mm Hg)      HA (hr)   \n")
	call fprintf (fdb,
	"##     %7.1f        %4.1f      %7.1f          %4.1f  \n")  
		call pargd (wave0 * 1.e4)
		call pargd (temp)
		call pargd (pres)
		call pargd (RADTODEG(ha0)/15.)

# Make sure that dec is properly formatted for Keck target list
	if (dec0 < 0.) {
		call fprintf (fdb,
			"##\n##  Field_Center.... %011.2h  -%010.1h  %7.2f\n")
	} else {
		call fprintf (fdb,
			"##\n##  Field_Center.... %011.2h   %010.1h  %7.2f\n")
	}
		call pargd (RADTODEG(ra0) / 15.)
		call pargd (abs (RADTODEG(dec0)))
		call pargd (equinox)

	if (guide_star) {
	    if (guide_dec < 0.) {
		call fprintf (fdb,
			"##\n##  Guide_Star...... %011.2h  -%010.1h  %7.2f\n")
	    } else {
		call fprintf (fdb,
			"##\n##  Guide_Star...... %011.2h   %010.1h  %7.2f\n")
	    }
		call pargd (guide_ra)
		call pargd (abs (guide_dec))
		call pargd (equinox)

		xtv = (TV_YOFF - ygs) / TV_MMPIX + 0.5
		ytv = (xgs - TV_XOFF) / TV_MMPIX + 0.5
		call fprintf (fdb,
		    "##\n##  Guide Star, Slit-Viewing TV: %dx %dy \n")
			call pargi (xtv)
			call pargi (ytv)
	} else {
		call fprintf (fdb, "##\n##  Guide_Star...... (none)\n")
	}

	call fprintf (fdb,
	"##\n##     Mask Position Angle: %8.3f deg in %7.2f   (%6.1f in %6.1f)\n##\n")
		call pargr (pa_mask + delpa)
		call pargd (epoch)
		call pargr (pa_mask)
		call pargd (equinox)

	call fprintf (fdb, "##\n## Xobj   Yobj     Ymin    Ymax   Rel_PA  CCDx   CCDy  Code  Mag  ID\n##\n")

# Loop through objects/slit and print:
	do i = 0, npts-1 {
		ndx = Memi[bufndx+i]
		call fprintf (fdb,
		    "%7.3f %7.3f %8.3f %7.3f %6.2f %7.1f %6.1f %4d %5.2f %s\n")
			call pargr (Memr[bufx+i])
			call pargr (Memr[bufy+i])
			call pargr (Memr[bufy1+i])
			call pargr (Memr[bufy2+i])
			call pargr (Memr[bufpa+ndx])
			call pargr (Memr[bufxccd+i])
			call pargr (Memr[bufyccd+i])
			call pargi (Memi[bufpc+i])
			call pargr (Memr[bufmag+ndx])
			call pargstr (Memc[bufid+ndx*SZ_ID], SZ_ID)
	}

	call close (fdb)

	call mfree (bufndx, TY_INT)
	call mfree (bufid, TY_CHAR)
	call mfree (bufpc, TY_INT)
	call mfree (bufyccd, TY_REAL)
	call mfree (bufxccd, TY_REAL)
	call mfree (bufmag, TY_REAL)
	call mfree (bufyt, TY_REAL)
	call mfree (bufyb, TY_REAL)
	call mfree (bufpa, TY_REAL)
	call mfree (bufpc, TY_REAL)
	call mfree (bufy2, TY_REAL)
	call mfree (bufy1, TY_REAL)
	call mfree (bufy, TY_REAL)
	call mfree (bufx, TY_REAL)

end



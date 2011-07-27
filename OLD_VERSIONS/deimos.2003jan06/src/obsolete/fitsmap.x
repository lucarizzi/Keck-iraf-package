# FITSMAP: This is a version of MAPMASK that outputs FITS tables
# MAPMASK: Map RA, DEC into X,Y on slitmask
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
include	<time.h>
include "keck.h"
include "deimos.h"
include	"ftab.h"

define	TV_XOFF		231.8 	# edge of TV CCD (mm, FPCS)
define	TV_YOFF		31.5	# edge of TV CCD (mm, FPCS)
define	TV_MMPIX	0.1835	# pixel scale: mm (FP) / pix (TV)

procedure	t_fitsmap()

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

#####
char	fmtstr[SZ_LINE]				# format for writing table row
int	fmlen
char	date[SZ_TIME]					# time/date of file
pointer	kdat				# FITS table data structure
# int	inul
real	rnul
# double	dnul
char	snul[1]
#####
double	ra_pnt, dec_pnt
pointer	bufa, bufd, bufp		# FITS vector pointers
char	creatask[60]
char	person[30]
char	unknown[3]

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
#	call strcpy ("-99999", inul, 6)
#	call strcpy ("-9.e99", rnul, 6)
#	call strcpy ("-9.e99", dnul, 6)

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

####### FOR FITS TABLES
	call malloc (bufa, npts, TY_DOUBLE)
	call malloc (bufd, npts, TY_DOUBLE)
	call malloc (bufp, npts, TY_REAL)
####### END

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
####
	ra_pnt = ra0
	dec_pnt = dec0
####
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

######## SAVE COORD's for FITS TABLE
		Memd[bufa+ndx] = ra
		Memd[bufd+ndx] = dec
		if (specific_pa)
			Memr[bufp+ndx] = pa
		else
			Memr[bufp+ndx] = INDEF
######## END

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

######## SET UP FITS FILE
        fdb = open (output, NEW_FILE, TEXT_FILE)

	fmlen = SZ_LINE
	call strcpy ("", snul, 1)
	call strcpy ("Phillips <phillips@ucolick.org>", person, 35)
	call strcpy ("hacked up version of mapmask", creatask, 35)
	call strcpy ("???", unknown, 3)

	call ft_date (date, SZ_TIME)

##### TMP!! Write the primary HDU
	call ftab_whdu0 (fdb)
########

########  WRITE ObjectCat TABLE
	call ftab_init (kdat, 18, fmtstr)

	call ftcol_def (kdat, "ObjectId",  "I6", snul, snul,  6, fmtstr, fmlen)
	call ftcol_def (kdat, "OBJECT",   "A68", snul, snul, 68, fmtstr, fmlen)
	call ftcol_def (kdat, "RA_OBJ", "F12.8","deg", rnul, 12, fmtstr, fmlen)
	call ftcol_def (kdat, "DEC_OBJ","F12.8","deg", rnul, 12, fmtstr, fmlen)
	call ftcol_def (kdat, "RADESYS",   "A8", snul, snul, 8, fmtstr, fmlen)
	call ftcol_def (kdat, "EQUINOX", "F8.3",   "a", rnul,  8, fmtstr, fmlen)
#WORK HERE...
	call ftcol_def (kdat, "MJD-OBS", "F11.3",   "d", 11, fmtstr, fmlen)
	call ftcol_def (kdat, "mag",      "F7.3", snul,  7, fmtstr, fmlen)
	call ftcol_def (kdat, "pBand",      "A6", snul,  6, fmtstr, fmlen)
	call ftcol_def (kdat, "RedShift", "F10.7",snul, 10, fmtstr, fmlen)
	call ftcol_def (kdat, "MajAxis", "F9.2", "arcsec", 9, fmtstr, fmlen)
	call ftcol_def (kdat, "MajAxPA", "F8.2",    "deg", 8, fmtstr, fmlen)
	call ftcol_def (kdat, "MinAxis", "F9.2", "arcsec", 9, fmtstr, fmlen)
	call ftcol_def (kdat, "PM_RA",  "F9.4", "arcsec/a", 9, fmtstr, fmlen)
	call ftcol_def (kdat, "PM_DEC", "F9.4", "arcsec/a", 9, fmtstr, fmlen)
	call ftcol_def (kdat, "Parallax", "F7.4", "arcsec", 7, fmtstr, fmlen)
	call ftcol_def (kdat, "ObjClass",  "A20", snul, 20, fmtstr, fmlen)
	call ftcol_def (kdat, "CatFilePK",  "I6", snul,  6, fmtstr, fmlen)


	call ftab_whead (fdb, kdat, npts)

	call pkwstr (fdb, kdat, "EXTNAME", "ObjectCat", snul)
	call pkwi   (fdb, kdat, "EXTVER", 0, snul)
	call pkwi   (fdb, kdat, "NPRIKEY", 1, snul)
	call pkwstr (fdb, kdat, "PKTYP1", "ObjectId", snul)

	call pkwi   (fdb, kdat, "NFORHDU", 1, snul)
	call pkwstr (fdb, kdat, "HDLOC1", snul, snul)
	call pkwstr (fdb, kdat, "HDXTN1", "TABLE", snul)
	call pkwstr (fdb, kdat, "HDNAM1", "CatFiles", snul)

	call pkwi   (fdb, kdat, "NFORKEY", 1, snul)
	call pkwi   (fdb, kdat, "FKHDU1", 1, snul)
	call pkwstr (fdb, kdat, "FKTYP1", "CatFilePK", snul)
	call pkwstr (fdb, kdat, "FKFYP1", "CatFilePK", snul)

	call pkwstr (fdb, kdat, "DATE", date, snul)
	call pkwstr (fdb, kdat, "AUTHOR", person, snul)
	call pkwstr (fdb, kdat, "CREATOR", creatask, snul)
	call pkwstr (fdb, kdat, "CATNAME", objfile, snul)

#call eprintf ("%d: %s\n")
#call pargi (TWUSED(kdat))
#call pargstr (fmtstr)
#call eprintf ("*%d %d*\n")
#call pargi (ndx)
#call pargi (TROWPT(kdat))

# Write out the OBJECT CAT TABLE
	do ndx = 0, npts-1 {
		call sprintf (BUFROW(kdat), TWUSED(kdat), fmtstr)
			call pargi (ndx+1)
			call pargstr (Memc[bufid+ndx*SZ_ID])
			call pargd   (Memd[bufa+ndx] * 15.)
			call pargd   (Memd[bufd+ndx])
			call pargstr (unknown)
			call pargd   (equinox)
			call pargr   (0.)	# (INDEF)
			call pargr   (Memr[bufmag+ndx])
			call pargstr (unknown)
			call pargr   (0.)	# (INDEF)
			call pargr   (0.)	# (INDEF)
			call pargr   (Memr[bufp+ndx])
			call pargr   (0.)	# (INDEF)
			call pargr   (0.)	# (INDEF)
			call pargr   (0.)	# (INDEF)
			call pargr   (0.)	# (INDEF)
			if (Memi[bufpc+ndx] == CODE_GS)
				call pargstr ("Guide_Star")
			else if (Memi[bufpc+ndx] == CODE_AS)
				call pargstr ("Alignment_Star")
			else
				call pargstr ("Program_Target")
			call pargi   (1)
			
		call ftab_wrow (fdb, kdat)
	}
	call ftab_free (fdb, kdat)

########

##### WRITE CatFiles TABLE
	call ftab_init (kdat, 2, fmtstr)

	call ftcol_def (kdat, "CatFilePK",     "I6", snul,   6, fmtstr, fmlen)
	call ftcol_def (kdat, "CatFileName", "A255", snul, 255, fmtstr, fmlen)

	call ftab_whead (fdb, kdat, 1)

	call pkwstr (fdb, kdat, "EXTNAME", "CatFiles", snul)
	call pkwi   (fdb, kdat, "EXTVER", 0, snul)
	call pkwi   (fdb, kdat, "NPRIKEY", 1, snul)
	call pkwstr (fdb, kdat, "PKTYP1", "CatFilePK", snul)

	call pkwstr (fdb, kdat, "DATE", date, snul)
	call pkwstr (fdb, kdat, "AUTHOR", person, snul)
	call pkwstr (fdb, kdat, "CREATOR", creatask, snul)
	call pkwstr (fdb, kdat, "CATNAME", objfile, snul)

#call eprintf ("%d: %s\n")
#call pargi (TWUSED(kdat))
#call pargstr (fmtstr)
#call eprintf ("*%d %d*\n")
#call pargi (ndx)
#call pargi (TROWPT(kdat))

# Write out the OBJECT CAT TABLE
	call sprintf (BUFROW(kdat), TWUSED(kdat), fmtstr)
		call pargi (1)
		call pargstr ("Some name of a catalog file????")
			
	call ftab_wrow (fdb, kdat)

	call ftab_free (fdb, kdat)


########  WRITE MaskDesign TABLE
	call ftab_init (kdat, 12, fmtstr)

	call ftcol_def (kdat, "DesId",     "I11", snul, 11, fmtstr, fmlen)
	call ftcol_def (kdat, "DesName",   "A68", snul, 68, fmtstr, fmlen)
	call ftcol_def (kdat, "DesAuth",   "A68", snul, 68, fmtstr, fmlen)
	call ftcol_def (kdat, "DesCreat",  "A68", snul, 68, fmtstr, fmlen)
	call ftcol_def (kdat, "DesDate",   "A19", snul, 19, fmtstr, fmlen)
	call ftcol_def (kdat, "RA_PNT",  "F12.8", "deg", 12, fmtstr, fmlen)
	call ftcol_def (kdat, "DEC_PNT", "F12.8", "deg", 12, fmtstr, fmlen)
	call ftcol_def (kdat, "RADEPNT",    "A8", snul,  8, fmtstr, fmlen)
	call ftcol_def (kdat, "EQUINPNT","F13.6",   "a", 13, fmtstr, fmlen)
	call ftcol_def (kdat, "PA_PNT",  "F12.7", "deg", 12, fmtstr, fmlen)
	call ftcol_def (kdat, "DATE_PNT",  "A19", snul, 19, fmtstr, fmlen)
	call ftcol_def (kdat, "LST_PNT",   "F8.3", "deg",  8, fmtstr, fmlen)

	call ftab_whead (fdb, kdat, 1)

	call pkwstr (fdb, kdat, "EXTNAME", "MaskDesign", snul)
	call pkwi   (fdb, kdat, "EXTVER", 0, snul)
	call pkwi   (fdb, kdat, "NPRIKEY", 1, snul)
	call pkwstr (fdb, kdat, "PKTYP1", "DesId", snul)

	call pkwstr (fdb, kdat, "DATE", date, snul)
	call pkwstr (fdb, kdat, "AUTHOR", person, snul)
	call pkwstr (fdb, kdat, "CREATOR", creatask, snul)

	call sprintf (BUFROW(kdat), TWUSED(kdat), fmtstr)
		call pargi (1)
		call pargstr (output)
		call pargstr (person)
		call pargstr ("DesCreat?")
		call pargstr (date)
		call pargd (ra_pnt * 15.)
		call pargd (dec_pnt)
		call pargstr (unknown)
		call pargd (equinox)
		call pargr (pa_mask)
		call pargstr (date)
		call pargd (RADTODEG(ha0))
			
	call ftab_wrow (fdb, kdat)

	call ftab_free (fdb, kdat)


########  WRITE MaskBlu TABLE
	call ftab_init (kdat, 16, fmtstr)

	call ftcol_def (kdat, "BluId",     "I11", snul, 11, fmtstr, fmlen)
	call ftcol_def (kdat, "DesId",     "I11", snul, 11, fmtstr, fmlen)
	call ftcol_def (kdat, "BluName",   "A68", snul, 68, fmtstr, fmlen)
	call ftcol_def (kdat, "BluObsvr",  "A68", snul, 68, fmtstr, fmlen)
	call ftcol_def (kdat, "BluCreat",  "A68", snul, 68, fmtstr, fmlen)
	call ftcol_def (kdat, "BluDate",   "A19", snul, 19, fmtstr, fmlen)
	call ftcol_def (kdat, "LST_Use",  "F8.3","deg",  8, fmtstr, fmlen)
	call ftcol_def (kdat, "Date_Use",  "A19", snul, 19, fmtstr, fmlen)
	call ftcol_def (kdat, "TELESCOP",  "A68", snul, 68, fmtstr, fmlen)
	call ftcol_def (kdat, "RefrAlg",   "A68", snul, 68, fmtstr, fmlen)
	call ftcol_def (kdat, "AtmTempC", "F5.1","degC", 5, fmtstr, fmlen)
	call ftcol_def (kdat, "AtmPres",  "F6.1","mmHg", 6, fmtstr, fmlen)
	call ftcol_def (kdat, "AtmHumid", "F5.3", "K/m", 5, fmtstr, fmlen)
	call ftcol_def (kdat, "AtmTTLap", "F7.5", snul,  7, fmtstr, fmlen)
	call ftcol_def (kdat, "RefWave",  "F7.2", "nm",  7, fmtstr, fmlen)
	call ftcol_def (kdat, "DistMeth",  "A68", snul, 68, fmtstr, fmlen)

	call ftab_whead (fdb, kdat, 1)

	call pkwstr (fdb, kdat, "EXTNAME", "MaskBlu", snul)
	call pkwi   (fdb, kdat, "EXTVER", 0, snul)
	call pkwi   (fdb, kdat, "NPRIKEY", 1, snul)
	call pkwstr (fdb, kdat, "PKTYP1", "BluId", snul)

	call pkwi   (fdb, kdat, "NFORHDU", 1, snul)
	call pkwstr (fdb, kdat, "HDLOC1", snul, snul)
	call pkwstr (fdb, kdat, "HDXTN1", "TABLE", snul)
	call pkwstr (fdb, kdat, "HDNAM1", "MaskDesign", snul)

	call pkwi   (fdb, kdat, "NFORKEY", 1, snul)
	call pkwi   (fdb, kdat, "FKHDU1", 1, snul)
	call pkwstr (fdb, kdat, "FKTYP1", "DesId", snul)
	call pkwstr (fdb, kdat, "FKFYP1", "DesId", snul)

	call pkwstr (fdb, kdat, "DATE", date, snul)
	call pkwstr (fdb, kdat, "OBSERVER", person, snul)
	call pkwstr (fdb, kdat, "CREATOR", creatask, snul)

	call sprintf (BUFROW(kdat), TWUSED(kdat), fmtstr)
		call pargi (1)
		call pargi (1)
		call pargstr (output)
		call pargstr (person)
		call pargstr (creatask)
		call pargstr (date)
		call pargd   (RADTODEG(ha0))
		call pargd   (date)		# TMP: should be for use
		call pargstr ("Keck II")
		call pargstr ("RefrAlg_TBD")
		call pargr   (temp)
		call pargr   (pres)
		call pargr   (0.)		# TMP: no clue
		call pargr   (0.)		# TMP: no clue
		call pargr   (0.1*wave0)
		call pargr   ("DistAlg_TBD")
			
	call ftab_wrow (fdb, kdat)

	call ftab_free (fdb, kdat)

##############*END*#####################################

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

########
########  WRITE DesiSlits TABLE
	call wftab_4 (fdb, kdat, date, person, creatask, pa_mask,
		Memr[bufy], Memr[bufy1], Memr[bufy2], Memr[bufpa],
		Memi[bufpc], Memi[bufndx], nslits, Memd[bufa], Memd[bufd], npts)

	call mfree (bufp, TY_REAL)
	call mfree (bufd, TY_DOUBLE)
	call mfree (bufa, TY_DOUBLE)
########


# Find appropriate slit limits for CODE_RF objects:
## FITS TABLE: At this point, ndx is no longer needed; reassign for SlitObjMap
#	do i = nslits, npts-1 {
	do i = 0, npts-1 {
		if (Memi[bufpc+i] == CODE_RF) {
			y = Memr[bufy+i]
			do j = 0, nslits-1 {
				if (Memr[bufy1+j] > y && Memr[bufy2+j] < y) {
					Memr[bufy1+i] = Memr[bufy1+j]
					Memr[bufy2+i] = Memr[bufy2+j]
					Memi[bufndx+i] = j + 1
					break
				}
			}
		} else if (Memi[bufpc+i] == CODE_GS) {
			Memi[bufndx+i] = 0
		} else {
			Memi[bufndx+i] = i + 1
		}
	}


# Now work out predicted CCD coordinates
#	call ccd_map (Memr[bufx], Memr[bufy], Memr[bufxccd], Memr[bufyccd], npts)

# One final kludge: precession will cause a change in PA -- until precession
# is worked into the whole package, calculate the updated PA for now.
# Precession approximation from Lang, using m,n for 1975.0
	delpa = RADTODEG (atan ((epoch - equinox) * 9.7157e-5 * sin (ra0) / cos (dec0)))

########
########  WRITE SlitObjMap TABLE
	call wftab_5 (fdb, kdat, date, person, creatask, Memi[bufndx], npts)

	call close (fdb)
########

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

## DESISLITS table:
# This is only APPROX!  Take size, find approx delra, deldec to center of slit;
# apply to ra,dec.  Also, PA is PA on mask surface, not sky.

procedure	wftab_4 (fd, kdat, date, person, creatask, maskpa, y, y1, y2, slitpa, pcode, index, nslits, ra, dec, npts)

pointer	fd
pointer	kdat
char	date[ARB]
char	person[ARB]
char	creatask[ARB]
real	maskpa
real	y[nslits]
real	y1[nslits], y2[nslits]
real	slitpa[nslits]
int	pcode[nslits]
int	index[nslits]
int	nslits
double	ra[npts], dec[npts]
int	npts

char	fmtstr[SZ_LINE]
int	fmlen

char	acode
real	dely, ylen
real	pa
real	raoff, decoff
double	racen, deccen

int	i

begin

########  WRITE DesiSlits TABLE
	fmlen = SZ_LINE

	call ftab_init (kdat, 9, fmtstr)

	call ftcol_def (kdat, "dSlitId",   "I11",    "", 11, fmtstr, fmlen)
	call ftcol_def (kdat, "DesId",     "I11",    "", 11, fmtstr, fmlen)
	call ftcol_def (kdat, "slitRA",  "F12.8", "deg", 12, fmtstr, fmlen)
	call ftcol_def (kdat, "slitDec", "F12.8", "deg", 12, fmtstr, fmlen)
	call ftcol_def (kdat, "slitTyp",    "A1",    "",  1, fmtstr, fmlen)
	call ftcol_def (kdat, "slitLen", "F11.3", "arcsec", 11, fmtstr, fmlen)
	call ftcol_def (kdat, "slitLPA",  "F8.3", "deg",  8, fmtstr, fmlen)
	call ftcol_def (kdat, "slitWid", "F11.3", "arcsec", 11, fmtstr, fmlen)
	call ftcol_def (kdat, "slitWPA",  "F8.3", "deg",  8, fmtstr, fmlen)

	call ftab_whead (fd, kdat, nslits)

	call pkwstr (fd, kdat, "EXTNAME", "DesiSlits", "")
	call pkwi   (fd, kdat, "EXTVER", 0, "")
	call pkwi   (fd, kdat, "NPRIKEY", 1, "")

	call pkwstr (fd, kdat, "PKTYP1", "DesId", "")
	call pkwstr (fd, kdat, "PKFYP1", "DesId", "")

	call pkwi   (fd, kdat, "NFORKEY", 1, "")
	call pkwi   (fd, kdat, "NFORHDU", 1, "")

	call pkwstr (fd, kdat, "FKTYP1", "DesId", "")
	call pkwstr (fd, kdat, "FKFYP1", "DesId", "")
	call pkwi   (fd, kdat, "FKHDU1", 1, "")
	call pkwstr (fd, kdat, "HDLOC1", "", "")
	call pkwstr (fd, kdat, "HDXTN1", "TABLE", "")
	call pkwstr (fd, kdat, "HDNAM1", "MaskDesign", "")

	call pkwstr (fd, kdat, "DATE", date, "")
	call pkwstr (fd, kdat, "AUTHOR", person, "")
	call pkwstr (fd, kdat, "CREATOR", creatask, "")

	do i = 1, nslits {
		ylen = y1[i] - y2[i]
		dely = 0.5 * ylen - (y[i] - y2[i])
		ylen = ylen / cos (DEGTORAD(slitpa[i])) / 0.7253
		dely = dely / cos (DEGTORAD(slitpa[i])) / 0.7253
		
		pa = slitpa[i] + maskpa
# very approx!
		raoff  = ylen * sin (DEGTORAD(pa)) / cos (DEGTORAD(dec[i]))
		decoff = ylen * cos (DEGTORAD(pa))

		racen  = ra[index[i]+1] * 15. + raoff / 3600.	#  0-indexed
		deccen = dec[index[i]+1] + decoff / 3600.	#  0-indexed

		if (pcode[i] == CODE_GS || pcode[i] == CODE_RF)
			acode = 'G'
		else if (pcode[i] == CODE_AS)
			acode = 'A'
		else
			acode = 'P'
		
		call sprintf (BUFROW(kdat), TWUSED(kdat), fmtstr)
			call pargi (i)
			call pargi (1)
			call pargd (racen)
			call pargd (deccen)
			call pargc (acode)
			call pargr (ylen)
			call pargr (pa)
			call pargr (0.95)
			call pargr (maskpa+90.)
			
		call ftab_wrow (fd, kdat)
	}

	call ftab_free (fd, kdat)
	
end


procedure	wftab_5 (fd, kdat, date, person, creatask, index, npts)

pointer	fd
pointer	kdat
char	date[ARB]
char	person[ARB]
char	creatask[ARB]
int	index[npts]
int	npts

char	fmtstr[SZ_LINE]
int	fmlen

int	i

begin

########  WRITE SlitObj TABLE
	fmlen = SZ_LINE

	call ftab_init (kdat, 3, fmtstr)

	call ftcol_def (kdat, "DesId",     "I11", "", 11, fmtstr, fmlen)
	call ftcol_def (kdat, "ObjectId",  "I11", "", 11, fmtstr, fmlen)
	call ftcol_def (kdat, "dSlitId",   "I11", "", 11, fmtstr, fmlen)

	call ftab_whead (fd, kdat, npts)

	call pkwstr (fd, kdat, "EXTNAME", "SlitObjMap", "")
	call pkwi   (fd, kdat,  "EXTVER", 0, "")
	call pkwi   (fd, kdat, "NPRIKEY", 2, "")

	call pkwstr (fd, kdat, "PKTYP1", "DesId", "")
	call pkwstr (fd, kdat, "PKTYP2", "ObjectId", "")

	call pkwi   (fd, kdat, "NFORKEY", 3, "")
	call pkwi   (fd, kdat, "NFORHDU", 3, "")

	call pkwstr (fd, kdat, "FKTYP1", "DesId", "")
	call pkwstr (fd, kdat, "FKFYP1", "DesId", "")
	call pkwi   (fd, kdat, "FKHDU1", 1, "")
	call pkwstr (fd, kdat, "HDLOC1", "", "")
	call pkwstr (fd, kdat, "HDXTN1", "TABLE", "")
	call pkwstr (fd, kdat, "HDNAM1", "MaskDesign", "")

	call pkwstr (fd, kdat, "FKTYP2", "dSlitId", "")
	call pkwstr (fd, kdat, "FKFYP2", "dSlitId", "")
	call pkwi   (fd, kdat, "FKHDU2",  2, "")
	call pkwstr (fd, kdat, "HDLOC2", "", "")
	call pkwstr (fd, kdat, "HDXTN2", "TABLE", "")
	call pkwstr (fd, kdat, "HDNAM2", "DesiSlits", "")

	call pkwstr (fd, kdat, "FKTYP3", "ObjectId", "")
	call pkwstr (fd, kdat, "FKFYP3", "ObjectId", "")
	call pkwi   (fd, kdat, "FKHDU3", 3, "")
	call pkwstr (fd, kdat, "HDLOC3", "", "")
	call pkwstr (fd, kdat, "HDXTN3", "TABLE", "")
	call pkwstr (fd, kdat, "HDNAM3", "ObjectCat", "")

	call pkwstr (fd, kdat, "DATE", date, "")
	call pkwstr (fd, kdat, "AUTHOR", person, "")
	call pkwstr (fd, kdat, "CREATOR", creatask, "")

	do i = 1, npts {
	    call sprintf (BUFROW(kdat), TWUSED(kdat), fmtstr)
		call pargi (1)
		call pargi (i)
		call pargi (index[i])
			
	    call ftab_wrow (fd, kdat)
	}

	call ftab_free (fd, kdat)
	
end




#
# MAP_SORT: Sort the map info, separating the file by SLIT objects, REF objects
# and the Guide Star.  This task could be a lot cleaner if we stored the
# table as an array rather than individual vectors.
#

int	procedure map_sort (y, npt, x, yb, yt, ndx, pcode, guide_star)

real	y[npt]					# y vector (to sort by)
int	npt					# length
real	x[npt]					# x vector (to sort by)
real	yb[npt], yt[npt]			# y-limit info
int	ndx[npt]				# index to ID
int	pcode[npt]				# priority code
bool	guide_star				# Is there a guide star?

int	i, i2
int	code
int	ndx2
int	nref					# number of reference objects
int	nobj					# number of program targets

begin
# We start by separating the reference objects and guide star from the others:

# if guide star present, adjust so it has space at end
	if (guide_star) {
		i2 = npt - 1
	} else {
		i2 = npt
	}


# now go through list backwards, switching reference,guide objects towards end
	nref = 0
	nobj = 0
	ndx2 = i2
	do i = i2, 1, -1 {
		code = pcode[i]
		if (code == CODE_GS) {		# If GS, switch to end; test new
		    	call rndxswitch (y, ndx, npt, i, npt)
			code = pcode[npt]
		}
		if (code == CODE_RF) {
# if in last position, or next also ref, don't switch
			if (i < ndx2)
				call rndxswitch (y, ndx, npt, i, ndx2)
			nref = nref + 1
			ndx2 = ndx2 - 1
		} else {
			nobj = nobj + 1
		}
	}
	ndx2 = ndx2 + 1				# first ref object


# Now sort the slits
	call ndxsortr (y, ndx, nobj)

# ... and reference objects
	if (nref > 1)
		call ndxsortr (y[ndx2], ndx[ndx2], nref)

	call re_orderr (ndx, x, x, npt, YES)
	call re_orderr (ndx, yb, yb, npt, YES)
	call re_orderr (ndx, yt, yt, npt, YES)
	call re_orderi (ndx, pcode, pcode, npt, YES)

# and return value of slits
	return (nobj)
end


# SLIT_ADJUST: lengthen/shorten slits to use maximum space
# Recall, these are based on delta values so we can treat the bar sensibly
# Only objects which fall on mask included in list, and guide star is at end
# At this point, all values in mm

procedure	slit_adjust (y, ybot, ytop, y1, y2, pcode, n, sep_slit)

real	y[n]				# y-object
real	ybot[n], ytop[n]		# min. y length, bottom and top
real	y1[n], y2[n]			# actual ymin, ymax
int	pcode[n]			# priority/type code
int	n				# length of object list
real	sep_slit			# spacing between slits

int	i
int	imax
real	ymax, del
real	wt1, wt2, sum			# weights to distribute del

begin
	ymax = y[1]
	imax = 1
	do i = 2, n {
		if (y[i] > ymax) {
			ymax = y[i]
			imax = i
		}
	}

	if (pcode[1] != CODE_AS)
		y2[1] = MIN_MY - MASK_Y0	# put at edge; ref'd to MASK_Y0
	if (pcode[imax] != CODE_AS)
		y1[imax] = MAX_MY - MASK_Y0

	do i = 2, imax {
		del = y2[i] - y1[i-1] - sep_slit

# check for alignment stars (guide star should come at end)
		if (pcode[i] == CODE_AS)
			wt2 = 0.
		else if (del < 0.)
			wt2 = 0.5
		else
			wt2 = ytop[i]

		if (pcode[i-1] == CODE_AS)
			wt1 = 0.
		else if (del < 0.)
			wt1 = 0.5
		else
			wt1 = ybot[i-1]

		sum = wt1 + wt2
		if (sum == 0.)				# both AS, no changes
			next

		y2[i]   = y2[i]   - del * wt2 / sum	# Note sign!
		y1[i-1] = y1[i-1] + del * wt1 / sum
	}

# Now adjust at the bar so that no slit crosses the mask center:
	do i = 2, imax {
		if (y[i] > 0. && y[i-1] < 0.) {
		    if (pcode[i] != CODE_AS && pcode[i] != CODE_GS)
			y2[i] = 0.
		    if (pcode[i-1] != CODE_AS && pcode[i-1] != CODE_GS)
			y1[i-1] = 0.
		}
	}
end
# Ideally, we should keep a buffer of "reserves" for each end, so that
# we can adjust for a "min_length" condition, too.


#
# RE_ORDERR: reorder a real vector by index (input and output may be same)
#
procedure re_orderr (ndx, in, out, n, zero_index)

int	ndx[n]				# index of positions
real	in[n]				# input vector to be reordered
real	out[n]				# output vector -- reodered by index
int	n				# length of vectors
bool	zero_index			# Is the vector zero-indexed?

int	i
int	ioff				# offset for zero indexing
pointer	sp, buf				# pointer for stack, work vector

begin
# allocate stack memory for work vector:
	call smark (sp)
	call salloc (buf, n, TY_REAL)

# do the reordering:
	if (zero_index)
		ioff = 1
	else
		ioff = 0

	do i = 1, n
		Memr[buf+i-1] = in[ndx[i]+ioff]

# copy to output vector
	call amovr (Memr[buf], out, n)

	call sfree (sp)
end

#
# RE_ORDERI: reorder an integer vector by index (input and output may be same)
#
procedure re_orderi (ndx, in, out, n, zero_index)

int	ndx[n]				# index of positions
int	in[n]				# input vector to be reordered
int	out[n]				# output vector -- reodered by index
int	n				# length of vectors
bool	zero_index			# Is the vector zero-indexed?

int	i
int	ioff				# offset for zero indexing
pointer	sp, buf				# pointer for stack, work vector

begin
# allocate stack memory for work vector:
	call smark (sp)
	call salloc (buf, n, TY_INT)

# do the reordering:
	if (zero_index)
		ioff = 1
	else
		ioff = 0

	do i = 1, n
		Memi[buf+i-1] = in[ndx[i]+ioff]

# copy to output vector
	call amovi (Memi[buf], out, n)

	call sfree (sp)
end

#
# NDXSORTR: Sort a list of real values (low to high) and carry index
# (see also re_orderX)

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
#
# RNDXSWITCH: Switch 2 lines of real items and index 
#

procedure	rndxswitch (r, ndx, n, i1, i2) 

real	r[n]				# Real vector
int	ndx[n]				# Index vector
int	n				# number of points
int	i1, i2				# indices of items to be switched

real	h
int	ih

begin
	h = r[i1]
	r[i1] = r[i2]
	r[i2] = h

	ih = ndx[i1]
	ndx[i1] = ndx[i2]
	ndx[i2] = ih
end


# DEIMOS.H -- specific parameters for DEIMOS  XXX NEEDS TO BE REVIEWED!!!

###
### TELESCOPE PARAMS THAT VARY BY INSTRUMENT:
###

# define FL_TEL	150393.95D0	# Foc Len of KII WITH 3" PULLBACK (T8_im.zmx)
# define FL_TEL	150100.4D0	# Foc Len of KII WITH 3" PULLBACK (REF)
# define FL_TEL	150310.5D0	# FUDGED FOC LEN to give 0.0014 more
define	FL_TEL	150327.0D0	# FUDGED FOC LEN to give 0.00011 more than above

define	R_IMSURF 2133.6D0	# radius of image surface (mm) 
					# NB - should match ray trace program!!
#define	 DIST_C0 4.535218e-4 	# distortion term in r (Sutin 3" Ray Tr)
define	DIST_C0	0.0e-4 		# distortion term in r scale should be 0
define	DIST_C2	-1.111311e-8	# distortion term in r**3  (Sutin 3" RT)

define	MM_ARCS	0.7277		# mm/arcsec in focal plane  (gnomonic)
					## XXX DOESN'T CONFORM!! NEEDS ATTN!
					## May be moot: obsolete??
define	PPLDIST	20018.4D0	# Distance to Telescope Exit Pupil (mm) (REF)

###
### SLITMASK PARAMETERS:
###

define	M_ANGLE	6.00D0		# Mask angle (tilt) in degrees
define	M_RCURV	2120.9D0	# Mask radius of curvature (mm)	[83.5in]
define	ZPT_YM	128.803		# Dist to tel.axis, in SMCS-XXX (mm) 5.071in

# MASK_HT0 is calculated as ZPT_YM * tan (M_ANGLE) - MASK_ZINTCPT (0.4+0.026in)
define	MASK_HT0 2.717	# Height (mm) above datum at SMCS y=0
# MASK_HT0 is calculated as ZPT_YM * tan (M_ANGLE) - MASK_ZINTCPT (-0.4in)
# define	MASK_HT0 3.378		# Height (mm) above datum at SMCS y=0


###
### DETECTOR PARAMETERS:
###

define	PIX_SZ		0.015		# pixel sz (mm)
define	CCDXPIX		 2048.		# 2048. pixels wide
define	CCDYPIX		 4096.		# 4096. pixels deep
define	CCDXEDG		0.154		# 84 micron to P3 + 70 um to sawline
define	CCDYEDG		0.070		# 40 micron to P3 + 30 um to sawline
define	NOMXGAP		   1.		# 1 mm gap
define	NOMYGAP		 0.1		# 0.1 mm gap
define	FCSYEDG		0.140		# 140 um active to sawline
define	FCSYPIX		  600.		# 600 pixels high as used


###
### MISC:
###

define		SZ_INLN	128	# Maximum length of input line
define		SZ_ID	32	# Maximum length of ID string
define		CODE_RF	0	# Priority code for reference (addn) objects
define		CODE_GS	-1	# Priority code for guide star
define		CODE_AS	-2	# Priority code for alignment star


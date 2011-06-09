# LRIS.H -- specific parameters for LRIS  XXX NEEDS TO BE REVIEWED!!!

###
### TELESCOPE PARAMS THAT VARY BY INSTRUMENT:
###

# define FL_TEL	150393.95D0	# Foc Len of KII WITH 3" PULLBACK (T8_im.zmx)
# define TEL_FOC	149583.D0		# focal length (mm)
define	FL_TEL	149600.4D0	# Nominal LRIS Tel foc -- UPDATE!!
define	R_IMSURF 2133.6D0	# radius of image surface (mm) 
					# NB - should match ray trace program!!
#define DIST_C0	4.535218e-4 	# distortion term in r (Sutin 3" Ray Tr)
define	DIST_C0	0.0e-4 		# distortion term in r scale should be 0
define	DIST_C2	-1.111311e-8	# distortion term in r**3  (Sutin 3" RT)

define	MM_ARCS	0.7253		# mm/arcsec in focal plane  (gnomonic)
					# shortcut from FL_TEL/206264.8
define	PPLDIST	19968.9D0	# Distance to Telescope Exit Pupil (mm) (REF)
## ZEMAX says 19948.0 w/o ADC
## ZEMAX says 19941.2 w/  ADC
## ZEMAX says R_IMSURF=2149.2 w/o ADC
## ZEMAX says R_IMSURF=2146.5 w 3" pullback
## ZEMAX says PPLDIST=20023.15 w 3" pullback
## ZEMAX says FL_TEL=150155.7 w 3" pullback


###
### SLITMASK PARAMETERS:
###

define	M_ANGLE	8.06D0		# Mask angle (tilt) in degrees
define	M_RCURV	2192.1D0	# Mask radius of curvature (mm)

## There is a confusion about whether ZPT_YM is truly in SMCS (tilted) XXX
# Let us define the zeropoint as 7.0 arcmin
define	ZPT_YM	307.309		# Dist to tel.axis, in SMCS-XXX (mm)

# MASK_HT0 is calculated as ZPT_YM * sin?(M_ANGLE) - MASK_ZINTCPT (20.9mm)
define	MASK_HT0 22.188		# Height (mm) above datum at SMCS y=0

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


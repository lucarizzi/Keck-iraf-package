# DASTROM: Derive astrometric positions from DEIMOS images.  Based on QTRMAPings

# STRUCTURE:
# 1. Read in list of RA,DEC,xpx ypx for reference stars
# 1a. refract and produce XYARCS
#   need to adjust tel axis and rotation to match
# 2. For list of pixel locations, produce XARCS, YARCS
# 3. Calculate sky coords using XYARCS, telax, PA  -- DONE
# 3a unrefract

# TBD:
# 1. Input DEIMOS pointing coords
# 2. Refract/unrefract
# 3. Entry of reference stars; xform, mod's to telax and PA
# 4. TEST IT

include	<math.h>
include	<math/gsurfit.h>
include	"deimos.h"
include	"instrument.h"
include	"dsimulator.h"

procedure	t_dastrom()

char	input[SZ_FNAME]				# input file of x,y pairs, ID
char	inref[SZ_FNAME]				# input file of ref stars
char	amap[SZ_FNAME]			# input alpha--sky map
char	bmap[SZ_FNAME]			# input ICS--beta map
pointer	fd
pointer	fda, fdb

double	sys[NPARAM]				# system parameters
real	ccd[NCCD,3]				# CCD geometry
double	a3[3,3]				# grating transform

real	xpix, ypix
int	n
char	id_obj[SZ_ID]			# object name
char	remain[SZ_LINE]			# remainder of line

real	xics, yics
real	xas, yas
pointer	asfx, asfy			# pointers to surface fits (amap)
pointer	bsfx, bsfy			# pointers to surface fits (bmap)

double	r			# radius of object from tel-axis
double	phi			# PA on sky from tel-axis
double	sind, sina		# sin of dec, delta-RA
double	ra0, dec0, pa0		# RA, Dec and PA on axis
double	ra_obj, dec_obj
double	x, y
pointer	indat

bool	strne()
int	fscan(), nscan()
real	clgetr()
pointer	open()

begin
	call clgstr ("input", input, SZ_FNAME)
	fd = open (input, READ_ONLY, TEXT_FILE)

	MU(sys) = DEGTORAD (clgetr ("mu"))		# XXX, ect
	GR_YERR(sys) = DEGTORAD (clgetr ("roll3"))
	GR_ZERR(sys) = 0.

	call clgstr ("amap", amap, SZ_FNAME)
	fda = open (amap, READ_ONLY, TEXT_FILE)

	call clgstr ("bmap", bmap, SZ_FNAME)
	fdb = open (bmap, READ_ONLY, TEXT_FILE)

# Initialize the maps
	call gs_ingest (fda, asfx, asfy)
	call gs_ingest (fdb, bsfx, bsfy)

# set up the grating transform
	call ccd_geom (ccd, sys)
	call gsetup (a3, sys)

# Allocate the indat data structure and TELDAT vector (lifted from dsim)
	call malloc (indat, NINDAT, TY_STRUCT)
	call malloc (PTTELDAT(indat), NTELPAR, TY_DOUBLE)

# Get the field center info:
	RA_FLD(indat) = DEGTORAD (clgetr ("ra0") / 15.)
	DEC_FLD(indat) = DEGTORAD (clgetr ("dec0"))
	PA_ROT(indat) = DEGTORAD(clgetr ("pa_fld"))

# need to worry about refraction here??  Probably
	call fld2telax (indat)

# These need defining:
	ra0  = RA_TEL(indat)
	dec0 = DEC_TEL(indat)
	pa0  = PA_ROT(indat)

	while (fscan (fd) != EOF) {
		call gargr (xpix)
		call gargr (ypix)
		call gargwrd (id_obj, SZ_ID)
		call gargstr (remain, SZ_LINE)

		if (nscan() < 3) {
			call eprintf ("Poorly formatted input line\n")
			next
		}

# get mapping into ICS pixels
		call pane_to_ics (xpix, ypix, ccd, n, xics, yics)

# get mapping into x,y arcsec
		call qrevmod0 (xics, yics, a3, sys, asfx, asfy, bsfx, bsfy, xas, yas)

# the following lifted from dsim.sky_coords ...
		x = xas
		y = yas
		r = sqrt (x*x + y*y)
		r = atan (r/206264.8D0)

		phi = pa0 - atan2 (y, x)	# WORK

		sind = sin (dec0) * cos (r) + cos (dec0) * sin (r) * cos (phi)

		sina = sin (r) * sin (phi) / sqrt (1. - sind*sind)

		dec_obj = asin (sind)
		ra_obj = ra0 + asin (sina)
		if (ra_obj < 0.)
			ra_obj = ra_obj + TWOPI
		else if (ra_obj >= TWOPI)
			ra_obj = ra_obj - TWOPI

# Print it all out
		call printf ("%-16s %13.3h %13.2h %-s ## %6.1f %6.1f\n")
			call pargstr (id_obj)
			call pargd (RADTODEG(ra_obj)/15.d0)
			call pargd (RADTODEG(dec_obj))
			call pargstr (remain)
			call pargr (xpix)
			call pargr (ypix)
	}

	call close (fdb)
	call close (fda)
	call close (fd)
end

#
# QREVMOD0: Do the reverse model, 0th order, to sky
#

procedure	qrevmod0 (xics, yics, a3, sys, asfx, asfy, bsfx, bsfy, xas, yas)

real	xics, yics			# pixel values in ICS
double	sys[NPARAM]			# system parameters
double	a3[3,3]				# grating transform
pointer	asfx, asfy, bsfx, bsfy		# pointers to surface fits
real	xas, yas			# X,Y in slitmask

double	r[3]
double	alpha, beta, gamma
real	tanx, tany			# should be double; check mapping evals

real	gseval()
begin
# Get mapping and convert to r[3]
	tanx = gseval (bsfx, xics, yics)
	tany = gseval (bsfy, xics, yics)

	r[3] = -1.d0 / sqrt (1. + tanx*tanx + tany*tany)
	r[1] = r[3] * tanx
	r[2] = r[3] * tany

# xform into grating system
	call gen_xfm (r, a3, YES)

# convert to beta,gamma (signs may not be right)
	beta = -atan2 (-r[2], -r[3])
	gamma = atan2 (r[1], sqrt (r[3]*r[3]+r[2]*r[2]))

# Apply the grating equation
	alpha = -beta

# convert alpha, gamma into x,y,z (cf Schroeder p259); note sign reversal of alpha
	r[1] = sin (gamma)
	r[2] = sin (-alpha) * cos (gamma)
	r[3] = cos (-alpha) * cos (gamma)

# xform out of grating system
	call gen_xfm (r, a3, NO)

# convert to tanx, tany
	tanx = (-r[1] / -r[3])
	tany = (-r[2] / -r[3])
# call eprintf ("tanx,tany: %5f %5f\n")
# call pargr (tanx)
# call pargr (tany)

# get mapping into Slitmask coord
	xas = gseval (asfx, tanx, tany)
	yas = gseval (asfy, tanx, tany)
end


procedure	pane_to_ics (xpix, ypix, ccd, n, xics, yics)

real	xpix, ypix			# Input PANE values
real	ccd[NCCD,3]			# CCD geometrical params
int	n				# returned Chip number
real	xics, yics			# returned ICS values

int	nx, ny				# zero-indexed loc. in array
real	xccd, yccd			# CCD values w/in pane, normal orient.
double	cosa, sina, bx, by
double	x, y
double	xp, yp

int	det_chip()
begin
# determine chip number
	n = det_chip (xpix, ypix)

# (the following does not take pix centers into account, but hopefully no one
#  is doing astrometry within 0.5 px of chip edge)
	nx = (xpix-1.) / CCDXPIX
	ny = (ypix-1.) / CCDYPIX
	xccd = xpix - nx * CCDXPIX
	yccd = ypix - ny * CCDYPIX
	
# get the CCD info
	cosa = cos (ccd[n,3])
	sina = sin (ccd[n,3])
	bx = ccd[n,1]			# Refer to chip center
	by = ccd[n,2]			# Refer to chip center

# calculate ICS
	x = xccd - 0.5 * (CCDXPIX + 1.)
	y = yccd - 0.5 * (CCDYPIX + 1.)

	xp =  cosa * x - sina * y + bx
	yp =  sina * x + cosa * y + by

	xics = xp
	yics = yp

end

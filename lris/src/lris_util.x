# LRIS utility routines, including:
#
# -- ATM_GEOM: Works out Horizon-based geometry and Atm. dispersion
# -- DEWAR_ELEV: Works out Horizon-based geometry of dewar; returns elevation
# -- CHK_ELEV: Checks the elevation of the telescope against limits


include <math.h>
include "lris.h"

#
# ATM_GEOM: Works out Horizon-based geometry and Atm. dispersion
# NB: all REAL inputs
#

procedure atm_geom (ha, dec, lambda1, lambda2, lat, zret, azret, pa, atmdisp)

real	ha, dec, lat
real	lambda1, lambda2
real	zret, azret, pa, atmdisp

double	sinlat, coslat
double	sinha, cosha, sindec, cosdec
double	sina, cosa, sinz, cosz
double	z, dec1, dec2, ha1, ha2, dha, ddec
real	q, delndx, r1, r2

begin
# Work out the atmospheric refraction constants for Keck (p=610, t=0)
	q = (1.e4/lambda1) ** 2
	delndx = (54.47 + 24976.8/(146.-q) + 216.3/(41.-q)) * 1.e-6
	r1 = delndx / (1. + 2.*delndx)		# approx.

	q = (1.e4/lambda2) ** 2
	delndx = (54.47 + 24976.8/(146.-q) + 216.3/(41.-q)) * 1.e-6
	r2 = delndx / (1. + 2.*delndx)		# approx.

	r1 = RADTODEG(r1) * 3600.	# Convert to arcsec
	r2 = RADTODEG(r2) * 3600.	# Convert to arcsec

# Work out the geometry:
	ha1 = DEGTORAD (ha*15.)
	dec1 = DEGTORAD (dec)

	sinlat = sin (DEGTORAD (lat))
	coslat = cos (DEGTORAD (lat))
	sinha = sin (ha1)
	cosha = cos (ha1)
	sindec = sin (dec1)
	cosdec = cos (dec1)

	cosz = sindec * sinlat + cosdec * coslat * cosha
	sinz = sqrt (1 - cosz * cosz)
	sina = (-cosdec * sinha) / sinz
	cosa = (sindec*coslat - cosdec*cosha*sinlat) / sinz

	z = acos (cosz)
	zret = RADTODEG(z)
	azret = RADTODEG(atan2 (sina, cosa))

# Estimate dispersion (first order is fine...)
	atmdisp = tan (z) * abs (r2 - r1)

# Now, subtract a small amount (1") from z:
	z = z - DEGTORAD (1./3600.)
	sinz = sin (z)
	cosz = cos (z)
 	sindec = cosz*sinlat + sinz*cosa*coslat
	cosdec = sqrt (1. - sindec*sindec)

	sinha = (-sinz * sina) / cosdec
	cosha = (cosz*coslat - sinz*cosa*sinlat) / cosdec

	ha2 = atan2 (sinha, cosha)
	dec2 = asin (sindec)
	
	dha = mod ((ha1 - ha2 + PI), double (TWOPI)) - PI
	dha = dha * cosdec
	ddec = dec1- dec2
	pa = atan2 (ddec, dha)
	pa = RADTODEG (pa) + 90.

end

#
# DEWAR_ELEV: Works out Horizon-based geometry of dewar; returns elevation
# NB: all REAL inputs
# NB: simple geometry -- no atmospheric corrections
#

real	procedure dewar_elev (ha, dec, pa, lat, axis)

real	ha, dec
real	pa			# PA of dewar (mask PA + 90)
real	lat, axis		# Latitude, angle of collimator-to-camera
real	elret			# returned dewar elevation

real	ha1, ha2
real	sinlat, coslat
real	sinpa, cospa, sindec, cosdec, cosax, sinax
real	sindec2, cosdec2, sindha, cosdha
real	dec1
real	cosz, z
# real	azret, sina, cosa, sinz, dec2 

begin

# Work out the geometry:
	ha1 = DEGTORAD (ha*15.)
	dec1 = DEGTORAD (dec)

	sinlat = sin (DEGTORAD (lat))
	coslat = cos (DEGTORAD (lat))
	sinpa = sin (DEGTORAD (pa))
	cospa = cos (DEGTORAD (pa))
	sinax = sin (DEGTORAD (axis))
	cosax = cos (DEGTORAD (axis))
	sindec = sin (dec1)
	cosdec = cos (dec1)

	sindec2 = cosax * sindec - sinax * cosdec * cospa
	cosdec2 = sqrt (1. - sindec2*sindec2)
	sindha = sinax * sinpa / cosdec2
	cosdha = (cosax - sindec * sindec2) / (cosdec * cosdec2)
	ha2 = ha1 + atan2 (sindha, cosdha)
	cosz = sindec2 * sinlat + cosdec2 * coslat * cos (ha2)

	z = acos (cosz)
	elret = 90. - RADTODEG(z)
	return (elret)

end

#
# CHK_ELEV: Checks the elevation of the telescope against limits
# Simple model to check elevation (only checks first and last)
#

real	procedure chk_elev (az1, z1, az2, z2)

real	az1, z1, az2, z2		# Azimuth, Z-dist at start and end

begin
	if (max (z1, z2) > SRV_ZMX) {
		if (SRV_AZ1 < max (az1, az2) && min (az1, az2) < SRV_AZ2) {
			return (1.)
		} else if (max (z1, z2) > GEN_ZMX) {
			return (2.)
		} else {
			return (0.)
		}
	}
end

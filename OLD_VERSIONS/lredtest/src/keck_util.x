# UNDER REVISION (first untested; remaining unreviewed)
# Note that GNOM_TO_DPROJ is used in dsim.x
#
# KECK utility routines, including:
#
# -- CHK_ELEV: Checks the elevation of the telescope against limits
# -- ATM_GEOM: Works out Horizon-based geometry and Atm. dispersion
# -- ATM_REFR: Returns atm refraction in arcsec for given lambda
# -- RAD_REFRACT:  calculate refracted celestial coordinates


include <math.h>
include "keck.h"

#
# CHK_ELEV: Checks the elevation of the telescope against limits.
# Real precision adequate.
#

real	procedure chk_elev (ha1, ha2, dec, ha_lim, z1, z2)

real	ha1, ha2			# Hour angle at start, end of exposure
real	dec				# declination
real	ha_lim				# HA at limit
real	z1, z2				# Z-dist at start, end

real	sinlat, coslat, sind, cosd	# sin, cos of lat, dec
real	sinh, cosh			# sin, cos HA
real	sinz, cosz, sinaz, cosaz	# sin, cos of zenith dist, azimuth
real	az				# azimuth at limit

begin
	sinlat = sin (DEGTORAD (OBS_LAT))
	coslat = cos (DEGTORAD (OBS_LAT))
	sind = sin (DEGTORAD (dec))
	cosd = cos (DEGTORAD (dec))
	z1 = RADTODEG (acos (sind*sinlat + cosd*coslat*cos (DEGTORAD(ha1*15.))))
	z2 = RADTODEG (acos (sind*sinlat + cosd*coslat*cos (DEGTORAD(ha2*15.))))

# first check for service platform; limiting HA (in hours) is:
# (if we ever return to K1, note that ha_lim is negative of this)
	sinz = sin (DEGTORAD (SRV_ZMX))
	cosz = cos (DEGTORAD (SRV_ZMX))
	cosh = (cosz - sind * sinlat) / (cosd*coslat)
	sinh = sqrt (1. - cosh * cosh)
	ha_lim = RADTODEG (acos (cosh)) / 15.

	sinaz = -cosd * sinh / sinz
	cosaz = (sind * coslat - cosd * sinlat * cosh) / sinz
	az = atan2 (sinaz, cosaz)

# is azimuth between platform limits?  If so, return...

	if (SRV_AZ2 > az && az > SRV_AZ1 && ha2 > ha_lim)
		return (LM_SRVZ)

# ... if not, then check for general elevation limit ...
	cosz = cos (DEGTORAD (GEN_ZMX))
	cosh = (cosz - sind * sinlat) / (cosd*coslat)
	ha_lim = RADTODEG (acos (cosh)) / 15.

	if (ha1 > 12.) {
		if (ha1 < ha_lim + 24.)
			return (LM_GENZ)
	} else {
		if (ha2 > ha_lim)
			return (LM_GENZ)
	}
		
# ... apparently OK!
	return (LM_OK)
end

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
#real	q, delndx, r1, r2
real	atm_refr()

begin

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
	atmdisp = tan (z) * abs (atm_refr (lambda2) - atm_refr (lambda1))
	

# Now, subtract a small amount (6") from z (refraction is always to decr. z):
	z = z - DEGTORAD (6./3600.)
	sinz = sin (z)
	cosz = cos (z)
 	sindec = cosz*sinlat + sinz*cosa*coslat
	cosdec = sqrt (1. - sindec*sindec)

	sinha = (-sinz * sina) / cosdec
	cosha = (cosz*coslat - sinz*cosa*sinlat) / cosdec

	ha2 = atan2 (sinha, cosha)
	dec2 = asin (sindec)

	dha = mod ((ha2 - ha1 + PI), double (TWOPI)) - PI
	dha = dha * cosdec
	ddec = dec2 - dec1
	pa = RADTODEG (atan2 (-dha, ddec))

end

#
# ATM_REFR: returns refraction (arcsec) at a particular wavelength;
# user must apply factor = TANZ * cos(ypa - adpa) / (arcsec/pix * Y_binning)
#

real	procedure atm_refr (lambda)

real	lambda				# central wavelength for slit (A)

real	adjust
double	q, delndx, r
real	p, t

begin
# From Allen, p.124. Note lambda0 should technically be _vacuum_

# adjustment from standard conditions:
	p = PRESSURE
	t = TYPTEMP
	adjust = p * (1. + p * (1.049e-6 - t * 0.0157e-6)) /
					(720.883 * (1. + t * 0.003661))

	q = (1.e4/lambda) ** 2
	delndx = (64.328 + 29498.1/(146.-q) + 255.4/(41.-q)) * 1.e-6
	delndx = adjust * delndx

	r = delndx / (1. + 2.*delndx)		# approx.
	r = r * 206265.				# Convert to arcsec

	return (r)
end

procedure	rad_refract (ra, dec, h, lat, pres, temp, wave, ara, adec)
#
#     Calculate the apparent declination and right ascention corrected
#     for atmospheric refraction, given the true RA and DEC, etc.
#     RA and DEC are in radians.
#
#     Parameters -    (">" input, "<" output)
#
#     (<) ara        (Real*8) Apparent Right Ascention in radians
#
#     (<) adec       (Real*8) Apparent Declination in radians
#  
#     (>) ra         (Real*8) True Right Ascention in radians
#
#     (>) dec        (Real*8) True Declination in radians
#
#     (>) h          (Real*8) Hour angle in radians
#  
#     (>) lat        (Real*8) Observer's latitude in radians
#
#     (>) pres       (Real*8) Local air pressure in mm Hg.
#
#     (>) temp       (Real*8) Local air temperature in Celsius.
#   
#     (>) wave       (Real*8) Wavelength in micrometers
#
#     (>) vp         (Real*8) Vapor pressure in mm Hg (not currently used.)
#                                         [ACP: because formula NOT in Allen?]
#     (<) zd         (Real*8) Zenith distance in degrees
#
#     (<) n          (Real*8) Atmospheric refractivity
#                                         [ACP: ie `R0'; approx. index(air)-1]
#     (<) r          (Real*8) Refraction coefficient
#
#     History:  CRO/CIT.  20 June 1988.  Original.
#  Converted to spp by A. Phillips, Nov 95
#
double	ara, adec, ra, dec, lat, pres, temp, wave
double	h, da, dd, num, tcorr, tanz
double	wavers, cosq, sinq, sina, sinz, cosz
double	zd, n, r
#

begin

# Altitude and Zenith Distance
	sina = sin (lat) * sin (dec) + cos (lat) * cos (dec) * cos (h)
	cosz = sina
	sinz = sqrt (1. - sina * sina)
	tanz = sinz / cosz

	if  (sinz != 0.) {
		sinq = cos (lat) * sin (h) / sinz
		cosq =  (cos (HALFPI - lat) - cosz * cos (HALFPI - dec)) / 
						(sinz * sin (HALFPI - dec))
	} else {
#  zenith distance is 0 so q doesn't matter. r=0 anyway.
		sinq = 0.
		cosq = 0.
	}

# refractive index at standard TP
	wavers = 1. / (wave * wave)
	n = 64.328 + (29498.1 / (146. - wavers)) +  (255.4 / (41. - wavers)) 
	n = 1.e-6 * n

# temperature and pressure correction
	tcorr = 1. + 0.003661 * temp
	num   = 720.88 * tcorr
	n = n * ((pres * (1. + (1.049 - 0.0157 * temp) * 1.e-6 * pres)) / num)
	
	zd = asin (sinz)
	zd = RADTODEG (zd)

# r, refraction
	r = n * 206265. * tanz

# changes in ra and dec in radians.
	da = DEGTORAD (r / 3600.) * sinq / cos (dec)
	dd = DEGTORAD (r / 3600.) * cosq

# corrected ra and dec
	ara =  ra  + da
	adec = dec + dd
end


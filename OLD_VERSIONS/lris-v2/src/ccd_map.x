include <math/gsurfit.h>

# CCD_MAP: map x,y (slitmask) into CCD coords. Assumes fit from geomap entered
# in def_map()

# To produce the coeffs, run geomap on a file that has x,y,X,Y (CCD,SM).
# The two columns of coeffs are for the x-fit and y-fit respectively.
# Note that the linear terms (surface1) must be added to the first 3
# coeffs of the distortion terms (surface2)

procedure	ccd_map (xmm, ymm, xccd, yccd, npts)

real	xmm[npts], ymm[npts]			# input slitmask coords
real	xccd[npts], yccd[npts]			# output ccd coords

pointer	sfx, sfy

begin
	call def_ccd_map (sfx, sfy)
	call gsvector (sfx, xmm, ymm, xccd, npts)
	call gsvector (sfy, xmm, ymm, yccd, npts)
	call gsfree (sfx)
	call gsfree (sfy)

end

# DEF_CCD_MAP: define the slit-ccd mapping

procedure	def_ccd_map (sfx, sfy)

# These coefficients were taken from a fit to the M71_D1 map.

pointer	sfx, sfy			# pointers to surface fits in x,y

int	ncoeff
pointer	xcoeff, ycoeff			# coeff's in x,y

begin
	ncoeff = 20
	call malloc (xcoeff, ncoeff, TY_REAL)
	call malloc (ycoeff, ncoeff, TY_REAL)

	Memr[xcoeff  ]  = 2.
	Memr[xcoeff+1]  = 4.
	Memr[xcoeff+2]  = 3.
	Memr[xcoeff+3]  = 1.
	Memr[xcoeff+4]  = 1.
	Memr[xcoeff+5]  = 2048.
	Memr[xcoeff+6]  = 1.
	Memr[xcoeff+7]  = 2048.
	Memr[xcoeff+8]  = -0.03770211951349326 + 98.4300681858445
	Memr[xcoeff+9]  = -0.4415477110265099  + 158.3504279124463
	Memr[xcoeff+10] = -0.03673128448546001 + 0.4223812500419837
	Memr[xcoeff+11] = -0.3239684789820529
	Memr[xcoeff+12] = -0.0124154686668477
	Memr[xcoeff+13] = -0.08896883339120662
	Memr[xcoeff+14] = -0.02163392758515836
	Memr[xcoeff+15] = 0.007576498655666576
	Memr[xcoeff+16] = -0.1595099437642113
	Memr[xcoeff+17] = -0.6104127399575367
	Memr[xcoeff+18] = -0.03350676226625679
	Memr[xcoeff+19] = 0.03939341278295458

	Memr[ycoeff  ]  = 2.
	Memr[ycoeff+1]  = 3.
	Memr[ycoeff+2]  = 4.
	Memr[ycoeff+3]  = 1.
	Memr[ycoeff+4]  = 1.
	Memr[ycoeff+5]  = 2048.
	Memr[ycoeff+6]  = 1.
	Memr[ycoeff+7]  = 2048.
	Memr[ycoeff+8]  = 0.02331238975350751  +  179.4542329239503
	Memr[ycoeff+9]  = -0.01657133908694772 +    0.3298482036463803
	Memr[ycoeff+10] = 0.01524631501641042  + -156.8215724169838
	Memr[ycoeff+11] = 0.2782043893954322
	Memr[ycoeff+12] = 0.5754596900630351
	Memr[ycoeff+13] = 0.4655777759949629
	Memr[ycoeff+14] = 0.1552991074183973
	Memr[ycoeff+15] = -0.05722616917061861
	Memr[ycoeff+16] = 0.1082981782594116
	Memr[ycoeff+17] = 0.2773043587328404
	Memr[ycoeff+18] = 0.0833984099970036
	Memr[ycoeff+19] = -0.07470371597586598

	call gsrestore (sfx, Memr[xcoeff], ncoeff)
	call gsrestore (sfy, Memr[xcoeff], ncoeff)

	call mfree (xcoeff, TY_REAL)
	call mfree (ycoeff, TY_REAL)
end

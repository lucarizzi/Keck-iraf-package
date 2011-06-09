include	<imhdr.h>
include	<error.h>
include <math/curfit.h>
 
# T_LPROCESS: Do a 1-pass bias/trim/flat-field, for double/single-amp mode -- 
# including optional flip in single amp. THIS IS EXTENDED l2process
#
#  Needs a check for bwidth, cskip inside bounds!
# Should allow optional windowing --- means redefining x1,x2,y1,y2
#
# 13jul98 -- added bias-structure frame (if now very cumbersome, not v. logical)
# 26apr00 -- (ver 0.1) added the scale factors to header; fixed to preserve
#             intscale effect (that is, if fact1 or fact2 > 0.5, cut down)

define		SCALE	0.5		# Scaling for use with SHORT format
define		SZ_FITS	60		# size of FITS char string
define		FIT_BK	2		# Backoff for bias fit

procedure t_lprocess()

char	image[SZ_FNAME]
char	flatim[SZ_FNAME]
char	biasim[SZ_FNAME]
char	imageout[SZ_FNAME], imorig[SZ_FNAME]
char	history[SZ_LINE]
bool	flatten					# flatten the image?
bool	intscale				# half and store as short?
int	x1, y1, x2, y2				# region of interest
int	l1, l2					# lines of bias strip
int	osa1, osa2				# 1st overscan/bias columns
int	osb1, osb2				# 2nd overscan/bias columns
int	nxlow					# last x of lower half of CCD
pointer	im1, im2, im3, im4
long	v1[IM_MAXDIM], v2[IM_MAXDIM], v3[IM_MAXDIM], v4[IM_MAXDIM]
pointer	line1, line2, line3, line4

bool	bias_sub				# subtract bias-structure image
int	nx, ny					# size of output image
real	biasa, biasb				# bias levels 1 and 2
int	xoff, i
int	bwid, bmid, bmed
int	nxf, nyf, nxb, nyb
real	bsum
real	scaling
pointer	bufb

# New variables
bool	reverse					# Flip x?
char	kwchar[SZ_FITS]
int	prepix, postpix				# default pre/postpix (per amp)
int	binx, biny				# binning in x,y
int	ival, ip, len
int	namps					# number of amps
int	nord					# number of spline segments
int	ndx, ier
real	ratio, fact1, fact2			# scaling factors for amp-gain
pointer	blev1, blev2, ypix			# pointers to bias level(s), y
pointer	wts					# pointers to weights (unused)
pointer	cv					# pointer to curvefit

bool	clgetb(), streq()
int	clgeti(), impnlr(), imgnlr(), asoki(), imgeti()
int	nscan(), strlen(), imaccf()
real	clgetr(), imgetr()
pointer	immap(), imgs2i()

begin
	call clgstr ("image", image, SZ_FNAME)
	call clgstr ("out_image", imageout, SZ_FNAME)
	flatten = clgetb ("flatten")
	intscale = clgetb ("intscale")
	reverse = clgetb ("reverse")
	nord = clgeti ("norder")
#	ratio = clgetr ("ratio2_1")

#	x1 = clgeti ("x1")
#	y1 = clgeti ("y1")
#	x2 = clgeti ("x2")
#	y2 = clgeti ("y2")
#	nxlow = clgeti ("nxlow")
#	low1 = clgeti ("low1")
#	low2 = clgeti ("low2")
#	upp1 = clgeti ("upp1")
#	upp2 = clgeti ("upp2")

#	nx = x2 - x1 + 1
#	ny = y2 - y1 + 1


# Open FLATFIELD IMAGE if requested:
	if (flatten) {
		call clgstr ("flatim", flatim, SZ_FNAME)
		if (streq (flatim, "")) {
			call eprintf ("WARNING: not flattened!\n")
			flatten = false
		} else {
			im2 = immap (flatim, READ_ONLY, 0)
			nxf = IM_LEN(im2,1)
			nyf = IM_LEN(im2,2)
		}
	}

# Open BIAS-STRUCTURE IMAGE if requested:
	call clgstr ("biasim", biasim, SZ_FNAME)
	if (streq (biasim, "")) {
		call eprintf ("WARNING: No bias structure!\n")
		bias_sub = false
	} else {
		bias_sub = true
		im4 = immap (biasim, READ_ONLY, 0)
		nxb = IM_LEN(im4,1)
		nyb = IM_LEN(im4,2)
	}

	im1 = immap (image, READ_ONLY, 0)

# test to see if processed already
	if (imaccf (im1, "LPROCESS") == YES)
		call fatal (0, "IMAGE ALREADY PROCESSED!")

# Get relevant keywords ...
# ... for the number of prescan columns     # ?? FIX - need string
	if (imaccf (im1, "NUMAMPS") == YES) {
		namps = imgeti (im1, "NUMAMPS")
	} else if (imaccf (im1, "AMPLIST") == YES) {
		call imgstr (im1, "AMPLIST", kwchar, SZ_FITS)
		namps = 0
		len = strlen (kwchar)
		for (ip = 1; ip <= len; ip = ip + 1) {
			if (kwchar[ip] == ',' || kwchar[ip] == '\042')
				kwchar[ip] = ' '
		}
		call sscan (kwchar)
		do i = 1, 4 {
			call gargi (ival)
			if (nscan() == i && ival > 0)
				namps = namps + 1
			else
				break
		}
	} else if (imaccf (im1, "IMTYPE") == YES) {
		call imgstr (im1, "IMTYPE", kwchar, SZ_FITS)
		if (streq (kwchar, "TWOAMPTOP")) {
			namps = 2
		} else {
			namps = 1
			call eprintf ("ONE AMP assumed (based on IMTYPE)!\n")
		}
	} else {
		call eprintf ("NUMAMPS missing; ONE AMP assumed!\n")
		namps = 1
	}

	if (imaccf (im1, "PREPIX") == YES) {
		prepix = imgetr (im1, "PREPIX")
	} else {
		prepix = clgeti ("pre_pix")
		call eprintf ("PREPIX missing; %2d from param file used!\n")
			call pargi (prepix)
	}

	if (imaccf (im1, "POSTPIX") == YES) {
		postpix = imgetr (im1, "POSTPIX")
	} else {
		postpix = clgeti ("post_pix")
		call eprintf ("POSTPIX missing; %2d from param file used!\n")
			call pargi (postpix)
	}

	if (prepix != 21) {
#		call eprintf ("Prepix != 21, NOT resetting!\n")
		call eprintf ("Prepix != 21, resetting!\n")
		prepix = 21
		postpix = postpix - 1
	}

	if (imaccf (im1, "BINNING") == YES) {
		call imgstr (im1, "BINNING", kwchar, SZ_FITS)
		len = strlen (kwchar)
		for (ip = 1; ip <= len; ip = ip + 1) {
			if (kwchar[ip] == ',' || kwchar[ip] == '\042')
				kwchar[ip] = ' '
		}
		call sscan (kwchar)
			call gargi (binx)
			call gargi (biny)
			if (nscan() < 2)
				call fatal (0, "Problem with decoding BINNING")
	} else {
		call eprintf ("BINNING missing; assume no binning!\n")
		binx = 1
		biny = 1
	}

# Now ready to calculate window:
# What happens to PREPIX, POSTPIX with binning?? Supposedly in BINNED pixels
	nx = 2048 / binx
	ny = 2048 / biny

# Check against other frames, if appropriate:
	if (flatten) {
		if (nx != nxf || ny != nyf)
			call fatal (0, "Incongruous images (flat)!\n")
	}
	if (bias_sub) {
		if (nx != nxb || ny != nyb)
			call fatal (0, "Incongruous images (bias)!\n")
	}

	x1 = prepix * namps + 1
	x2 = x1 + nx - 1
	nxlow = (x2 + x1 - 1) / 2
	y1 = 1
	y2 = ny
	osa1 = prepix * namps + nx + 1
	osa2 = prepix * namps + nx + postpix
	if (namps == 2) {
		osb1 = prepix * namps + nx + postpix + 1
		osb2 = prepix * namps + nx + postpix + postpix
	} else {
		osb1 = 0
		osb2 = 0
	}
if (osa2 != IM_LEN(im1,1) && osb2 != IM_LEN(im1,1)) {
call eprintf ("x sections don't add up! %d:%d  %d:%d+%d:%d\n")
call pargi (x1)
call pargi (x2)
call pargi (osa1)
call pargi (osa2)
call pargi (osb1)
call pargi (osb2)
}

if (binx > 1)
call eprintf ("Binning in x > 1 is untested!!\n")


# Allocate the bias buffers:
	call malloc (blev1, ny, TY_REAL)
	if (namps == 2)
		call malloc (blev2, ny, TY_REAL)

# Calculate the bias levels:
	bwid = osa2 - osa1 + 1
	if (mod (bwid, 2) != 1) {
		call eprintf ("WARNING: Overscan osa1 increased by 1\n")
		osa1 = osa1 + 1
		bwid = bwid - 1
	}
	bmid = (bwid + 1) / 2

	if (nord == 0) {
		l1 = clgeti ("l1")
		l2 = clgeti ("l2")
	} else {
		l1 = 1 + FIT_BK			# NB: these limits assumed below
		l2 = ny - FIT_BK
	}

	bufb = imgs2i (im1, osa1, osa2, l1, l2)

	bsum = 0.
	do i = 0, l2-l1 {
		bmed = asoki (Memi[bufb+i*bwid], bwid, bmid)
		Memr[blev1+i] = bmed
		bsum = bsum + bmed
	}
	biasa = bsum / (l2 - l1 + 1)

# Do the spline fit if requested
	if (nord > 0) {
		call malloc (ypix, ny, TY_REAL)
		call malloc (wts, ny, TY_REAL)
		do i = 0, ny-1
			Memr[ypix+i] = i + 1.
		call cvinit (cv, SPLINE3, 1, 1., real (ny))
		call cvfit (cv, Memr[ypix+FIT_BK], Memr[blev1], Memr[wts],
						ny-2*FIT_BK, WTS_UNIFORM, ier)
		if (ier != OK) {
			call eprintf ("Error in fit; average used\n")
			call amovkr (biasa, Memr[blev1], ny)
		}
		call cvvector (cv, Memr[ypix], Memr[blev1], ny)
	} else {
		call amovkr (biasa, Memr[blev1], ny)
	}


	if (namps == 2) {
		bwid = osb2 - osb1 + 1
		if (mod (bwid, 2) != 1) {
			call eprintf ("WARNING: Overscan osb1 increased by 1\n")
			osb1 = osb1 + 1
			bwid = bwid - 1
		}
		bmid = (bwid + 1) / 2
		bufb = imgs2i (im1, osb1, osb2, l1, l2)
		bsum = 0.
		do i = 0, l2-l1 {
			bmed = asoki (Memi[bufb+i*bwid], bwid, bmid)
			Memr[blev2+i] = bmed
			bsum = bsum + bmed
		}
		biasb = bsum / (l2 - l1 + 1)

	    if (nord > 0) {
		call malloc (ypix, ny, TY_REAL)
		call malloc (wts, ny, TY_REAL)
		do i = 0, ny-1
			Memr[ypix+i] = i + 1.
		call cvinit (cv, SPLINE3, 1, 1., real (ny))
		call cvfit (cv, Memr[ypix+FIT_BK], Memr[blev2], Memr[wts],
						ny-2*FIT_BK, WTS_UNIFORM, ier)
		if (ier != OK) {
			call eprintf ("Error in fit; average used\n")
			call amovkr (biasb, Memr[blev2], ny)
		}
		call cvvector (cv, Memr[ypix], Memr[blev2], ny)
	    } else {
		call amovkr (biasb, Memr[blev2], ny)
	    }
	}

# open output image and set pixel type:
	call xt_mkimtemp (image, imageout, imorig, SZ_FNAME) # ref pkg$xtools
	im3 = immap (imageout, NEW_COPY, im1)
	IM_LEN(im3,1) = nx
	IM_LEN(im3,2) = ny

# figure out scaling
	if (intscale) {
		IM_PIXTYPE(im3) = TY_SHORT
		scaling = SCALE
	} else {
		IM_PIXTYPE(im3) = TY_REAL
		scaling = 1.
	}

	fact1 = scaling
	if (namps == 2) {
		ratio = clgetr ("ratio2_1")
		fact1 = ratio / (0.5 * (1. + ratio)) * scaling
		fact2 = 1. / (0.5 * (1. + ratio)) * scaling
		if (intscale) {
		    if (fact1 > 0.5) {
			fact2 = fact2 * 0.5 / fact1
			fact1 = 0.5
			call eprintf ("inscale factor modified to < 0.5! \n")
		    } else if (fact2 > 0.5) {
			fact1 = fact1 * 0.5 / fact2
			fact2 = 0.5
			call eprintf ("inscale factor modified to < 0.5! \n")
		    }
		}
	}

# Add info to image header
	call sprintf (history, SZ_LINE, ": lproc bias%6.1f/%6.1f; bstrct=%s; flat=%s")
		call pargr (biasa)
		call pargr (biasb)
		call pargstr (biasim, SZ_FNAME)
		call pargstr (flatim, SZ_FNAME)

	if (intscale)
		call strcat ("; SCALED\n", history, SZ_LINE)
	else
		call strcat ("\n", history, SZ_LINE)
	call bksp_strcat (history, IM_HISTORY(im3), SZ_IMHIST)

# flag as processed:
	call imaddr (im3, "LPROCESS", 0.1)

# add ORIGNAME to header:
	call imastr (im3, "ORIGNAME", image, SZ_FITS)
	
# ... and DATASEC
	call sprintf (kwchar, SZ_FITS, "[%d:%d,%d:%d]")
		call pargi (x1)
		call pargi (x2)
		call pargi (y1)
		call pargi (y2)
	call imastr (im3, "DATASEC", kwchar, SZ_FITS)
	call printf ("DATASEC = %s \n")
		call pargstr (kwchar)

# ... and BIASSEC
	if (nord == 0)
		call sprintf (kwchar, SZ_FITS, "[%d:%d,%d:%d]")
	else
		call sprintf (kwchar, SZ_FITS, "[%d:%d,%d:%d] -- SPLINE3 FIT")
		call pargi (osa1)
		call pargi (osa2)
		call pargi (y1)
		call pargi (y2)
	call imastr (im3, "BIASSEC1", kwchar, SZ_FITS)
	call printf ("BIASSEC1 = %s \n")
		call pargstr (kwchar)

	if (namps == 2) {
		call sprintf (kwchar, SZ_FITS, "[%d:%d,%d:%d]")
			call pargi (osb1)
			call pargi (osb2)
			call pargi (y1)
			call pargi (y2)
		call imastr (im3, "BIASSEC2", kwchar, SZ_FITS)
		call printf ("BIASSEC2 = %s \n")
			call pargstr (kwchar)
	}

# ... and bias levels:
	call imaddr (im3, "BIASLEV1", biasa)
	if (namps == 2)
		call imaddr (im3, "BIASLEV2", biasb)

# ... and bias-structure info:
	if (bias_sub)
		call imastr (im3, "BIASIM", biasim, SZ_FITS)

# ... and scaling info:
	call imaddr (im3, "SCALFAC1", fact1)
	if (namps == 2)
		call imaddr (im3, "SCALFAC2", fact2)

# ... and flat info:
	if (flatten)
		call imastr (im3, "FLATIM", flatim, SZ_FITS)


# Setup for processing
	call amovkl (long(1), v1, IM_MAXDIM)
	call amovkl (long(1), v2, IM_MAXDIM)
	call amovkl (long(1), v3, IM_MAXDIM)
	call amovkl (long(1), v4, IM_MAXDIM)
	v1[2] = y1
	ndx = y1 - 1

# Do the processing:
# double-amp mode
	if (bias_sub) {			# Start  bias-sub if
	if (namps == 2) {
	  xoff = x1 - 1			# offset to first pixel
call eprintf ("SCALING: fact1,2 = %f %f\n")
call pargr (fact1)
call pargr (fact2)
	  if (flatten) {
	    while (impnlr(im3, line3, v3) != EOF ) {
		if (imgnlr (im1, line1, v1) == EOF)
			call eprintf (" READ ERROR (input)\n")
		if (imgnlr (im2, line2, v2) == EOF)
			call eprintf (" READ ERROR (flat)\n")
		if (imgnlr (im4, line4, v4) == EOF)
			call eprintf (" READ ERROR (bias)\n")
		do i = 0, nxlow-x1
			Memr[line3+i] = (Memr[line1+xoff+i] - Memr[blev1+ndx] -
					Memr[line4+i]) * fact1 / Memr[line2+i]
		do i = nxlow-xoff, nx-1
			Memr[line3+i] = (Memr[line1+xoff+i] - Memr[blev2+ndx] -
					Memr[line4+i]) * fact2 / Memr[line2+i]
		ndx = ndx + 1
	    }
	  } else {
	    while (impnlr(im3, line3, v3) != EOF ) {
		if (imgnlr (im1, line1, v1) == EOF)
			call eprintf (" READ ERROR (input)\n")
		if (imgnlr (im4, line4, v4) == EOF)
			call eprintf (" READ ERROR (bias)\n")
		do i = 0, nxlow-x1
			Memr[line3+i] = (Memr[line1+xoff+i] - Memr[blev1+ndx] -
							Memr[line4+i]) * fact1
		do i = nxlow-xoff, nx-1
			Memr[line3+i] = (Memr[line1+xoff+i] - Memr[blev2+ndx] -
							Memr[line4+i]) * fact2
		ndx = ndx + 1
	    }
	  }
	} else if (reverse) {
	  xoff = x1 + (nx - 1) - 1		# offset to first pixel
	  if (flatten) {
	    while (impnlr(im3, line3, v3) != EOF ) {
		if (imgnlr (im1, line1, v1) == EOF)
			call eprintf (" READ ERROR (input)\n")
		if (imgnlr (im2, line2, v2) == EOF)
			call eprintf (" READ ERROR (flat)\n")
		if (imgnlr (im4, line4, v4) == EOF)
			call eprintf (" READ ERROR (bias)\n")
		do i = 0, nx-1
			Memr[line3+i] = (Memr[line1+xoff-i] - Memr[blev1+ndx] -
					Memr[line4+i]) * fact1 / Memr[line2+i]
		ndx = ndx + 1
	    }
	  } else {
	    while (impnlr(im3, line3, v3) != EOF ) {
		if (imgnlr (im1, line1, v1) == EOF)
			call eprintf (" READ ERROR (input)\n")
		if (imgnlr (im4, line4, v4) == EOF)
			call eprintf (" READ ERROR (bias)\n")
# flip input line1 here
		do i = 0, nx-1
			Memr[line3+i] = (Memr[line1+xoff-i] - Memr[blev1+ndx] -
							Memr[line4+i]) * fact1
		ndx = ndx + 1
	    }
	  }
	} else {
# single-amp, non-flip
	  xoff = x1 - 1			# offset to first pixel
	  if (flatten) {
	    while (impnlr(im3, line3, v3) != EOF ) {
		if (imgnlr (im1, line1, v1) == EOF)
			call eprintf (" READ ERROR (input)\n")
		if (imgnlr (im2, line2, v2) == EOF)
			call eprintf (" READ ERROR (flat)\n")
		if (imgnlr (im4, line4, v4) == EOF)
			call eprintf (" READ ERROR (bias)\n")
		do i = 0, nx-1
			Memr[line3+i] = (Memr[line1+xoff+i] - Memr[blev1+ndx] -
					Memr[line4+i]) * fact1 / Memr[line2+i]
		ndx = ndx + 1
	    }
	  } else {
	    while (impnlr(im3, line3, v3) != EOF ) {
		if (imgnlr (im1, line1, v1) == EOF)
			call eprintf (" READ ERROR (input)\n")
		if (imgnlr (im4, line4, v4) == EOF)
			call eprintf (" READ ERROR (bias)\n")
		do i = 0, nx-1
			Memr[line3+i] = (Memr[line1+xoff+i] - Memr[blev1+ndx] -
							Memr[line4+i]) * fact1
		ndx = ndx + 1
	    }
	  }
	}

	} else {		# continue  bias-sub if

	if (namps == 2) {
	  xoff = x1 - 1			# offset to first pixel
call eprintf ("SCALING: fact1,2 = %f %f\n")
call pargr (fact1)
call pargr (fact2)
	  if (flatten) {
	    while (impnlr(im3, line3, v3) != EOF ) {
		if (imgnlr (im1, line1, v1) == EOF)
			call eprintf (" READ ERROR (input)\n")
		if (imgnlr (im2, line2, v2) == EOF)
			call eprintf (" READ ERROR (flat)\n")
		do i = 0, nxlow-x1
			Memr[line3+i] = (Memr[line1+xoff+i] - Memr[blev1+ndx]) *
							fact1 / Memr[line2+i]
		do i = nxlow-xoff, nx-1
			Memr[line3+i] = (Memr[line1+xoff+i] - Memr[blev2+ndx]) *
							fact2 / Memr[line2+i]
		ndx = ndx + 1
	    }
	  } else {
	    while (impnlr(im3, line3, v3) != EOF ) {
		if (imgnlr (im1, line1, v1) == EOF)
			call eprintf (" READ ERROR (input)\n")
		do i = 0, nxlow-x1
			Memr[line3+i] = (Memr[line1+xoff+i] - Memr[blev1+ndx]) *
									fact1
		do i = nxlow-xoff, nx-1
			Memr[line3+i] = (Memr[line1+xoff+i] - Memr[blev2+ndx]) *
									fact2
		ndx = ndx + 1
	    }
	  }
	} else if (reverse) {
	  xoff = x1 + (nx - 1) - 1		# offset to first pixel
	  if (flatten) {
	    while (impnlr(im3, line3, v3) != EOF ) {
		if (imgnlr (im1, line1, v1) == EOF)
			call eprintf (" READ ERROR (input)\n")
		if (imgnlr (im2, line2, v2) == EOF)
			call eprintf (" READ ERROR (flat)\n")
		do i = 0, nx-1
			Memr[line3+i] = (Memr[line1+xoff-i] - Memr[blev1+ndx]) *
							fact1 / Memr[line2+i]
		ndx = ndx + 1
	    }
	  } else {
	    while (impnlr(im3, line3, v3) != EOF ) {
		if (imgnlr (im1, line1, v1) == EOF)
			call eprintf (" READ ERROR (input)\n")
# flip input line1 here
		do i = 0, nx-1
			Memr[line3+i] = (Memr[line1+xoff-i] - Memr[blev1+ndx]) *
									fact1
		ndx = ndx + 1
	    }
	  }
	} else {
# single-amp, non-flip
	  xoff = x1 - 1			# offset to first pixel
	  if (flatten) {
	    while (impnlr(im3, line3, v3) != EOF ) {
		if (imgnlr (im1, line1, v1) == EOF)
			call eprintf (" READ ERROR (input)\n")
		if (imgnlr (im2, line2, v2) == EOF)
			call eprintf (" READ ERROR (flat)\n")
		do i = 0, nx-1
			Memr[line3+i] = (Memr[line1+xoff+i] - Memr[blev1+ndx]) *
							fact1 / Memr[line2+i]
		ndx = ndx + 1
	    }
	  } else {
	    while (impnlr(im3, line3, v3) != EOF ) {
		if (imgnlr (im1, line1, v1) == EOF)
			call eprintf (" READ ERROR (input)\n")
		do i = 0, nx-1
			Memr[line3+i] = (Memr[line1+xoff+i] - Memr[blev1+ndx]) *
									fact1
		ndx = ndx + 1
	    }
	  }
	}
	} # end bias-sub if

	if (nord > 0) {
		call cvfree (cv)
		call mfree (wts, TY_REAL)
		call mfree (ypix, TY_REAL)
	}

	call imunmap (im1)
	call imunmap (im3)
	call xt_delimtemp (imageout, imorig)

	if (namps == 2)
		call mfree (blev2, TY_REAL)
	call mfree (blev1, TY_REAL)

	if (bias_sub)
		call imunmap (im4)
	if (flatten)
		call imunmap (im2)
end

#############################################################
#############################################################

real	osa1[8,2]			# overscan area x1
real	osa2[8,2]			# overscan area x2
real	blev[8,2]			# constant bias level?
real	gain[8,2]			# fixed gain level


	if (imaccf (im1, "PREPIX") == YES) {
		prepix = imgetr (im1, "PREPIX")
	} else {
		prepix = clgeti ("pre_pix")
		call eprintf ("PREPIX missing; %2d from param file used!\n")
			call pargi (prepix)
	}

	if (imaccf (im1, "POSTPIX") == YES) {
		postpix = imgetr (im1, "POSTPIX")
	} else {
		postpix = clgeti ("post_pix")
		call eprintf ("POSTPIX missing; %2d from param file used!\n")
			call pargi (postpix)
	}

	if (imaccf (im1, "BINNING") == YES) {
		call imgstr (im1, "BINNING", kwchar, SZ_FITS)
		len = strlen (kwchar)
		for (ip = 1; ip <= len; ip = ip + 1) {
			if (kwchar[ip] == ',' || kwchar[ip] == '\042')
				kwchar[ip] = ' '
		}
		call sscan (kwchar)
			call gargi (binx)
			call gargi (biny)
			if (nscan() < 2)
				call fatal (0, "Problem with decoding BINNING")
	} else {
		call eprintf ("BINNING missing; assume no binning!\n")
		binx = 1
		biny = 1
	}
if (binx > 1)
call eprintf ("Binning in x > 1 is untested!!\n")


# Now ready to calculate window:
# What happens to PREPIX, POSTPIX with binning?? Supposedly in BINNED pixels
	nx = 2048 / binx
	ny = 2048 / biny

# Check against other frames, if appropriate:
	if (flatten) {
		if (nx != nxf || ny != nyf)
			call fatal (0, "Incongruous images (flat)!\n")
	}
	if (bias_sub) {
		if (nx != nxb || ny != nyb)
			call fatal (0, "Incongruous images (bias)!\n")
	}

if (namps == 2)
	call fatal (0, "Whoa-- 2 amp mode needs to be implemented!")





	x1 = prepix * namps + 1
	x2 = x1 + nx - 1
	nxlow = (x2 + x1 - 1) / 2
	y1 = 1
	y2 = ny
	osa1 = prepix * namps + nx + 1
	osa2 = prepix * namps + nx + postpix


#1. determine number of amps, chips, pre-postpix, binning etc
#2. determine geometry what are overscan regions, etc.
#3. fit bias level, fill vector
#4. do the actual vector ops:
#	fill const bias value (either single or from fit)
#	read image line
#	subtract from line (nxchip times)
#	subtract bias-struct image line
#	multi. by gain factor (nxchip times)
#	divide by flat image line
#	write out image line
#	switch at midpoint-y, if needed

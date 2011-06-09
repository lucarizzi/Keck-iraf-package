include	<imhdr.h>
include	<error.h>
 
# T_L2PROCESS: Do a 1-pass bias/trim/flat-field, aloowing for two amplifiers
#
#  Needs a check for bwidth, cskip inside bounds!

define		SCALE	0.5		# Scaling for use with SHORT format

procedure t_l2process()

char	image[SZ_FNAME]
char	flatim[SZ_FNAME]
char	imageout[SZ_FNAME], imorig[SZ_FNAME]
char	history[SZ_LINE]
bool	flatten					# flatten the image?
bool	intscale				# half and store as short?
int	x1, y1, x2, y2				# region of interest
int	l1, l2					# lines of bias strip
int	low1, low2				# col. of lower bias strip
int	upp1, upp2				# col. of upper bias strip
int	nxlow					# last x of lower half of CCD
pointer	im1, im2, im3
long	v1[IM_MAXDIM], v2[IM_MAXDIM], v3[IM_MAXDIM]
pointer	line1, line2, line3

int	nx, ny					# size of output image
real	bias1, bias2				# bias levels 1 and 2
int	xoff, i
int	bwid, bmid
real	bsum
real	scaling
pointer	bufb

bool	clgetb(), streq()
int	clgeti(), impnlr(), imgnlr(), asoki()
pointer	immap(), imgs2i()

begin
	call clgstr ("image", image, SZ_FNAME)
	call clgstr ("out_image", imageout, SZ_FNAME)
	flatten = clgetb ("flatten")
	intscale = clgetb ("intscale")

	x1 = clgeti ("x1")
	y1 = clgeti ("y1")
	x2 = clgeti ("x2")
	y2 = clgeti ("y2")
	nxlow = clgeti ("nxlow")
	l1 = clgeti ("l1")
	l2 = clgeti ("l2")
	low1 = clgeti ("low1")
	low2 = clgeti ("low2")
	upp1 = clgeti ("upp1")
	upp2 = clgeti ("upp2")

	nx = x2 - x1 + 1
	ny = y2 - y1 + 1

	im1 = immap (image, READ_ONLY, 0)

# Calculate the bias levels:
	bwid = low2 - low1 + 1
	if (mod (bwid, 2) != 1) {
		call eprintf ("WARNING: low1 decreased by 1\n")
		low1 = low1 - 1
		bwid = bwid + 1
	}
	bmid = (bwid + 1) / 2

	bufb = imgs2i (im1, low1, low2, l1, l2)
	bsum = 0.
	do i = 0, l2-l1 {
		bsum = bsum + asoki (Memi[bufb+i*bwid], bwid, bmid)
	}
	bias1 = bsum / (l2 - l1 + 1)

	bwid = upp2 - upp1 + 1
	if (mod (bwid, 2) != 1) {
		call eprintf ("WARNING: upp2 increased by 1\n")
		upp2 = upp2 + 1
		bwid = bwid + 1
	}
	bmid = (bwid + 1) / 2
	bufb = imgs2i (im1, upp1, upp2, l1, l2)
	bsum = 0.
	do i = 0, l2-l1 {
		bsum = bsum + asoki (Memi[bufb+i*bwid], bwid, bmid)
	}
	bias2 = bsum / (l2 - l1 + 1)

	if (flatten) {
		call clgstr ("flatim", flatim, SZ_FNAME)
		if (streq (flatim, "")) {
			call eprintf ("WARNING: not flattened!\n")
			flatten = false
		} else {
			im2 = immap (flatim, READ_ONLY, 0)

			if (nx != IM_LEN(im2,1) || ny != IM_LEN(im2,2))
				call fatal ("Incongruous image sizes!\n")
		}
	}

# open output image and write
	call xt_mkimtemp (image, imageout, imorig, SZ_FNAME) # ref pkg$xtools
	im3 = immap (imageout, NEW_COPY, im1)
	IM_LEN(im3,1) = nx
	IM_LEN(im3,2) = ny
	if (intscale) {
		IM_PIXTYPE(im3) = TY_SHORT
		scaling = SCALE
	} else {
		IM_PIXTYPE(im3) = TY_REAL
		scaling = 1.
	}

# add info to image header
	if (flatten) {
	    call sprintf (history, SZ_LINE, ": lproc bias%6.1f/%6.1f; flat=%s")
		call pargr (bias1)
		call pargr (bias2)
		call pargstr (flatim, SZ_FNAME)
	} else {
	    call sprintf (history, SZ_LINE, ": lproc bias%6.1f/%6.1f; NO FLAT")
		call pargr (bias1)
		call pargr (bias2)
	}
	if (intscale)
		call strcat ("; SCALED\n", history, SZ_LINE)
	else
		call strcat ("\n", history, SZ_LINE)
	call bksp_strcat (history, IM_HISTORY(im3), SZ_IMHIST)

# do the copy
	call amovkl (long(1), v1, IM_MAXDIM)
	call amovkl (long(1), v2, IM_MAXDIM)
	call amovkl (long(1), v3, IM_MAXDIM)
	v1[2] = y1

	xoff = x1 - 1
	if (flatten) {
	    while (impnlr(im3, line3, v3) != EOF ) {
		if (imgnlr (im1, line1, v1) == EOF)
			call eprintf (" READ ERROR (input)\n")
		if (imgnlr (im2, line2, v2) == EOF)
			call eprintf (" READ ERROR (flat)\n")
		do i = 0, nxlow-x1
			Memr[line3+i] = (Memr[line1+xoff+i] - bias1) *
							scaling / Memr[line2+i]
		do i = nxlow-xoff, nx-1
			Memr[line3+i] = (Memr[line1+xoff+i] - bias2) *
							scaling / Memr[line2+i]
	    }
	} else {
	    while (impnlr(im3, line3, v3) != EOF ) {
		if (imgnlr (im1, line1, v1) == EOF)
			call eprintf (" READ ERROR (input)\n")
		do i = 0, nxlow-x1
			Memr[line3+i] = (Memr[line1+xoff+i] - bias1) * scaling
		do i = nxlow-xoff, nx-1
			Memr[line3+i] = (Memr[line1+xoff+i] - bias2) * scaling
	    }
	}

	if (flatten)
		call imunmap (im2)
	call imunmap (im1)
	call imunmap (im3)
	call xt_delimtemp (imageout, imorig)
end

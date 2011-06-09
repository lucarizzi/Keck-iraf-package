include	<imhdr.h>
include <error.h>

#
# MKMOSAIC: Pack a list of subrasters into a larger image.
# Still some shortcomings -- eg no test for overflow at top

# THERE IS SOME KIND OF BUG; the imgs2r(im ...) seems to help ...
#

procedure t_mkmosaic()

char	image[SZ_FNAME]			# name of mosaic image
char	list[SZ_FNAME]			# name of list inp,x,y,sky,scl,nx,ny
char	def_image[SZ_FNAME]		# default input image
bool	mk_image			# create a mosaic image?
bool	verbose				# print out actions?
bool	nice_form			# truncate image to look nice?
int	ncol,nline			# dimensions of mosaic image
real	background			# background value (including border)
real	fill				# fill value
bool	stipple				# fill empty areas with stipple?
int	border				# width of border
char	title[SZ_IMTITLE]		# title string for output mosaic
int	xga1, xga2, yga1, yga2		# limits of optically active region
int	nx, ny				# default size of subraster
real	frac1, frac2			# low, upp fract to reject in autoscale
real	normflux			# flux value in output with auto-scaling
bool	auto_scale			# should the output be auto-scaled?

pointer	fda				# file descriptor for input list
pointer	im, ima				# image, subrast descriptors
				# For Each Subraster:
char	sr_image[SZ_FNAME]		# name of input image for subraster
real	xcen, ycen			# centers of fields
real	sky				# sky to be subtracted (if any)
real	scl				# scaling factor (if any)
int	i, j				# "window" position
int	szx, szy			# size of subraster (boxes)
int	xa1, xa2, ya1, ya2		# limits of optically active region
real	cnts				# summed counts in sky-sub. subraster

int	mx, my				# size of subraster (pixels)
int	ix, iy				# actual center of subraster copied
int	nbx, nby			# number of windows in x, y
int	ndx				# index to subraster position
int	k, narg, toggle
int	xoff, yoff, noff, moff, npix, mpix, lpix
int	i1, i2, j1, j2
int	ia1, ia2, ja1, ja2
int	skip

long	v[IM_MAXDIM]
pointer	buffill
pointer	bufin, bufout

# These variables for multi-HDU files
#char	ccdprop[SZ_FNAME]	# file name of CCD properties
pointer	im0			# image pointer to primary HDU
pointer	mmap			# pointer to mosaic vectors
int	nccd
# int	nx, ny, i1, i2, j1, j2
# pointer	im2
pointer	chip_sect()

bool	clgetb(), streq()
int	fscan(), nscan()
int	clgeti()
int	impnlr()
pointer	imgs2r(), imps2r()
pointer	immap()
pointer	open()
real	clgetr()
errchk	immap()

begin
	call clgstr ("mosaic", image, SZ_FNAME)
	call clgstr ("list", list, SZ_FNAME)
	call clgstr ("def_image", def_image, SZ_FNAME)
	border = clgeti ("border")
	background = clgetr ("background")
	fill = clgetr ("fill")
	stipple = clgetb ("stipple")
	nx = clgeti ("nx")
	xga1 = clgeti("xa1")
	xga2 = clgeti("xa2")
	yga1 = clgeti("ya1")
	yga2 = clgeti("ya2")
	ny = clgeti ("ny")
	auto_scale = clgetb ("auto_scale")
	normflux = clgetr ("normflux")
	frac1 = clgetr ("f1")
	frac2 = clgetr ("f2")
	verbose = clgetb ("verbose")
	mk_image = clgetb ("mk_image")
	if (mk_image) {
## Test for nonexistence (existence below)
		ncol = clgeti ("ncols")
		nline = clgeti ("nlines")
		call clgstr ("title", title, SZ_IMTITLE)
		nice_form = clgetb ("nice_form")
		if (nice_form) {
			nbx = ncol / nx
			nby = nline / ny
			ncol = nbx * nx
			nline = nby * ny
			call eprintf ("Size: %d x %d, for %d x %d subrasters\n")
				call pargi (ncol)
				call pargi (nline)
				call pargi (nbx)
				call pargi (nby)
		}
		
		im = immap (image, NEW_IMAGE, 0)
		IM_NDIM(im) = 2
		IM_LEN(im,1) = ncol
		IM_LEN(im,2) = nline
		IM_PIXTYPE(im) = TY_REAL
		call strcpy (title, IM_TITLE(im), SZ_IMTITLE)   
		call amovkl (long(1), v, IM_MAXDIM)
#	fill with background
		while ( impnlr(im, bufout, v) != EOF )
			call amovkr (background, Memr[bufout], ncol)
		call imunmap (im)
	}

	im = immap (image, READ_WRITE, 0)
	ncol = IM_LEN(im,1)
	nline = IM_LEN(im,2)
call eprintf ("nc,nl=%d,%d\n")
call pargi (ncol)
call pargi (nline)

# Open file containing list of images to copy

	if (list[1] == EOS) {
		call eprintf ("no coord. list\n")
		call erract (EA_FATAL)
	}

	iferr (fda = open (list, READ_ONLY, TEXT_FILE)) {
		call eprintf ("open error on list file\n")
		call erract (EA_FATAL)
	}

# Make stipple pattern
	call calloc (buffill, ncol+2, TY_REAL)
	if (stipple) {
		do k = 1, ncol, 4
			Memr[buffill+k] = fill
	} else {
		do k = 0, ncol+1
			Memr[buffill+k] = fill
	}

# Work out mosaic geometry
	nbx = ncol / nx
	nby = nline / ny

#
call eprintf ("Need to research and fix problem; see source\n")
# Loop through and copy subrasters
	ndx = 1
	while (fscan (fda) != EOF) {
		call gargwrd (sr_image, SZ_FNAME)
		call gargr (xcen)
		call gargr (ycen)
		call gargr (sky)
		call gargr (scl)
		call gargi (i)
		call gargi (j)
		call gargi (szx)
		call gargi (szy)
		call gargi (xa1)
		call gargi (xa2)
		call gargi (ya1)
		call gargi (ya2)
		narg = nscan()

#	ckeck arg list, fill in defaults as needed
		if (narg < 3) {
		    call eprintf ("poorly formatted list, line skipped\n")
		    ndx = ndx + 1
		    next
		}
#	open input image
		if (streq (sr_image, ""))
			call strcpy (def_image, sr_image, SZ_FNAME)
#		iferr (ima = immap (sr_image, READ_ONLY, 0)) {
#			call eprintf ("Cannot open image %s -- skipped!\n")
#				call pargstr (sr_image, SZ_FNAME)
#			ndx = ndx + 1
#			next
#		}
# Open mosaic image:
		call mos_init (srimage, "deimos$test/prop.txt", im0, mmap, nccd)
		if (narg < 4)
			sky = 0.
		if (narg < 5)
			scl = 1.
		if (narg == 6) {
		    call eprintf ("Incomplete position, line skipped\n")
		    ndx = ndx + 1
		    next
		}
		if (narg < 7) {		# put into sequential position
			i = mod ((ndx-1), nbx) + 1
			j = (ndx - 1) / nbx + 1
		}
		if (narg < 8)
			szx = 1
		if (narg < 9)
			szy = 1
		if (narg < 10)
			xa1 = xga1
		if (narg < 11)
			xa2 = xga2
		if (narg < 10)
			ya1 = max (1, yga1)
		if (narg < 10)
			ya2 = yga2

		if (xcen == INDEF || ycen == INDEF) {
			skip = YES
			xcen = 1
			ycen = 1
		} else {
			skip = NO
		}

		ix = xcen + 0.5
		iy = ycen + 0.5
		mx = szx * nx
		my = szy * ny
		xa2 = min (IM_LEN(ima,1), xga2)
		ya2 = min (IM_LEN(ima,2), yga2)

		if (j > nby || j < 1 || i > nbx || i < 1) {
		    call eprintf ("window position outside bounds, skipped\n")
		    ndx = ndx + 1
		    next
		}

# Copy the subraster:
#	work out coords, etc; recall nx, ny, mx, my INCLUDE border width
		ix = xcen + 0.5
		noff = (mx - 1) / 2
		i1 = ix - noff + border/2
		i2 = i1 + (mx - border) - 1
		iy = ycen + 0.5
		noff = (my - 1) / 2
		j1 = iy - noff + border/2
		j2 = j1 + (my - border) - 1
#	checks for xa1, etc:
		ia1 = max (min (i1, xa2), xa1)
		ia2 = max (min (i2, xa2), xa1)
		ja1 = max (min (j1, ya2), ya1)
		ja2 = max (min (j2, ya2), ya1)
		xoff = ia1 - i1
		yoff = ja1 - j1
#	get subraster, scale
		bufin = imgs2r (ima, ia1, ia2, ja1, ja2)
		npix = (ia2 - ia1 + 1) * (ja2 - ja1 + 1)
		if (auto_scale) {
			call autoscl (Memr[bufin], npix, frac1,frac2, sky, cnts)
			scl = normflux / cnts
			call printf ("Sky,flux: %6f %6f\n")
				call pargr (sky)
				call pargr (cnts)
		}
		call aaddkr (Memr[bufin], -sky, Memr[bufin], npix)
		call amulkr (Memr[bufin],  scl, Memr[bufin], npix)
		
#	set up the output buffer, fill with stipple or const.
		i1 = (i-1) * nx + 1 + border/2
		i2 = i1 + (mx - border) - 1
		j1 = (nby-j) * ny + 1 + border/2	# top down
#		j1 = (j-1) * ny + 1 + border/2
		j2 = j1 + (my - border) - 1
		i1 = max (i1, 1)
		i2 = min (i2, ncol)
		j1 = max (j1, 1)
		j2 = min (j2, nline)
		npix = i2 - i1 + 1
		bufout = imps2r (im, i1, i2, j1, j2)
		noff = 0
		do k = j1, j2 {
		    toggle = mod (k, 2) * 2
		    call amovr (Memr[buffill+toggle], Memr[bufout+noff], npix)
		    noff = noff + npix
		}

#	do the copy
		if (skip == NO) {
		    mpix = ia2 - ia1 + 1
		    lpix = min (mpix, npix-xoff)
		    noff = yoff * npix + xoff
		    moff = 0
		    do k = ja1, ja2 {
		    	call amovr (Memr[bufin+moff], Memr[bufout+noff], lpix)
		    	moff = moff + mpix
		    	noff = noff + npix
		    }
		}
		
		if (verbose) {
			call printf ("Copying %s[%d:%d,%d:%d] at %d,%d\n")
				call pargstr (sr_image, SZ_FNAME)
				call pargi (ia1)
				call pargi (ia2)
				call pargi (ja1)
				call pargi (ja2)
				call pargi (i1+xoff)
				call pargi (j1+yoff)
			call eprintf ("")
		}
#	unmap the input image, increment index and go on
		call imunmap (ima)
		ndx = ((j - 1) * nbx) + i + 1
### HERE IS PROBLEM; need 2 buffers, not just one
#		bufout = imgs2r (im, i1, i2, j1, j2)
	}

# Put border around left edge and bottom
	if (border >= 1) {
		npix = border * nline
		bufout = imps2r (im, 1, border, 1, nline)
		call amovkr (background, Memr[bufout], npix)
		npix = border * ncol
		bufout = imps2r (im, 1, ncol, 1, border)
		call amovkr (background, Memr[bufout], npix)
	}
call eprintf ("DEBUG: looks, closeup %d %d %d\n")
call pargi (im)
call pargi (bufout)
call pargi (fda)

# Close up
	call mfree (buffill, TY_REAL)
	call close (fda)
call eprintf (" other closes OK -- IMUNMAP??\n")
	call imunmap (im)
end


procedure autoscl (vect, n, frac1, frac2, sky, flux)

real	vect[n]
int	n
real	frac1, frac2
real	sky, flux

int	i1, ilen, midpt
pointer	buf1

real	asokr(), asumr()
begin
	call malloc (buf1, n, TY_REAL)

	call asrtr (vect, Memr[buf1], n)
	i1 = (n-1) * frac1
	ilen = n * (frac2 - frac1)

	midpt = ilen / 2
	sky = asokr (Memr[buf1+i1], ilen, midpt)
	call aaddkr (vect, -sky, Memr[buf1], n)
	flux = asumr (Memr[buf1], n)
# call eprintf ("Sky level = %5f (%d--%d)\n")
# call pargr (sky)
# call pargi (i1)
# call pargi (ilen)

	call mfree (buf1, TY_REAL)
end

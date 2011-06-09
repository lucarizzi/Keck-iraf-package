include	<imhdr.h>
include <error.h>

define	DEF_PROP_FILE	"deimos$lib/prop/redprop.txt"	# default properties

#
# MKMONTAGE: Sample subrasters, in this case from a multi-HDU image.
# Hacked-up version of mkmosaic


procedure t_mkmontage()

char	in_image[SZ_FNAME]		# default input image
char	list[SZ_FNAME]			# name of list x,y,sky,scl,nx,ny
char	out_image[SZ_FNAME]		# name of output mosaic image
bool	mk_image			# create a mosaic image?
bool	verbose				# print out actions?
bool	nice_form			# truncate image to look nice?
int	ncol,nline			# dimensions of mosaic image
real	background			# background value (including border)
real	fill				# fill value
bool	stipple				# fill empty areas with stipple?
int	border				# width of border
char	title[SZ_IMTITLE]		# title string for output mosaic
int	nx, ny				# default size of subraster
real	frac1, frac2			# low, upp fract to reject in autoscale
real	normflux			# flux value in output with auto-scaling
bool	auto_scale			# should the output be auto-scaled?

pointer	fda				# file descriptor for input list
pointer	im				# output image descriptors
				# For Each Subraster:
real	xcen, ycen			# centers of fields
real	sky				# sky to be subtracted (if any)
real	scl				# scaling factor (if any)
int	i, j				# "window" position
int	szx, szy			# size of subraster (boxes)
real	cnts				# summed counts in sky-sub. subraster

int	mx, my				# size of subraster (pixels)
int	ix, iy				# actual center of subraster copied
int	nbx, nby			# number of windows in x, y
int	ndx				# index to subraster position
int	k, narg, toggle
int	xoff, yoff, noff, moff, npix, mpix, lpix
int	i1, i2, j1, j2
int	skip

long	v[IM_MAXDIM]
pointer	buffill
pointer	bufin, bufout

# These variables for multi-HDU files
#char	ccdprop[SZ_FNAME]	# file name of CCD properties
pointer	im0			# image pointer to primary HDU
pointer	mmap			# pointer to mosaic vectors
int	nccd
int	chip			# chip number for coords
# int	nx, ny, i1, i2, j1, j2
int	x1, x2, y1, y2
# pointer	im2
pointer	chip_sect()
int	det_chip()

bool	clgetb()
int	fscan(), nscan()
int	clgeti()
int	impnlr()
pointer	imps2r()
pointer	immap()
pointer	open()
real	clgetr()

begin
	call clgstr ("mosaic", out_image, SZ_FNAME)
	call clgstr ("list", list, SZ_FNAME)
	call clgstr ("in_image", in_image, SZ_FNAME)
	border = clgeti ("border")
	background = clgetr ("background")
	fill = clgetr ("fill")
	stipple = clgetb ("stipple")
	nx = clgeti ("nx")
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
		
		im = immap (out_image, NEW_IMAGE, 0)
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

	im = immap (out_image, READ_WRITE, 0)
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

call eprintf ("DEBUG1\n")

# Open input image:
 	call mos_init (in_image, DEF_PROP_FILE, im0, mmap, nccd)

call eprintf ("DEBUG2\n")
#
call eprintf ("Need to research and fix problems; see source\n")
# Loop through and copy subrasters
	ndx = 1
	while (fscan (fda) != EOF) {
		call gargr (xcen)
		call gargr (ycen)
		call gargr (sky)
		call gargr (scl)
		call gargi (i)
		call gargi (j)
		call gargi (szx)
		call gargi (szy)
		narg = nscan()

#	ckeck arg list, fill in defaults as needed
		if (narg < 2) {
		    call eprintf ("poorly formatted list, line skipped\n")
		    ndx = ndx + 1
		    next
		}
		if (narg < 3)
			sky = 0.
		if (narg < 4)
			scl = 1.
		if (narg == 5) {
		    call eprintf ("Incomplete position, line skipped\n")
		    ndx = ndx + 1
		    next
		}
		if (narg < 6) {		# put into sequential position
			i = mod ((ndx-1), nbx) + 1
			j = (ndx - 1) / nbx + 1
		}
		if (narg < 7)
			szx = 1
		if (narg < 8)
			szy = 1

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
		xoff = 0
		yoff = 0
#	get subraster, scale
# Get the image section (limit checking and zero-fill included)
		x1 = i1
		x2 = i2
		y1 = j1
		y2 = j2
		chip = det_chip (xcen, ycen)
		bufin = chip_sect (i1, i2, j1, j2, mmap, chip, YES)

		npix = (i2 - i1 + 1) * (j2 - j1 + 1)
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
		    mpix = i2 - i1 + 1
		    lpix = min (mpix, npix-xoff)
		    noff = yoff * npix + xoff
		    moff = 0
		    do k = j1, j2 {
		    	call amovr (Memr[bufin+moff], Memr[bufout+noff], lpix)
		    	moff = moff + mpix
		    	noff = noff + npix
		    }
		}
		
		if (verbose) {
			call eprintf ("Copying [%d:%d,%d:%d](%d) at %d,%d\n")
				call pargi (x1)
				call pargi (x2)
				call pargi (y1)
				call pargi (y2)
				call pargi (chip)
				call pargi (i1+xoff)
				call pargi (j1+yoff)
			call eprintf ("\n")
		}
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
	call mos_free (im0, mmap)
call eprintf (" other closes OK -- IMUNMAP??\n")
	call imunmap (im)
call eprintf ("ends ...\n")

end



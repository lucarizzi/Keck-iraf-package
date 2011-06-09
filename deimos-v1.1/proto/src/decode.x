# DECODE: Attempt to unpack MDs IDL structures into FITS images
# This requires the CFITSIO library, and makes use of BJW's similar
# FORTRAN routine.
####### Definitions for 

# Define the struct for the target data
define	DATLEN		7
define	PTFLUX		Memi[$1]
define	PTIVAR		Memi[$1+1]
define	PTMASK		Memi[$1+2]
define	PTWAVE		Memi[$1+3]
define	PTLCOEF		Memi[$1+4]
define	PTTCOEF		Memi[$1+5]
define	PTDLAM		Memi[$1+6]

# Note that the number of targets can also be described here.

define	FLUX		Memr[PTFLUX($1)+$2]		# flux
define	IVAR		Memr[PTIVAR($1)+$2]		# inverse variance
define	MASK		Mems[PTMASK($1)+$2]		# mask
define	WAVE		Memr[PTWAVE($1)+$2]		# wave (A)
define	LCOEF		Memd[PTLCOEF($1)+$2]		# 
define	TCOEF		Memd[PTTCOEF($1)+$2]		# 
define	DLAM		Memr[PTDLAM($1)+$2]		# 

define	MAXDIM		2		# maximum table-entry array dimensions

## NB: local "fitsio.h" comes from all the defines in cfitsio.h (?) with
## appropriate characters replaced.  XXX NOTE IN MAINTENANCE DOC

###  %	include "deimos$cfitsio/f77.inc"

include <imhdr.h>
include	<math/curfit.h>
include "fitsio.h"

procedure	t_decode ()

char	input[SZ_FNAME]			# Input bintable
char	output[SZ_FNAME]		# output test file
int	hdu				# number of hdu to decode

int	nx, ny
pointer	dat
pointer	dim			# Pointer to length array
pointer	cv1, cv2

pointer	im
# long	v[IM_MAXDIM]

double	x
int	i, ndx
real	y, yc, w		# TMP
real	w0, dw

double	dcveval()
int	clgeti()
pointer	immap(), impl2r(), impl1r()
begin
	call clgstr ("input", input, SZ_FNAME)
	hdu = clgeti ("extn")
	call clgstr ("output", output, SZ_FNAME)

	if (clgeti ("ndim") == 2) {
		call decode2 (input, hdu, dat, dim, nx, ny, cv1, cv2)

		call eprintf ("Input, size: %s: %dx %d\n")
			call pargstr (input)
			call pargi (nx)
			call pargi (ny)
	yc = 0.5 * (ny+1)
	x = 1.
	y = 1.
	w = dcveval (cv1, x)# * (1. + (y - yc) * dcveval (cv2, x))+ DLAM(dat,0)
	call eprintf ("First pix wave = %7.2f + %7.2f\n")
		call pargr (w)
		call pargr (DLAM(dat,0))

	x = 4096.
	y = ny
	w = dcveval (cv1, x)# * (1. + (y - yc) * dcveval (cv2, x))+ DLAM(dat,ny-1)
	call eprintf ("Last  pix wave = %7.2f + %7.2f\n")
		call pargr (w)
		call pargr (DLAM(dat,ny-1))

		im = immap (output, NEW_IMAGE, 0)
		IM_NDIM(im) = 2
		IM_LEN(im,1) = nx
		IM_LEN(im,2) = ny
		IM_PIXTYPE(im) = TY_REAL

		ndx = 0
		do i = 1, ny {
			call amovr (FLUX(dat,ndx), Memr[impl2r(im, i)], nx)
			ndx = ndx + nx
		}
		call imunmap (im)

	} else {

		call decode1 (input, hdu, dat, dim, nx, ny, cv1, cv2)

		call eprintf ("Input, size: %s: %dx %d\n")
			call pargstr (input)
			call pargi (nx)
			call pargi (ny)

		im = immap (output, NEW_IMAGE, 0)

		dw = (WAVE(dat,nx-1) - WAVE(dat,0)) / 4095.	# XXX TEMP
		w0 = WAVE(dat,0)				# XXX TEMP
		w0 = w0 - dw * (w0/dw - int (w0/dw))

		call spec_out (im, dat, nx, w0, dw)

		call imunmap (im)
	}

end

	

#
# DECODE1:  decode a UCB 1-D bintable
#

procedure	decode1 (fname, hdu, dat, dim, nx, ny)

char	fname[ARB]		# FITS BINTABLE file name
int	hdu			# HDU/extn number
pointer	dat			# Pointer to array structure
pointer	dim			# Pointer to length array
int	nx, ny			# Size of array

char	xname[SZ_FNAME]		# FITS name with extention
%	character*80 f77nam
# %	character*80 f77lin
# %	character*80 f77nul

int	lu
int	stat
int	nc, nr
int	n1, n2

int	col_decode()
begin
# Read in desired table: build name
	call sprintf (xname, SZ_FNAME, "%-s[%d]")
		call pargstr (fname)
		call pargi (hdu)
	call f77pak (xname, f77nam, SZ_FNAME)

	stat = 0
# open the file and check format
	call ftgiou (lu, stat)
	call ftnopn (lu, f77nam, READONLY, stat)   # XXX ftopen

	call ftgncl (lu, nc, stat)
	call ftgnrw (lu, nr, stat)

	if (nc != 8 || nr != 1) {
		call eprintf ("%s, nc=%d nr=%d stat=%d\n")
		call pargstr (xname)
		call pargi (nc)
		call pargi (nr)
		call pargi (stat)
		call fatal (0, "Unexpected format in input FITS table")
	}

# OK, we seem to have a real file here; let's try to unpack it

# First, allocate vectors
	call malloc (dat, DATLEN, TY_STRUCT)
	call malloc (dim, DATLEN, TY_INT)

# Unpack FLUX, LAMBDA, IVAR, (MASK),
	stat = col_decode (lu, dat, TY_REAL, PTFLUX(dat), PTFLUX(dim), 1, n1, n2)
	nx = n1
	ny = n2

	stat = col_decode (lu, dat, TY_REAL, PTWAVE(dat), PTWAVE(dim), 2, n1, n2)
	stat = col_decode (lu, dat, TY_REAL, PTIVAR(dat), PTIVAR(dim), 3, n1, n2)

	call ftclos (lu, stat)
	call ftfiou (lu, stat)

end

#
# DECODE2:  decode a UCB "SLITS" bintable
#

procedure	decode2 (fname, hdu, dat, dim, nx, ny, cv1, cv2)

char	fname[ARB]		# FITS BINTABLE file name
int	hdu			# HDU/extn number
pointer	dat			# Pointer to array structure
pointer	dim			# Pointer to length array
int	nx, ny			# Size of array
pointer	cv1, cv2		# curfit pointers; probably should go elsewhere

char	xname[SZ_FNAME]		# FITS name with extention
%	character*80 f77nam
# %	character*80 f77lin
# %	character*80 f77nul

int	lu
int	stat
int	nc, nr
int	n1, n2

int	col_decode()
begin
# Read in desired table: build name
	call sprintf (xname, SZ_FNAME, "%-s[%d]")
		call pargstr (fname)
		call pargi (hdu)
	call f77pak (xname, f77nam, SZ_FNAME)

	stat = 0
# open the file and check format
	call ftgiou (lu, stat)
	call ftnopn (lu, f77nam, READONLY, stat)   # XXX ftopen

	call ftgncl (lu, nc, stat)
	call ftgnrw (lu, nr, stat)

	if (nc != 10 || nr != 1) {
		call eprintf ("%s, nc=%d nr=%d stat=%d\n")
		call pargstr (xname)
		call pargi (nc)
		call pargi (nr)
		call pargi (stat)
		call fatal (0, "Unexpected format in input FITS table(2)")
	}

#	call f77pak ("???", f77nul, 80)
#	call ftgcvs (lu, 1, 1, 1, 1, f77nul, f77lin, isnul, stat)
#	call f77upk (f77lin, DESNAME(indat), 80)
#	call ftgcvs (lu, 2, 1, 1, 1, f77nul, f77lin, isnul, stat)
#	call f77upk (f77lin, PROJNAME(indat), 80)

# OK, we seem to have a real file here; let's try to unpack it

# First, allocate vectors
	call malloc (dat, DATLEN, TY_STRUCT)
	call malloc (dim, DATLEN, TY_INT)

# Unpack FLUX, IVAR, MASK, .., LCOEF, TCOEF, .., DLAM
	stat = col_decode (lu, dat, TY_REAL, PTFLUX(dat), PTFLUX(dim), 1, n1, n2)
	nx = n1
	ny = n2

	stat = col_decode (lu, dat, TY_REAL, PTIVAR(dat), PTIVAR(dim), 2, n1, n2)
	stat = col_decode (lu, dat, TY_CHAR,  PTMASK(dat), PTMASK(dim), 3, n1, n2)
	stat = col_decode (lu, dat, TY_DOUBLE, PTLCOEF(dat), PTLCOEF(dim), 7, n1, n2)
	stat = col_decode (lu, dat, TY_DOUBLE, PTTCOEF(dat), PTTCOEF(dim), 8, n1, n2)
	stat = col_decode (lu, dat, TY_REAL, PTDLAM(dat), PTDLAM(dim), 10, n1, n2)
	

# set up the polynomial evaluations
	nc = PTLCOEF(dim)
	call dcvset (cv1, LEGENDRE, 1.D0, 4096.D0, LCOEF(dat,0), nc)

	nc = PTTCOEF(dim)
	call dcvset (cv2, LEGENDRE, 1.D0, 4096.D0, TCOEF(dat,0), nc)


	call ftclos (lu, stat)
	call ftfiou (lu, stat)

end

#
# This is the general routine for decoding and unpacking a table-element.  It
# is currently resticted to 2-D, primarily by the listing of only 2 axes in args
#

int	procedure col_decode (lu, dat, data_type, ptdat, ptdim, icol, n1, n2)

int	lu
pointer	dat
int	data_type
pointer	ptdat
int	ptdim
int	icol
int	n1, n2

int	data_code
int	nax, naxes[MAXDIM]
int	n, nn, nw
int	stat

bool	nuls			# Null value(s) present?

begin
	stat = 0
	call amovki (1, naxes, MAXDIM)
	call ftgtdm (lu, icol, MAXDIM, nax, naxes, stat)
	n1 = naxes[1]
	n2 = naxes[2]
	n = n1 * n2
	call malloc (ptdat, n, data_type)
	ptdim = n

# Unpack vector; first, get type:
	call ftgtcl (lu, icol, data_code, nn, nw, stat)

call eprintf ("DEBUG: icol=%d, %d x %d; data_code=%d\n")
call pargi (icol)
call pargi (n1)
call pargi (n2)
call pargi (data_code)
	if (n != nn)
		call fatal (0, "COL_DECODE: dimension mismatch!")
	
# unpack (note: we are _not_ supporting NULL values here; use FTGCVt for that)
# also, note that INT assumes long

	switch (data_code) {
	    case 11:
		if (data_type != TY_CHAR)
			call fatal (0, "COL_DECODE: data_type mismatch!")
		call ftgcvi (lu, icol, 1, 1, n, 0, Memc[ptdat], nuls, stat)

	    case 21:
		if (data_type != TY_SHORT)
			call fatal (0, "COL_DECODE: data_type mismatch!")
		call ftgcvi (lu, icol, 1, 1, n, INDEFS, Mems[ptdat], nuls, stat)

	    case 41:
		if (data_type != TY_INT)
			call fatal (0, "COL_DECODE: data_type mismatch!")
		call ftgcvj (lu, icol, 1, 1, n, INDEFI, Memi[ptdat], nuls, stat)

	    case 42:
		if (data_type != TY_REAL)
			call fatal (0, "COL_DECODE: data_type mismatch!")
		call ftgcve (lu, icol, 1, 1, n, INDEFR, Memr[ptdat], nuls, stat)

	    case 82:
		if (data_type != TY_DOUBLE)
			call fatal (0, "COL_DECODE: data_type mismatch!")
		call ftgcvd (lu, icol, 1, 1, n, INDEFD, Memd[ptdat], nuls, stat)

	    default:
		call fatal (0, "COL_DECODE: unsupported data_code")
	}

	if (nuls) {
		call eprintf ("NULLS PRESENT -- NOT SUPPORTED?")
	}
	if (stat != 0) {
		call eprintf ("ERROR ON READ: stat=%d\n")
		call pargi (stat)
	}

	return (stat)
end


#  Currently HARD-CODED to assume reals

procedure	spec_out (im, dat, nvect, w0, dw)

pointer	im			# pointer to output image
pointer	dat			# pointer to data structure
int	nvect			# length of input
real	w0, dw

int	nx
int	i
pointer	mbuf, wbuf

pointer	impl3r()
begin
	nx = (WAVE(dat,nvect-1) - w0) / dw + 1

# Allocate working buffers
	call malloc (mbuf, nx, TY_REAL)
	call malloc (wbuf, nx, TY_REAL)

# Set up image
	IM_NDIM(im) = 3				# XXX
	IM_LEN(im,1) = nx			# XXX
	IM_LEN(im,2) = 1			# XXX
	IM_LEN(im,3) = 5			# XXX
	IM_PIXTYPE(im) = TY_REAL

# Add relevant keywords to output header
	call imastr (im, "CTYPE1", "LINEAR")
	call imaddr (im, "CRPIX1", 1.)
	call imaddr (im, "CRVAL1", w0)
	call imaddr (im, "CD1_1", dw)		# Marc uses angstroms
	call imaddr (im, "CDELT1", dw)

# Axis 2 wcs should be OK, add axis 3
#	call imastr (im, "WCSDIM", 3)
#	call imastr (im, "CTYPE3", "LINEAR")
#	call imaddr (im, "CRPIX3", 1.)
#	call imaddr (im, "CRVAL3", 1.)
#	call imaddr (im, "CDELT3", 1.)
#	call imaddr (im, "CD3_3", 1.)
#	call imaddr (im, "LTM3_3", 1.)

# Change system for splot
	call imastr (im, "WAT0_001", "system=equispec")

call eprintf ("no flux conservation applied\n")

	call lin_spec (WAVE(dat,0), FLUX(dat,0), nvect, Memr[impl3r (im, 1, 1)], w0, dw, nx, NO)
#	call amovkr (1., Memr[impl3r (im, 1, 2)], nx)
#	call amovkr (0., Memr[impl3r (im, 1, 3)], nx)
#	call lin_spec (WAVE(dat,0), IVAR(dat,0), nvect, Memr[impl3r (im, 1, 4)], w0, dw, nx, NO)

	call lin_spec (WAVE(dat,0), IVAR(dat,0), nvect, Memr[wbuf], w0, dw, nx, NO)
	do i = 0, nx-1 {
		if (Memr[wbuf+i] <= 1.e-9) {
			Memr[wbuf+i] = 1.e6
			Memr[mbuf+i] = 0.
		} else {
			Memr[wbuf+i] = 1. / Memr[wbuf+i] 
			Memr[mbuf+i] = 1.
		}
	}
	call amovr  (Memr[mbuf], Memr[impl3r (im, 1, 2)], nx)
	call amovkr (0., Memr[impl3r (im, 1, 3)], nx)
	call amovr  (Memr[wbuf], Memr[impl3r (im, 1, 4)], nx)
	call amovr  (Memr[wbuf], Memr[impl3r (im, 1, 5)], nx)

	call mfree (wbuf, TY_REAL)
	call mfree (mbuf, TY_REAL)

end


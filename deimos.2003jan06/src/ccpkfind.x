# CCPKFIND: Find the peak of 2-D CC function

include <imhdr.h>

procedure t_ccpkfind ()

char	image[SZ_FNAME]			# input image
int	xref, yref		# x,y of corner for xc region
pointer	im

int	ned			# edge around pk
int	nx, ny
int	i, j
int	ipk, jpk
int	inext, jnext		# second highest
real	xoff, yoff
real	maxval
pointer	buflin, bufsect

int	clgeti()
pointer	immap(), imgs2r(), imgl2r()

begin

# Read in parameters:
	call clgstr ("image", image, SZ_FNAME)
	xref = clgeti ("xref")
	yref = clgeti ("yref")
	ned = clgeti ("border")

	im = immap (image, READ_ONLY, 0)

	nx = IM_LEN(im,1)
	ny = IM_LEN(im,2)

	maxval = 0.

	do j = 1, ny {
		buflin = imgl2r (im, j)
		do i = 0, nx-1 {
			if (Memr[buflin+i] > maxval) {
				inext = ipk
				jnext = jpk
				ipk = i + 1
				jpk = j
				maxval = Memr[buflin+i]
			}
		}
	}

# OK, we now have peak pixel; extract box around peak
	bufsect = imgs2r (im, ipk-ned, ipk+ned, jpk-ned, jpk+ned)

# fit parabola to pk, using weighted averages
	xoff = 0.
	yoff = 0.

	call ppf (Memr[bufsect], ned*2+1, xoff, yoff)
	call ppf (Memr[bufsect], ned*2+1, xoff, yoff)

#	call printf ("%20s %8.2f %8.2f\n")
	call printf ("%8.2f %8.2f\n")
#		call pargstr (image)
		call pargr (ipk+xoff+xref-1)
		call pargr (jpk+yoff+yref-1)
	
	if (abs (ipk-inext) > 1 || abs (jpk-jnext) > 1)
		call printf ("Poorly defined peak!\n")
	call imunmap (im)
end

procedure	ppf (zarray, n, xoff, yoff)

real	zarray[n,n]
int	n
real	xoff, yoff

real	z[3]		# 1-D work array
real	wt[3]		# 1-D weights array

int	i, j
int	nc
int	iv, jv

begin
	nc = n / 2 + 1

	wt[1] = (0.5 - xoff)
	wt[2] = 1.
	wt[3] = (0.5 + xoff)

	do j = -1, 1 {
	    jv = j + 2
	    z[jv] = 0.
	    do i = -1, 1 {
		z[jv] = z[jv] + wt[i+2]*zarray[nc+i,nc+j]
	    }
	}

	yoff = 0.5 * (z[1] - z[3]) / (z[1] + z[3] - 2.*z[2])



	wt[1] = (0.5 - yoff)
	wt[2] = 1.
	wt[3] = (0.5 + yoff)

	do i = -1, 1 {
	    iv = i + 2
	    z[iv] = 0.
	    do j = -1, 1 {
		z[iv] = z[iv] + wt[j+2]*zarray[nc+i,nc+j]
	    }
	}
	
	xoff = 0.5 * (z[1] - z[3]) / (z[1] + z[3] - 2.*z[2])
end
	

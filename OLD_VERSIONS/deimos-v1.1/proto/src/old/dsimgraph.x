#
# These are the graphical routines for DSIMULATOR
#

include	<math.h>
include	<gset.h>
include	<gim.h>
include	"dsimulator.h"
include	"deimos.h"

define	FAST	10.
define	MEDIUM	1.
define	SLOW	0.1
define	LBL1	"-X"
define	LBL2	"+Y"
define	LBL3	"Seen From Front"

procedure	select_targ (tdat, ntarg, indat)

pointer	tdat
int	ntarg
pointer	indat

pointer	mdat			# pointer to data structure for mask outlines

#### XXX TOTALLY TMP!!  -- used for drawing slits only
pointer	sdat
int	nslit

int	i
int	nlist			# number of the active list

real	anginc, xyinc, speed
real	xoff, yoff, aoff, doff
real	cosp, sinp		# could be double
real	cosdec			# could be double

char	command[32]			# not sure if 32 is good
int	wcs, key
real	wx, wy
pointer	gp

pointer	gopen()
int	clgcur(), get_sw_nearest0()

begin
	nlist = 1		## XXX TMP DEBUG

# Read in the mask outlines
	call get_outlines (mdat)

# Open the graphics stream
	gp = gopen ("stdgraph", NEW_FILE, STDGRAPH)

	call fp_layout (gp, mdat, indat, tdat, ntarg, nlist, LBL1, LBL2, LBL3)

# start the interactive stuff

	speed = FAST
	anginc = DEGTORAD(1.) * speed		# TMP_HARDCODE (was PA_INCR)
	xyinc = 6. * speed		# TMP_HARDCODE (was TR_INCR)

	cosdec = cos (DEC_FLD(indat))
	cosp = cos (PA_ROT(indat))
	sinp = sin (PA_ROT(indat))

	while ( clgcur("coord", wx, wy, wcs, key, command, 32) != EOF ) {

	    if (key == 'q') {
		break
	    }

	    switch (key) {			# NB note TWO switch cases!

# x = cosp
# y = sinp
# translations  XXX there's a lot of excess calc's going on here

		case 'h':
			xoff = DEGTORAD(xyinc/3600.)
			yoff = 0.
			doff =  xoff * cosp - yoff * sinp
			aoff =  xoff * sinp + yoff * cosp
			RA_FLD(indat) = RA_FLD(indat) + aoff / cosdec
			DEC_FLD(indat) = DEC_FLD(indat) + doff
			
		case 'j':
			xoff = 0.
			yoff = DEGTORAD(xyinc/3600.)
			doff =  xoff * cosp - yoff * sinp
			aoff =  xoff * sinp + yoff * cosp
			RA_FLD(indat) = RA_FLD(indat) + aoff / cosdec
			DEC_FLD(indat) = DEC_FLD(indat) + doff
			
		case 'k':
			xoff = 0.
			yoff = -DEGTORAD(xyinc/3600.)
			doff =  xoff * cosp - yoff * sinp
			aoff =  xoff * sinp + yoff * cosp
			RA_FLD(indat) = RA_FLD(indat) + aoff / cosdec
			DEC_FLD(indat) = DEC_FLD(indat) + doff
			
		case 'l':
			xoff = -DEGTORAD(xyinc/3600.)
			yoff = 0.
			doff =  xoff * cosp - yoff * sinp
			aoff =  xoff * sinp + yoff * cosp
			RA_FLD(indat) = RA_FLD(indat) + aoff / cosdec
			DEC_FLD(indat) = DEC_FLD(indat) + doff
			
		case 'p':
			PA_ROT(indat) = PA_ROT(indat) + anginc

		case 'n':
			PA_ROT(indat) = PA_ROT(indat) - anginc

		case 's':
			do nlist = 1, 4
				call selector (tdat, ntarg, nlist, 8.)	# XXX
			nlist = 1
			do i = 1, ntarg
				call mark_obj (gp, tdat, i, nlist)

		case 'a':
			i = get_sw_nearest0 (gp, XARCS(tdat,0), YARCS(tdat,0),
				SEL(tdat,0), ntarg, NO, wx, wy, wcs)
			SEL(tdat,i) = YES
			call mark_obj (gp, tdat, i, SAMPL(tdat,i))

		case 'd':
			i = get_sw_nearest0 (gp, XARCS(tdat,0), YARCS(tdat,0),
				SEL(tdat,0), ntarg, YES, wx, wy, wcs)
			SEL(tdat,i) = NO
			call mark_obj (gp, tdat, i, nlist)

		case 'i':
			call amovkr (NO, SEL(tdat,0), ntarg)
# XXX and call slit_free()

	    }

	    switch (key) {

		case 'h','j','k','l','p','n':
			cosdec = cos (DEC_FLD(indat))
			cosp = cos (PA_ROT(indat))
			sinp = sin (PA_ROT(indat))
		    call fld2telax (indat)
		    call tel_coords (tdat, ntarg, indat)

	    }

	    switch (key) {

		case 'r','c','h','j','k','l','p','n','f','i','w':
			call fp_layout (gp, mdat, indat, tdat, ntarg, nlist, LBL1, LBL2, LBL3)

		case '.':
			if (speed == FAST) {
				speed = MEDIUM
				call printf (" medium")
			} else if (speed == MEDIUM) {
				speed = SLOW
				call printf (" slow")
			} else if (speed == SLOW) {
				speed = FAST
				call printf (" fast")
			}
			anginc = DEGTORAD(1.) * speed	# TMP -- see above
			xyinc = 6. * speed	# TMP -- see above

## TMP XXX
		case 'x':
			call gen_slits (tdat, ntarg, sdat, nslit, indat)
#			call len_slits (tdat, ntarg, sdat, nslit, indat)
			do i = 0, nslit-1
				call mark_slit (gp, sdat, i)

		case '?':
		    call gpagefile (gp, KEYSFILE, "simulator cursor commands")

		case 'I':
			call fatal (0, "INTERRUPT")
	    }
	}

	call gclose (gp)
end

#
# GET_OUTLINES: read in the outlines of mask, TV in Focal plane
# 

procedure	get_outlines (mdat)

pointer	mdat				# pointer to Foc Plane struct

char	tchar
int	nact
int	ndx
pointer	fdf

int	line_count()
int	fscan(), nscan()
pointer	open()
begin
	call malloc (mdat, MDATLEN, TY_STRUCT)

	fdf = open (FP_FILE, READ_ONLY, TEXT_FILE)

	nact = line_count (fdf)

	NFPDES(mdat) = nact
	call malloc (PTFPLX(mdat), nact, TY_REAL)
	call malloc (PTFPLY(mdat), nact, TY_REAL)
	call malloc (PTFPLZ(mdat), nact, TY_INT)

# read in the descriptions (x,y,pen)
	ndx = 0
	while (fscan (fdf) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#' || nscan() == 0) {
			next
		}
		call reset_scan()
		call gargr (FPLX(mdat,ndx))
		call gargr (FPLY(mdat,ndx))
		call gargi (FPLZ(mdat,ndx))
		if (nscan() < 3)
			FPLZ(mdat,ndx) = -1
		else if (nscan() < 2)
			call fatal (0, "Problem with foc_plane file!")
		ndx = ndx + 1
	}

# do ndx = 0, nact-1 {
# call eprintf (" %6.0f %6.0f  %2d\n")
# call pargr (FPLX(mdat,ndx))
# call pargr (FPLY(mdat,ndx))
# call pargi (FPLZ(mdat,ndx))
# }

	call close (fdf)
end

#
# FP_LAYOUT: graphically describe the Focal Plane layout
#
	

procedure	fp_layout (gp, mdat, indat, tdat, ntarg, nlist, lab1, lab2, lab3)

pointer	gp
pointer	mdat
pointer	indat
pointer	tdat
int	ntarg
int	nlist			# number of the "active" list
char	lab1[ARB], lab2[ARB], lab3[ARB]

int	i
int	pen				# code for "pen" action

real	gx1, gx2, gy1, gy2
real	gszx, gszy
real	hbox				# box half-width (arcsec)
real	x0, y0, x, y
real	cosp, sinp

real	aspect

real	ggetr()
# pointer	gopen()
begin
#	hbox = PL_FWID(limit) / 2.
	hbox = 1100. / 2.

	gx1 = -hbox / 1.	# PL_XMAG(limit)	# these are for LongSlit
	gx2 =  hbox / 1.	# PL_XMAG(limit)
	gy1 = -hbox*0.8
	gy2 =  hbox*1.2
	gszx = gx2 - gx1
	gszy = gx2 - gx1
	

# gp = gopen ("stdgraph", NEW_FILE, STDGRAPH)
# get the aspect ratio for the device, so that squares are square.
	aspect = ggetr (gp, "ar")

	gy1 = gy1 * aspect
	gy2 = gy2 * aspect

# set windows
	call gclear (gp)
#	call gsview (gp, (0.5-0.49*aspect), (0.5+0.49*aspect), 0.01, 0.99)
	call gswind (gp, gx1, gx2, gy1, gy2)

# Draw the FP boundaries
	call gseti (gp, G_PLCOLOR, YELLOW)
	do i = 0, NFPDES(mdat)-1 {
		pen = FPLZ(mdat,i)
		if (pen == 0) {
			call gamove (gp, FPLX(mdat,i), FPLY(mdat,i))
			next
		}
		if (pen > 0) {
			if (pen == 1)
				call gseti (gp, G_PLTYPE, GL_SOLID)
			else if (pen == 2)
				call gseti (gp, G_PLTYPE, GL_DASHED)
			else if (pen == 3)
				call gseti (gp, G_PLTYPE, GL_DOTTED)
		}
		call gadraw (gp, FPLX(mdat,i), FPLY(mdat,i))
	}
# put in axes and labels
	call gseti (gp, G_PLTYPE, GL_SOLID)
	x0 = 0.12 * gszx + gx1
	y0 = 0.28 * gszy + gy1
	call gseti (gp, G_PLCOLOR, YELLOW)
	call gamove (gp, x0, y0)
	x =  0.05*gszy 
	y =  0
	call grdraw (gp, x, y)
	call gtext (gp, x0+1.2*x, y0, lab1, "h=l;v=c;q=h;s=0.9")
	call gamove (gp, x0, y0)
	x =  0
	y =  0.05*gszy 
	call grdraw (gp, x, y)
	call gtext (gp, x0, y0+1.2*y, lab2, "h=c;v=c;q=h;s=0.9")

	call gtext (gp, x0, y0-0.2*y, lab3, "h=c;v=t;q=h;s=0.9")


# Now mark the slits
	call gseti (gp, G_PLTYPE, GL_SOLID)

	do i = 0, ntarg-1
		call mark_obj (gp, tdat, i, nlist)



# Auxilliary plots:

# Draw the compass rose:
	x0 = 0.33 * gszx + gx1
	y0 = 0.23 * gszy + gy1
	call gseti (gp, G_PLCOLOR, RED)
	call gamove (gp, x0, y0)
	cosp = cos (PAR_ANG(indat)-PA_ROT(indat)) 
	sinp = sin (PAR_ANG(indat)-PA_ROT(indat))
	x =  0.07*gszy * cosp / 1.		# PL_XMAG(limit)
	y =  0.07*gszy * sinp
	call grdraw (gp, x, y)

	call gseti (gp, G_PLCOLOR, WHITE)
	call gamove (gp, x0, y0)
	cosp = cos (-PA_ROT(indat)) 
	sinp = sin (-PA_ROT(indat))
	x =  0.07*gszy * cosp / 1.		# PL_XMAG(limit)
	y =  0.07*gszy * sinp
	call grdraw (gp, x, y)
	call gtext (gp, x0+1.4*x, y0+1.4*y, "N", "h=c;v=c;q=h;s=0.7")
	call gamove (gp, x0, y0)
	x = -0.07*gszy * sinp / 1.		# PL_XMAG(limit)
	y =  0.07*gszy * cosp
	call grdraw (gp, x, y)


end

#
# MARK_OBJ: mark an object with relevant details
#

procedure	mark_obj (gp, tdat, i, nlist)

pointer	gp
pointer	tdat
int	i
int	nlist			# number of the "active" list

real	x, y, a, l1, l2
real	sina, cosa
begin
	if (STAT(tdat,i) == YES) {
	    if (SAMPL(tdat,i) == nlist)
		call gseti (gp, G_PLCOLOR, GREEN)
	    else if (PCODE(tdat,i) == CODE_AS)		# TMP! XXX
		call gseti (gp, G_PLCOLOR, MAGENTA)	# TMP! XXX
	    else
		call gseti (gp, G_PLCOLOR, BLUE)
	} else {
		call gseti (gp, G_PLCOLOR, RED)
	}


	if (SEL(tdat,i) == YES)
		call gseti (gp, G_PLCOLOR, WHITE)	# TMP XXX TEST

	x = XARCS(tdat,i)
	y = YARCS(tdat,i)
	a = RELPA(tdat,i)
	l1 = LEN1(tdat,i)	# XXX Need to work out the sense of these
	l2 = LEN2(tdat,i)
	if (a == INDEF) {
		sina = 0.
		cosa = 1.
	} else {
		sina = sin (a)
		cosa = cos (a)
	}

	if (PCODE(tdat,i) == CODE_AS) {
		call gmark (gp, x, y, GM_BOX, -2.*l1, -2.*l1)	# XXX
	} else {
		call gamove (gp, x-l1*cosa, y-l1*sina)
		call gadraw (gp, x+l2*cosa, y+l2*sina)
		call gmark (gp, x, y, GM_PLUS, -1.*SLWID(tdat,i), -1.*SLWID(tdat,i))	# TMP
	}
end

#
# MARK_SLIT: mark an object with relevant details
#

procedure	mark_slit (gp, sdat, i)

pointer	gp
pointer	sdat
int	i

real	x, y, a, l1, l2
real	yoff
real	sina, cosa
begin
	call gseti (gp, G_PLCOLOR, CYAN)

	yoff = 0.5*SLWID(sdat,i)

	x = X1(sdat,i)
	y = Y1(sdat,i)
	call gamove (gp, x, y+yoff)
	call gadraw (gp, x, y-yoff)

	x = X2(sdat,i)
	y = Y2(sdat,i)
	call gadraw (gp, x, y-yoff)
	call gadraw (gp, x, y+yoff)

	x = X1(sdat,i)
	y = Y1(sdat,i)
	call gadraw (gp, x, y+yoff)

end

#
# GET_SW_NEAREST0: get nearest point with appropriate value of an int switch
# "0" indicates zero-indexed
#

int procedure get_sw_nearest0 (gp, xdata, ydata, isw, ndata, sw_val, wx, wy, wcs)

pointer	gp
real 	xdata[ARB], ydata[ARB]
int	isw[ARB]
int	ndata
int	sw_val
real	wx, wy
int	wcs

int	nearest, i
real	ycorr
real	rsq, rsq_min
real	xndc, yndc, xgcndc, ygcndc

real	ggetr()

begin

# need to put in INDEF check 

# Get aspect ratio 
	ycorr = ggetr (gp, "ar")
	if (ycorr == 0.)
		ycorr = 1.

	rsq_min = 2.			# by def'n larger than NDC possible
	nearest = 0

	call gctran (gp, wx, wy, xgcndc, ygcndc, wcs, 0)

	if (sw_val == INDEF) {
	    do i = 1, ndata {
		call gctran (gp, xdata[i], ydata[i], xndc, yndc, wcs, 0)
		rsq = (xndc - xgcndc) ** 2 + ( (yndc - ygcndc) * ycorr) ** 2
		if (rsq < rsq_min) {
			rsq_min = rsq
			nearest = i
		}
	    }
	} else {
	    do i = 1, ndata {
		if (isw[i] != sw_val)
			next
		call gctran (gp, xdata[i], ydata[i], xndc, yndc, wcs, 0)
		rsq = (xndc - xgcndc) ** 2 + ( (yndc - ygcndc) * ycorr) ** 2
		if (rsq < rsq_min) {
			rsq_min = rsq
			nearest = i
		}
	    }
	}

	if (nearest != 0)
		call gscur (gp, xdata[nearest], ydata[nearest])
	else
		call eprintf ("no appropriate points")
	
	return (nearest - 1)		# Zero indexed
end

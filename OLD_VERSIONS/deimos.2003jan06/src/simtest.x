# SIMTEST: A special version of simulator for designing TEST slitmasks
# in this case, CENTER is used to refer to the telescope axis.
# There is no check for overlapping slits

include	<math.h>
include	<gset.h>
include	<gim.h>
include	"deimos.h"
include	"dsimulator.h"

procedure	t_simtest()

char	objfile[SZ_FNAME]			# ID, RA, Dec, equ, prior,
char	output[SZ_FNAME]			# output file name
char	mdf[SZ_FNAME]				# mask design file name
char	plotfile[SZ_FNAME]			# optional plotfile
char	maskid[SZ_FNAME]			# mask name
bool	surfcam					# ascii surfcam format?

pointer	fda, fdc				# in/out file descriptors
pointer	fdp			# pointer to optional plot file

pointer	indat			# pointer to instr/tel data structure

int	ntarg
pointer	tdat			# pointer to data structure for targets

int	nslit
pointer	sdat			# pointer to data structure for slits

int	i	## TMP
real	offset, stdwid	## TMP?
pointer	mdat	## TMP?

bool	clgetb(), strne()
pointer	open()

begin
	call clgstr ("objfile", objfile, SZ_FNAME)
#	call clgstr ("output", output, SZ_FNAME)
	call clgstr ("mdf", mdf, SZ_FNAME)
	call clgstr ("plotfile", plotfile, SZ_FNAME)
	call clgstr ("maskid", maskid, SZ_FNAME)
	surfcam = clgetb ("surfcam")

# Read in telescope, default data:
	call data_init (indat)

# Open the input list of targets:
        fda = open (objfile, READ_ONLY, TEXT_FILE)

        fdc = open (mdf, NEW_FILE, TEXT_FILE)

# Read in the target data:
#	call targ_init (fda, tdat, ntarg, indat)

# FOR GRID, USE THIS LINE INSTEAD of above AND COMMENT OUT call tel_coords
 	call grid_init (fda, tdat, ntarg, indat)
#	call calg_init (fda, tdat, ntarg, indat)

# ADD SPECIAL FOR COHU MASK:
#	call eprintf ("SPECIALS FOR COHU MASK!\n")
#	do i = 9, 14 {
#		SLWID(tdat,i) = SLWID(tdat,i) + 0.5
#	}
#	do i = 15, 20 {
#		SLWID(tdat,i) = SLWID(tdat,i) + 1.0
#	}
#	do i = 21, 26 {
#		SLWID(tdat,i) = SLWID(tdat,i) + 1.5
#	}

## ADD SPECIAL FOR LOH1
#	call eprintf ("SPECIALS FOR LOH1 MASK!\n")
#	SLWID(tdat,0) = SLWID(tdat,0) + 0.5



# Assign CENTER to telescope axis
	RA_TEL(indat) =  RA_FLD(INDAT)
	DEC_TEL(indat) = DEC_FLD(INDAT)
	PA_ROT(indat) = 0.

# Calculate position wrt telescope axis:
#	call tel_coords (tdat, ntarg, indat)

# Turn on ALL stat flags (special case)
	call amovki (YES, STAT(tdat,0), ntarg)


## ADD SPECIAL FOR COLLISION TEST (shift YARCS by fixed amount)
#	call eprintf ("SPECIALS FOR COLLISION TEST!\n")
#	call aaddkr (YARCS(tdat,0),   0., YARCS(tdat,0), ntarg)


# Ready to display and be interactive
	call select_targ (tdat, ntarg, indat)

# Save output?

# Calculate the SLITS in celestial coordinates.
# It is here that slit lengths, edges must be resolved.  Does binding go here?
# Or is this a step that should actually be included above?
#

	call gen_slits (tdat, ntarg, sdat, nslit, indat)
#	call len_slits (tdat, ntarg, sdat, nslit, indat)
# TMP
call eprintf ("nslit=%d\n")
call pargi (nslit)

call eprintf (" calling sky_coord \n")
# Convert adopted slits into coordinated on the sky:
	call sky_coords (sdat, nslit, indat)


## For symmetry, we now recalculate the telescope coords and then mask coords.
## In the end, this operation may be performed at Keck immediately prior to
## mask fabrication.
#  Note that _refraction_ correction should probably go here, as a differential
#  correction, as it is tiny and specific to actual time and date of observation

call eprintf (" calling tel_coord \n")

# Calculate position wrt telescope axis:
	call tel_coords (sdat, nslit, indat)

 
 call eprintf (" calling mask_coords \n")
# Generate the mask coordinates for the slits.

 	call mask_coords (sdat, nslit)
 

# ADD SPECIAL FOR COHU MASK:
# NB: must be extremely careful with offset indices -- eg (i+2) would probably
# work, i+2 does not!!
#
#	call eprintf ("SPECIALS FOR COHU MASK!\n")
#	stdwid = YMM3(sdat,1) - YMM2(sdat,1)
#	do i = 9, 24, 3 {
#		offset = 0.5 * (YMM3(sdat,i) - YMM2(sdat,i) - stdwid)
#		YMM1(sdat,i) = YMM1(sdat,i) - offset
#		YMM2(sdat,i) = YMM2(sdat,i) - offset
#		YMM3(sdat,i) = YMM3(sdat,i) - offset
#		YMM4(sdat,i) = YMM4(sdat,i) - offset
#	}
#	do i = 11, 26, 3 {
#		offset = 0.5 * (YMM3(sdat,i) - YMM2(sdat,i) - stdwid)
#		YMM1(sdat,i) = YMM1(sdat,i) + offset
#		YMM2(sdat,i) = YMM2(sdat,i) + offset
#		YMM3(sdat,i) = YMM3(sdat,i) + offset
#		YMM4(sdat,i) = YMM4(sdat,i) + offset
#	}


 
    if (!(surfcam)) {
	call eprintf (" calling write_design ...\n")
	call write_design (indat, tdat, ntarg, sdat, nslit, fdc, objfile)

    } else {
	call eprintf (" ASCII Surfcam Format ...\n")
	call write_surfcam (sdat, nslit, fdc, maskid)
## SPECIAL FOR GOH w/ round holes
#	call write_surfcam_holes (sdat, nslit, fdc, maskid)
    }

	if (strne (plotfile, "")) {
            fdp = open (plotfile, NEW_FILE, TEXT_FILE)
	    call write_mongo (sdat, nslit, fdp, maskid, plotfile)
	    call close (fdp)
	}

	call close (fdc)



#
call eprintf (" [recall: still must check proper x,y and tan mapping]\n")
#
end


#
# GRID_INIT: initialize data structure for targets; fill
# 

procedure	grid_init (fda, tdat, ntarg, indat)

pointer	fda		# file descriptor; file has zptx,y and spacing
pointer	tdat
int	ntarg
pointer	indat

char	idstr[SZ_ID]
double	alpha, delta
double	zptx, zpty, xoff, yoff

int	ndx
int	i, j
int	j1, j2, i1, i2
char	jid

int	fscan(), nscan()

begin
	call eprintf ("Calling grid_init; SPECIAL FOR GRID\n")

# Count the entries
	i = fscan (fda)
	call gargd (zptx)
	call gargd (zpty)
	call gargd (xoff)
	call gargd (yoff)
	if (nscan() < 4)
		call fatal (0, "problem reading input file!!")
call eprintf ("%f %f %f %f\n")
call pargd (zptx)
call pargd (zpty)
call pargd (xoff)
call pargd (yoff)
	i1 = (-500.0 - zptx) / xoff
	i2 = (500.0 - zptx) / xoff
	j1 = (180.0 - zpty) / yoff
	j2 = (480.0 - zpty) / yoff
call eprintf ("%d %d %d %d\n")
call pargi (i1)
call pargi (i2)
call pargi (j1)
call pargi (j2)
	ntarg = (abs (j2 - j1) + 1) * (abs (i2 - i1) + 1)

# Allocate the vectors
	call malloc (tdat, TDATLEN, TY_STRUCT)
	call malloc (PTINDEX(tdat), ntarg, TY_INT)
	call malloc (PTRA(tdat), ntarg, TY_DOUBLE)
	call malloc (PTDEC(tdat), ntarg, TY_DOUBLE)
	call malloc (PTPA(tdat), ntarg, TY_REAL)
	call malloc (PTLEN1(tdat), ntarg, TY_REAL)
	call malloc (PTLEN2(tdat), ntarg, TY_REAL)
	call malloc (PTWID(tdat), ntarg, TY_REAL)
	call malloc (PTPCODE(tdat), ntarg, TY_INT)
	call malloc (PTSAMPL(tdat), ntarg, TY_INT)

	call malloc (PTSTAT(tdat), ntarg, TY_INT)
	call malloc (PTSEL(tdat), ntarg, TY_INT)
	call malloc (PTXARCS(tdat), ntarg, TY_REAL)
	call malloc (PTYARCS(tdat), ntarg, TY_REAL)
	call malloc (PTRELPA(tdat), ntarg, TY_REAL)

# Allocate slit-index -- cross-reference -- NOTE the zeroing of this array
	call calloc (PTSLNDX(tdat), ntarg, TY_INT)

# Allocate the string array
	call malloc (PTLINE(tdat), ntarg*SZ_LINE, TY_CHAR)


# Generate the data
	ndx = 0
	do j = j1, j2 {
		jid = 'a'
		do i = i1, i2 {
			call sprintf (idstr, SZ_ID, "h%1s%03d")
				call pargc (jid)
				call pargi (i)
			alpha = 0.
			delta = 0.
			INDEX(tdat,ndx) = ndx
			RA(tdat,ndx) = DEGTORAD(alpha*15.)
			DEC(tdat,ndx) = DEGTORAD(delta)
			call strcpy (idstr, DATLINE(tdat,ndx), SZ_LINE-1)
			PA(tdat,ndx) = INDEF
			PCODE(tdat,ndx) = 100
			LEN1(tdat,ndx) = 0.25	# DEF_HLEN(indat)
			LEN2(tdat,ndx) = 0.25	# DEF_HLEN(indat)
			SLWID(tdat,ndx) = 0.5	# DEF_SLWID(indat)
			SAMPL(tdat,ndx) = PRIMARY
			SEL(tdat,ndx) = 1
			RELPA(tdat,ndx) = INDEF
			XARCS(tdat,ndx) = i * xoff + zptx
			YARCS(tdat,ndx) = j * yoff + zpty
			if (sqrt (XARCS(tdat,ndx)**2+YARCS(tdat,ndx)**2) > 600.)
				next
# check for 10th
			if ((j/10.)-(j/10)==0. || (i/10.)-(i/10)==0.) {
				SLWID(tdat,ndx) = SLWID(tdat,ndx) * 1.5
				LEN1(tdat,ndx) = LEN1(tdat,ndx) * 1.5
				LEN2(tdat,ndx) = LEN2(tdat,ndx) * 1.5
			}
			ndx = ndx + 1
		}
		jid = jid + 1
	}
	ntarg = ndx

call eprintf ("ntarg=%d\n")
call pargi (ntarg)

end

#
# CALG_INIT: initialize data structure for CALGRID; fill
# 

procedure	calg_init (fda, tdat, ntarg, indat)

pointer	fda		# file descriptor; file has zptx,y and spacing
pointer	tdat
int	ntarg
pointer	indat

char	idstr[SZ_ID]
double	alpha, delta
double	zptx, zpty, xoff, yoff

int	ndx
int	i, j
int	j1, j2, i1, i2
char	jid

int	fscan(), nscan()

begin
	call eprintf ("Calling calg_init\n")

# Count the entries
	i = fscan (fda)
	call gargd (zptx)
	call gargd (zpty)
	call gargd (xoff)
	call gargd (yoff)
	if (nscan() < 4)
		call fatal (0, "problem reading input file!!")
call eprintf ("%f %f %f %f\n")
call pargd (zptx)
call pargd (zpty)
call pargd (xoff)
call pargd (yoff)
	i1 = (-500.0 - zptx) / xoff
	i2 = (500.0 - zptx) / xoff
	j1 = (180.0 - zpty) / yoff
	j2 = (480.0 - zpty) / yoff
call eprintf ("%d %d %d %d\n")
call pargi (i1)
call pargi (i2)
call pargi (j1)
call pargi (j2)
	ntarg = (abs (j2 - j1) + 1) * (abs (i2 - i1) + 1)
## WORK ABOVE!

# Allocate the vectors
	call malloc (tdat, TDATLEN, TY_STRUCT)
	call malloc (PTINDEX(tdat), ntarg, TY_INT)
	call malloc (PTRA(tdat), ntarg, TY_DOUBLE)
	call malloc (PTDEC(tdat), ntarg, TY_DOUBLE)
	call malloc (PTPA(tdat), ntarg, TY_REAL)
	call malloc (PTLEN1(tdat), ntarg, TY_REAL)
	call malloc (PTLEN2(tdat), ntarg, TY_REAL)
	call malloc (PTWID(tdat), ntarg, TY_REAL)
	call malloc (PTPCODE(tdat), ntarg, TY_INT)
	call malloc (PTSAMPL(tdat), ntarg, TY_INT)

	call malloc (PTSTAT(tdat), ntarg, TY_INT)
	call malloc (PTSEL(tdat), ntarg, TY_INT)
	call malloc (PTXARCS(tdat), ntarg, TY_REAL)
	call malloc (PTYARCS(tdat), ntarg, TY_REAL)
	call malloc (PTRELPA(tdat), ntarg, TY_REAL)

# Allocate slit-index -- cross-reference -- NOTE the zeroing of this array
	call calloc (PTSLNDX(tdat), ntarg, TY_INT)

# Allocate the string array
	call malloc (PTLINE(tdat), ntarg*SZ_LINE, TY_CHAR)


# Generate the trace data
	ndx = 0
	do i = i1, i2 {
		call sprintf (idstr, SZ_ID, "trace%04.0f")
			call pargd (i*xoff)
		alpha = 0.
		delta = 0.
		INDEX(tdat,ndx) = ndx
		RA(tdat,ndx) = DEGTORAD(alpha*15.)
		DEC(tdat,ndx) = DEGTORAD(delta)
		call strcpy (idstr, DATLINE(tdat,ndx), SZ_LINE-1)
		PA(tdat,ndx) = INDEF
		PCODE(tdat,ndx) = 100
		LEN1(tdat,ndx) = 0.5	# DEF_HLEN(indat)
		LEN2(tdat,ndx) = 0.5	# DEF_HLEN(indat)
		SLWID(tdat,ndx) = 2.0	# DEF_SLWID(indat)
		SAMPL(tdat,ndx) = PRIMARY
		SEL(tdat,ndx) = 1
		RELPA(tdat,ndx) = INDEF
		XARCS(tdat,ndx) = i * xoff + zptx
		YARCS(tdat,ndx) = zpty
		ndx = ndx + 1
	}
# Now add the wavelength boxes (spacing 3.1, 4.1, 5.1, 6.1, 7.1 arcmin)
	yoff = 60.
	do i = -7, 7, 2 {
	    jid = 'a'
	    zpty = 192.0
	    do j = 0, 4 {
		call sprintf (idstr, SZ_ID, "box_%s%04.0f")
			call pargc (jid)
			call pargr (i*xoff)
		alpha = 0.
		delta = 0.
		INDEX(tdat,ndx) = ndx
		RA(tdat,ndx) = DEGTORAD(alpha*15.)
		DEC(tdat,ndx) = DEGTORAD(delta)
		call strcpy (idstr, DATLINE(tdat,ndx), SZ_LINE-1)
		PA(tdat,ndx) = INDEF
		PCODE(tdat,ndx) = 100
		LEN1(tdat,ndx) = 0.5	# DEF_HLEN(indat)
		LEN2(tdat,ndx) = 0.5	# DEF_HLEN(indat)
		SLWID(tdat,ndx) = 1.0	# DEF_SLWID(indat)
		SAMPL(tdat,ndx) = PRIMARY
		SEL(tdat,ndx) = 1
		RELPA(tdat,ndx) = INDEF
		XARCS(tdat,ndx) = i * yoff + zptx + i/abs(i)*(j*3.75+xoff/2.)
		YARCS(tdat,ndx) = zpty + j * 60.
		ndx = ndx + 1
		jid = jid + 1
		
	    }
	}
	ntarg = ndx

end


procedure	write_surfcam (sdat, nslit, fdc, maskid)

pointer	sdat
int	nslit
pointer	fdc
char	maskid[ARB]			# mask ID name

int	i

begin
	call fprintf (fdc, "# (     Mask_ID:  %-28s)\n")
		call pargstr (maskid)
	call fprintf (fdc, "# (     X-origin: 14.778 inch from mask edge   )\n")
	call fprintf (fdc, "# (     Y-origin: 0.0    inch from mask edge   )\n")
	call fprintf (fdc, "# (     Max Tool Sz: X.XXX inch                )\n")
	call fprintf (fdc, "# (     Units are:   MM                        )\n")
	call fprintf (fdc, "# (                                            )\n")
	do i = 0, nslit-1 {
		call fprintf (fdc, "%10.4f %10.4f 0.00000\n")
			call pargd (XMM1(sdat,i))
			call pargd (YMM1(sdat,i))
		call fprintf (fdc, "%10.4f %10.4f 0.00000\n")
			call pargd (XMM2(sdat,i))
			call pargd (YMM2(sdat,i))
		call fprintf (fdc, "%10.4f %10.4f 0.00000\n")
			call pargd (XMM3(sdat,i))
			call pargd (YMM3(sdat,i))
		call fprintf (fdc, "%10.4f %10.4f 0.00000\n")
			call pargd (XMM4(sdat,i))
			call pargd (YMM4(sdat,i))
		call fprintf (fdc, "%10.4f %10.4f 0.00000\n")
			call pargd (XMM1(sdat,i))
			call pargd (YMM1(sdat,i))
		call fprintf (fdc, "newrow\n")
	}
end


procedure	write_surfcam_holes (sdat, nslit, fdc, maskid)

pointer	sdat
int	nslit
pointer	fdc
char	maskid[ARB]			# mask ID name

int	i
double	xhole, yhole

begin
	call fprintf (fdc, "# (     Mask_ID:  %-28s)\n")
		call pargstr (maskid)
	call fprintf (fdc, "# (     X-origin: 14.778 inch from mask edge   )\n")
	call fprintf (fdc, "# (     Y-origin: 0.0    inch from mask edge   )\n")
	call fprintf (fdc, "# (     Max Tool Sz: X.XXX inch DRILL HOLES    )\n")
	call fprintf (fdc, "# (     Units are:   MM                        )\n")
	call fprintf (fdc, "# (                                            )\n")
	do i = 0, nslit-1 {
		xhole = 0.25 * (XMM1(sdat,i) + XMM2(sdat,i) +
						XMM3(sdat,i) + XMM4(sdat,i))
		yhole = 0.25 * (YMM1(sdat,i) + YMM2(sdat,i) +
						YMM3(sdat,i) + YMM4(sdat,i))
		call fprintf (fdc, "%10.4f %10.4f 0.00000\n")
			call pargd (xhole)
			call pargd (yhole)
		call fprintf (fdc, "newrow\n")
	}
end


procedure	write_mongo (sdat, nslit, fdp, maskid, plotfile)

pointer	sdat
int	nslit
pointer	fdp
char	maskid[ARB]			# mask ID name
char	plotfile[ARB]			# plotfile name

int	i
int	pen
pointer	mdat

char	rootname[SZ_FNAME]	# attempt to contruct meaningful filename

int	strldx()
long	clktime()
begin
	i = strldx (".", plotfile)
	if (i == 0)
		i = SZ_FNAME
	call strcpy (plotfile, rootname, i-1)
	
	call fprintf (fdp, "psland %-s.ps\n")
		call pargstr (rootname)
	call fprintf (fdp, "color 1 0 0 0 \n")
	call fprintf (fdp, "limits -375.4 375.4 0. 260.\n")
	call fprintf (fdp, "square -1 -1 -1 -1\n")
	do i = 0, nslit-1 {
		call fprintf (fdp, "relocate %8.4f %8.4f\n")
			call pargd (XMM1(sdat,i))
			call pargd (YMM1(sdat,i))
		call fprintf (fdp, "draw   %8.4f %8.4f\n")
			call pargd (XMM2(sdat,i))
			call pargd (YMM2(sdat,i))
		call fprintf (fdp, "draw   %8.4f %8.4f\n")
			call pargd (XMM3(sdat,i))
			call pargd (YMM3(sdat,i))
		call fprintf (fdp, "draw   %8.4f %8.4f\n")
			call pargd (XMM4(sdat,i))
			call pargd (YMM4(sdat,i))
		call fprintf (fdp, "draw   %8.4f %8.4f\n")
			call pargd (XMM1(sdat,i))
			call pargd (YMM1(sdat,i))
	}

	call get_outlines (mdat)
	call fprintf (fdp, "color 2 0 1 0 \n")
	do i = 0, NFPDES(mdat)-1 {
		pen = FPLZ(mdat,i)
		if (pen == 0) {
			call fprintf (fdp, "relocate %8.4f %8.4f\n")
				call pargr (FPLX(mdat,i)*0.7277)
				call pargr (FPLY(mdat,i)*0.7277-ZPT_YM)
			next
		}
		if (pen > 0) {
			if (pen == 1)
				call fprintf (fdp, "ltype 0\n")
			else if (pen == 2)
				call fprintf (fdp, "ltype 2\n")
			else if (pen == 3)
				call fprintf (fdp, "ltype 1\n")
		}
		call fprintf (fdp, "draw   %8.4f %8.4f\n")
			call pargr (FPLX(mdat,i)*0.7277)
			call pargr (FPLY(mdat,i)*0.7277-ZPT_YM)
	}

	call fprintf (fdp, "limits 0 1 0 1 \n")
	call fprintf (fdp, "color 1 0 0 0 \n")

	call fprintf (fdp, "expand 1.3 \n")
	call fprintf (fdp, "rel 0.02 1.25 \n")
	call fprintf (fdp, "putlabel 6 ID: %s \n")
		call pargstr (maskid)

	call fprintf (fdp, "expand 1.1 \n")
	call fprintf (fdp, "rel 0.02 1.16 \n")
	call fprintf (fdp, "putlabel 6 (plotfile: %s) \n")
		call pargstr (plotfile)

	call fprintf (fdp, "expand 0.7 \n")
	call fprintf (fdp, "rel 0.02 1.08 \n")
	call cnvtime (clktime (0), rootname, SZ_FNAME)
	call fprintf (fdp, "putlabel 6 (%s) \n")
		call pargstr (rootname)

	call fprintf (fdp, "hard \n")

end


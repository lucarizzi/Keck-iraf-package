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

pointer	fda, fdc				# in/out file descriptors
pointer	fdp			# pointer to optional plot file

pointer	indat			# pointer to instr/tel data structure

int	ntarg
pointer	tdat			# pointer to data structure for targets

int	nslit
pointer	sdat			# pointer to data structure for slits

int	i	## TMP
int	pen	## TMP
real	offset, stdwid	## TMP?
pointer	mdat	## TMP?

bool	clgetb(), strne()
pointer	open()

begin
	call clgstr ("objfile", objfile, SZ_FNAME)
#	call clgstr ("output", output, SZ_FNAME)
	call clgstr ("mdf", mdf, SZ_FNAME)
	call clgstr ("plotfile", plotfile, SZ_FNAME)

# Read in telescope, default data:
	call data_init (indat)

# Open the input list of targets:
        fda = open (objfile, READ_ONLY, TEXT_FILE)

        fdc = open (mdf, NEW_FILE, TEXT_FILE)

# Read in the target data:
# Target is a grid!  Just fill the arrays
#	call targ_init (fda, tdat, ntarg, indat)

# Assign CENTER to telescope axis
	RA_TEL(indat) =  RA_FLD(INDAT)
	DEC_TEL(indat) = DEC_FLD(INDAT)
	PA_ROT(indat) = 0.

# Calculate position wrt telescope axis:
	call tel_coords (tdat, ntarg, indat)

# Turn on ALL stat flags (special case)
	call amovki (YES, STAT(tdat,0), ntarg)

# Ready to display and be interactive
	call select_targ (tdat, ntarg, indat)

# Save output?

# Calculate the SLITS in celestial coordinates.
# It is here that slit lengths, edges must be resolved.  Does binding go here?
# Or is this a step that should actually be included above?
#

	call gen_slits (tdat, ntarg, sdat, nslit, indat)
#	call len_slits (tdat, ntarg, sdat, nslit, indat)

call eprintf (" calling sky_coord \n")
# Convert adopted slits into coordinated on the sky:
	call sky_coords (sdat, nslit, indat)

# ADD SPECIAL FOR COHU MASK:
	do i = 9, 14 {
		SLWID(sdat,i) = SLWID(sdat,i) + 0.5
	}
	do i = 15, 20 {
		SLWID(sdat,i) = SLWID(sdat,i) + 1.0
	}
	do i = 21, 26 {
		SLWID(sdat,i) = SLWID(sdat,i) + 1.5
	}




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
	stdwid = YMM3(sdat,1) - YMM2(sdat,1)
	do i = 9, 24, 3 {
		offset = 0.5 * (YMM3(sdat,i) - YMM2(sdat,i) - stdwid)
		YMM1(sdat,i) = YMM1(sdat,i) - offset
		YMM2(sdat,i) = YMM2(sdat,i) - offset
		YMM3(sdat,i) = YMM3(sdat,i) - offset
		YMM4(sdat,i) = YMM4(sdat,i) - offset
	}
	do i = 11, 26, 3 {
		offset = 0.5 * (YMM3(sdat,i) - YMM2(sdat,i) - stdwid)
		YMM1(sdat,i) = YMM1(sdat,i) + offset
		YMM2(sdat,i) = YMM2(sdat,i) + offset
		YMM3(sdat,i) = YMM3(sdat,i) + offset
		YMM4(sdat,i) = YMM4(sdat,i) + offset
	}


 
#
call eprintf (" calling write_design ...\n")
	call write_design (indat, tdat, ntarg, sdat, nslit, fdc, objfile)

	if (strne (plotfile, "")) {
            fdp = open (plotfile, NEW_FILE, TEXT_FILE)
	    call fprintf (fdp, "psland mask.ps\n")
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
	    do i = 0, NFPDES(mdat)-1 {
		pen = FPLZ(mdat,i)
		if (pen == 0) {
			call fprintf (fdp, "relocate %8.4f %8.4f\n")
				call pargr (FPLX(mdat,i)*0.722)
				call pargr (FPLY(mdat,i)*0.722-ZPT_YM)
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
			call pargr (FPLX(mdat,i)*0.722)
			call pargr (FPLY(mdat,i)*0.722-ZPT_YM)
	    }
	    call fprintf (fdp, "hard\n")
	    call close (fdp)
	}

i = -1
# i = 0
if (i==0){
# SPECIAL FOR MILLING:
	call fprintf (fdc, "# (     Mask_ID:  COHU_0                       )\n")
	call fprintf (fdc, "# (     X-origin: 14.778 inch from mask edge   )\n")
	call fprintf (fdc, "# (     Y-origin: 0.0    inch from mask edge   )\n")
	call fprintf (fdc, "# (     Max Tool Sz: 0.010 inch                )\n")
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
}

	call close (fdc)



#
call eprintf (" [recall: still must check proper x,y and tan mapping]\n")
#
end



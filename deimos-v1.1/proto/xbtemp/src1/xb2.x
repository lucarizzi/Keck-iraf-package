# Check existence of file name.
# todo:  general routine to read in table names EXTNAME

######################### CFITSIO routines ##########################
## NB: local "fitsio.h" comes from all the defines in cfitsio.h (?) with
## appropriate characters replaced.  XXX NOTE IN MAINTENANCE DOC

include	<math.h>
include "instrument.h"
include "align.h"
include "fitsio.h"


procedure	t_get_boxes ()

int	pos			# Slot number

char	outname[SZ_FNAME]
char	guiname[SZ_GUINAM]
int	i
int	nbox
pointer	bdat
pointer	fd			# output file descriptor

bool	streq()
int	clgeti()
pointer	open()
begin
	pos = clgeti ("SLMSKPOS")		# eventually we read SLMSKPOS

	call get_boxes (pos, bdat, nbox, guiname)

	if (nbox == 0)
		call fatal (0, "No Alignment Boxes Found!")

# Construct output file name; open
	if (streq (guiname, "")) {
		call sprintf (outname, SZ_FNAME, "box.%02d")
			call pargi (pos)
	} else {
		call sprintf (outname, SZ_FNAME, "box.%s")
			call pargstr (guiname)
	}
	fd = open (outname, NEW_FILE, TEXT_FILE)
	call eprintf ("Writing file %s \n")
		call pargstr (outname)

	do i = 0, nbox-1 {
		if (XPX(bdat,i) == INDEF || YPX(bdat,i) == INDEF) {
		    call fprintf (fd,
			"# OFF_CHIP %s   # Mag_%s=%5.2f %8.3f %7.3f\n")
			call pargstr (OBJNAME(bdat,i))
			call pargc (PBAND(bdat,i))
			call pargr (MAG(bdat,i))
			call pargr (XMM(bdat,i))
			call pargr (YMM(bdat,i))
		} else {
		    call fprintf (fd,
			"%6.1f %6.1f  %s   # Mag_%s=%5.2f %8.3f %7.3f\n")
			call pargr (XPX(bdat,i))
			call pargr (YPX(bdat,i))
			call pargstr (OBJNAME(bdat,i))
			call pargc (PBAND(bdat,i))
			call pargr (MAG(bdat,i))
			call pargr (XMM(bdat,i))
			call pargr (YMM(bdat,i))
		}
	}

	call close (fd)

end


procedure	get_boxes (pos, bdat, nbox, guiname)

int	pos			# Position from SLMSKPOS
pointer	bdat
int	nbox
char	guiname[SZ_GUINAM]

char	fname[SZ_FNAME]

int	i

# needed for qxfm
char	amap[SZ_FNAME], bmap[SZ_FNAME]		# input mappings

pointer	fda, fdb
pointer	map[4]				# pointers to surf fits (1,2=amap;3,4=b)

double	sys[NPARAM]			# system parameters
real	ccd[NCCD,3]			# CCD geometry
double	a3[3,3]				# grating transform
real	xics, yics			# pixel values in ICS
real	scaling
int	stat
int	chip

int	qxfm(), qxfm_init()

real	clgetr()
pointer	open()

begin
	if (pos < 2 || pos > 12)
		call fatal (0, "Illegal pos value")

	call sprintf (fname, SZ_FNAME, "%sSlot.%d.fits")
		call pargstr (DYNA_DIR)
		call pargi (pos)

call eprintf ("Looking for %s\n")
call pargstr (fname)

	call get_all_info (fname, bdat, nbox, guiname)

# Now get qmodel values:
	MU(sys) = DEGTORAD (clgetr ("qmodel.mu"))
	GR_YERR(sys) = DEGTORAD (clgetr ("qmodel.roll3"))
	GR_ZERR(sys) = DEGTORAD (clgetr ("qmodel.o3"))

	scaling = 1. + clgetr ("qmodel.scale_adj")

	ORDER(sys) = 0		# always zeroth order
	GRLINES(sys) = 1	# always zeroth order

	call clgstr ("qmodel.amap", amap, SZ_FNAME)
	fda = open (amap, READ_ONLY, TEXT_FILE)

	call clgstr ("qmodel.bmap", bmap, SZ_FNAME)
	fdb = open (bmap, READ_ONLY, TEXT_FILE)

# Initialize the mappings
	stat = qxfm_init (fda, fdb, map, a3, sys, ccd)


	do i = 0, nbox-1 {

# calculate the mapping
		stat = qxfm (map, a3, sys, ccd, XMM(bdat,i), YMM(bdat,i), 1.,
			scaling, xics, yics, XPX(bdat,i), YPX(bdat,i), chip)

		if (stat == OFF_CHIP) {
			XPX(bdat,i) = INDEF
			YPX(bdat,i) = INDEF
		}
	}
end



procedure	get_all_info (fname, bdat, nbox, guiname)

char	fname[ARB]		# FITS file name w/o extn
pointer	bdat			# box data struct
int	nbox			# number of boxes
char	guiname[ARB]		# GI Name for slitmask

char	xname[SZ_FNAME]		# FITS name with extention
char	xline[80]		# XXX hardcode
%	character*80 f77nam
%	character*80 f77lin
%	character*80 f77nul

int	lu
int	stat
int	nc, nr

int	i, j, k
int	n, ndx
int	nslit
int	nobj
real	x, y, xsum, ysum

bool	isnul			# Null value(s) present?
int	inul			# dummy null value
real	rnul			# dummy null value

bool	strne()
begin
	rnul = INDEF
	call f77pak ("???", f77nul, 80)


# Get DesiSlits table (3) for slitTyp=A (col 5)
# call eprintf ("DEBUG: table[3]:\n")

	call strcpy (fname, xname, SZ_FNAME)
	call strcat ("[3]", xname, SZ_FNAME)
	call f77pak (xname, f77nam, SZ_FNAME)

	call ftgiou (lu, stat)
	call ftnopn (lu, f77nam, READONLY, stat)

	call ftgncl (lu, nc, stat)
	call ftgnrw (lu, nr, stat)

# check format:
	if (nc != 10) {
		call eprintf ("%s, nc=%d nr=%d stat=%d\n")
		call pargstr (xname)
		call pargi (nc)
		call pargi (nr)
		call pargi (stat)
		call fatal (0, "Unexpected format in input FITS table")
	}

	nslit = nr
	call box_alloc (bdat, nslit)

	ndx = 0
	do i = 0, nslit-1 {
		j = i + 1
		call ftgcvs (lu, 5, j, 1, 1, f77nul, f77lin, isnul, stat)
		call f77upk (f77lin, xline, 1)
		if (strne (xline, "A"))
			next

# we have an alignment box ... get dSlitId (col 1)
		call ftgcvj (lu, 1, j, 1, 1, inul, SLITNO(bdat,ndx), isnul, stat)
# XXX and box size (col 6,8)?
# increment counter
		ndx = ndx + 1
	}
	nbox = ndx
	call ftclos (lu, stat)



# Get BluSlits table (6);earch dSlitId (Col 3), get average of xmm, ymm (col 4+)
# call eprintf ("DEBUG: table[6]:\n")

	call strcpy (fname, xname, SZ_FNAME)
	call strcat ("[6]", xname, SZ_FNAME)
	call f77pak (xname, f77nam, SZ_FNAME)

	call ftgiou (lu, stat)
	call ftnopn (lu, f77nam, READONLY, stat)

	call ftgncl (lu, nc, stat)
	call ftgnrw (lu, nr, stat)

# check format:
	if (nc != 11 || nr != nslit) {
		call eprintf ("%s, nc=%d nr=%d stat=%d\n")
		call pargstr (xname)
		call pargi (nc)
		call pargi (nr)
		call pargi (stat)
		call fatal (0, "Unexpected format in input FITS table")
	}

	do i = 0, nbox-1 {
		ndx = SLITNO(bdat,i)
		do k = 0, nslit-1 {
			j = k + 1
			call ftgcvj (lu, 3, j, 1, 1, inul, n, isnul, stat)
			if (n != ndx)
				next

			xsum = 0.
			ysum = 0.

			call ftgcve (lu, 4, j, 1, 1, rnul, x, isnul, stat)
			xsum = xsum + x
			call ftgcve (lu, 5, j, 1, 1, rnul, y, isnul, stat)
			ysum = ysum + y

			call ftgcve (lu, 6, j, 1, 1, rnul, x, isnul, stat)
			xsum = xsum + x
			call ftgcve (lu, 7, j, 1, 1, rnul, y, isnul, stat)
			ysum = ysum + y

			call ftgcve (lu, 8, j, 1, 1, rnul, x, isnul, stat)
			xsum = xsum + x
			call ftgcve (lu, 9, j, 1, 1, rnul, y, isnul, stat)
			ysum = ysum + y

			call ftgcve (lu, 10, j, 1, 1, rnul, x, isnul, stat)
			xsum = xsum + x
			call ftgcve (lu, 11, j, 1, 1, rnul, y, isnul, stat)
			ysum = ysum + y

			XMM(bdat,i) = xsum / 4.
			YMM(bdat,i) = ysum / 4.
		}
	}
	call ftclos (lu, stat)
			


# Get SlitObjMap table (4); search dSlitId (Col 3), get ObjectId (col 2)
# call eprintf ("DEBUG: table[4]:\n")

	call strcpy (fname, xname, SZ_FNAME)
	call strcat ("[4]", xname, SZ_FNAME)
	call f77pak (xname, f77nam, SZ_FNAME)

	call ftgiou (lu, stat)
	call ftnopn (lu, f77nam, READONLY, stat)

	call ftgncl (lu, nc, stat)
	call ftgnrw (lu, nr, stat)

# check format:
	if (nc != 5 || nr < nslit) {
		call eprintf ("%s, nc=%d nr=%d stat=%d\n")
		call pargstr (xname)
		call pargi (nc)
		call pargi (nr)
		call pargi (stat)
		call fatal (0, "Unexpected format in input FITS table")
	}

	nobj = nr
	do i = 0, nbox-1 {
		ndx = SLITNO(bdat,i)
		do k = 0, nobj-1 {
			j = k + 1
			call ftgcvj (lu, 3, j, 1, 1, inul, n, isnul, stat)
			if (n != ndx)
				next

			call ftgcvj (lu, 2, j, 1, 1, inul, OBJNO(bdat,i), isnul, stat)
			break
		}
	}
	call ftclos (lu, stat)
			


# Get ObjectCat table (1); search ObjectId (Col 1), get OBJECT (col 2), etc
# call eprintf ("DEBUG: table[1]:\n")

	call strcpy (fname, xname, SZ_FNAME)
	call strcat ("[1]", xname, SZ_FNAME)
	call f77pak (xname, f77nam, SZ_FNAME)

	call ftgiou (lu, stat)
	call ftnopn (lu, f77nam, READONLY, stat)

	call ftgncl (lu, nc, stat)
	call ftgnrw (lu, nr, stat)

# check format:
	if (nc != 12 || nr < nslit) {
		call eprintf ("%s, nc=%d nr=%d stat=%d\n")
		call pargstr (xname)
		call pargi (nc)
		call pargi (nr)
		call pargi (stat)
		call fatal (0, "Unexpected format in input FITS table")
	}

	do i = 0, nbox-1 {
		ndx = OBJNO(bdat,i)
		do k = 0, nobj-1 {
			j = k + 1
			call ftgcvj (lu, 1, j, 1, 1, inul, n, isnul, stat)
			if (n != ndx)
				next

			call ftgcvs (lu, 2, j, 1, 1, f77nul, f77lin, isnul, stat)
			call f77upk (f77lin, OBJNAME(bdat,i), SZ_ID-1)
			call ftgcve (lu, 8, j, 1, 1, -99., MAG(bdat,i), isnul, stat)
			call ftgcvs (lu, 9, j, 1, 1, f77nul, f77lin, isnul, stat)
			call f77upk (f77lin, xline, 1)
			PBAND(bdat,i) = xline[1]
			if (xline[2] != EOS)
				call eprintf ("Warning; pband truncated\n")
			break
		}
	}
	call ftclos (lu, stat)


# Get MaskBlu table (5) for guiname (col 17)
# call eprintf ("DEBUG: table[5]:\n")

	call strcpy (fname, xname, SZ_FNAME)
	call strcat ("[5]", xname, SZ_FNAME)
	call f77pak (xname, f77nam, SZ_FNAME)

	call ftgiou (lu, stat)
	call ftnopn (lu, f77nam, READONLY, stat)

	call ftgncl (lu, nc, stat)
	call ftgnrw (lu, nr, stat)

# check format:
	if (nc != 18 || nr!= 1) {
		call eprintf ("%s, nc=%d nr=%d stat=%d\n")
		call pargstr (xname)
		call pargi (nc)
		call pargi (nr)
		call pargi (stat)
		call fatal (0, "Unexpected format in input FITS table")
	}

	call ftgcvs (lu, 17, 1, 1, 1, f77nul, f77lin, isnul, stat)
	if (isnul)
		call strcpy ("", guiname, SZ_GUINAM)
	else
		call f77upk (f77lin, guiname, SZ_GUINAM)

call eprintf ("DEBUG: guiname= %s\n")
call pargstr (guiname)

	call ftclos (lu, stat)


end

#
# BOX_ALLOC: allocate arrays for Boxes (stolen from targ_init)
#
procedure	box_alloc (bdat, nbox)

pointer	bdat
int	nbox

begin
# Allocate the vectors
	call malloc (bdat, BDATLEN, TY_STRUCT)
	call malloc (PTSLTNO(bdat), nbox, TY_INT)
	call malloc (PTOBJNO(bdat), nbox, TY_INT)
	call malloc (PTMAG(bdat), nbox, TY_REAL)
	call malloc (PTPBAND(bdat), nbox, TY_CHAR)	## XXX should be general
	call malloc (PTXMM(bdat), nbox, TY_REAL)
	call malloc (PTYMM(bdat), nbox, TY_REAL)
	call malloc (PTXPX(bdat), nbox, TY_REAL)
	call malloc (PTYPX(bdat), nbox, TY_REAL)

# Allocate the string name array
	call malloc (PTNAME(bdat), nbox*SZ_ID, TY_CHAR)

end


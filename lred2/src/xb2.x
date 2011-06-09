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

char	fname[SZ_FNAME]
char	outname[SZ_FNAME]
char	guiname[SZ_GUINAM]
int	i
int	nbox
pointer	bdat
pointer	fd			# output file descriptor

real	x0ics, y0ics

bool	streq()
int	clgeti()
pointer	open()
begin
	call clgstr ("image", fname, SZ_FNAME)

	if (streq ("", fname)) {
		pos = clgeti ("SLMSKPOS")	# eventually we read SLMSKPOS
		if (pos < 2 || pos > 12) {
			call fatal (0, "Illegal pos value")
		} else {
# build file name
			call sprintf (fname, SZ_FNAME, "%sSlot.%d.fits")
				call pargstr (DYNA_DIR)
				call pargi (pos)
		}
	}


call eprintf ("Looking for %s\n")
call pargstr (fname)

# Extract the info:
	call get_boxes (fname, bdat, nbox, guiname, x0ics, y0ics)

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
			"# OFF_CHIP   %s   %1s %5.2f # %8.3f %7.3f\n")
			call pargstr (OBJNAME(bdat,i))
			call pargc (PBAND(bdat,i))
			call pargr (MAG(bdat,i))
			call pargr (XMM(bdat,i))
			call pargr (YMM(bdat,i))
		} else {
		    call fprintf (fd,
			"%6.1f %6.1f  %s   %1s %5.2f # %8.3f %7.3f\n")
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


#
# GET_BOXES: figure out locations of alignment boxes
#

procedure	get_boxes (fname, bdat, nbox, guiname, x0ics, y0ics)

char	fname[ARB]
pointer	bdat
int	nbox
char	guiname[SZ_GUINAM]
real	x0ics, y0ics			# ICS coords of PO  XXX


int	i

# needed for qxfm
char	amap[SZ_FNAME], bmap[SZ_FNAME]		# input mappings

pointer	fda, fdb
pointer	map[4]				# pointers to surf fits (1,2=amap;3,4=b)

double	sys[NPARAM]			# system parameters
real	ccd[NCCD,3]			# CCD geometry
double	a3[3,3]				# grating transform
real	xics, yics			# pixel values in ICS
real	xrot, yrot
real	scaling
int	stat
int	chip

int	fits_ize()
int	qxfm(), qxfm_init()

real	clgetr()
pointer	open()

begin

call eprintf ("Looking for %s\n")
call pargstr (fname)

# resolve FITS file name:
	if (fits_ize (fname, SZ_FNAME) != OK)
		call fatal (0, "No FITS file")

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
		scaling, xics, yics, XPX(bdat,i), YPX(bdat,i), chip, YES, YES)

		if (stat != ON_CHIP) {
			XPX(bdat,i) = INDEF
			YPX(bdat,i) = INDEF
		}
	}

# get ICS of PO:    XXX -- HARCODED
	stat = qxfm (map, a3, sys, ccd, 0., 68., 1.,
			scaling, x0ics, y0ics, xrot, yrot, chip)
end

#
# GET_ALL_INFO: gets info on alignment stars from Mask Design FITS tables.
# Generic (I think);  the only hardcode is that the XMM,YMM must be in order
# following XMM1.

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
int	nr

int	i, j, k
int	n, ndx
int	nslit
int	nobj
int	ic1, ic2, ic3, ic4
real	x, y, xsum, ysum

bool	isnul			# Null value(s) present?
int	inul			# dummy null value
real	rnul			# dummy null value

int	posn_hdu(), get_col_ndx()
bool	strne()
begin
	rnul = INDEF
	call f77pak ("???", f77nul, 80)

	stat = 0		# apparently needed for FITSIO

	call f77pak (fname, f77nam, SZ_FNAME)

	call ftgiou (lu, stat)
	call ftnopn (lu, f77nam, READONLY, stat)

# Get DesiSlits table (3) for slitTyp=A (col 5)

	if (posn_hdu (lu, "DesiSlits") != OK)
		call fatal (0, "")

	if (get_col_ndx (lu, "slitTyp", ic1) != OK)
		call fatal (0, "column must exist!")
	if (get_col_ndx (lu, "dSlitId", ic2) != OK)
		call fatal (0, "column must exist!")

	call ftgnrw (lu, nr, stat)
	nslit = nr

	call box_alloc (bdat, nslit)

	ndx = 0
	do i = 0, nslit-1 {
		j = i + 1
		call ftgcvs (lu, ic1, j, 1, 1, f77nul, f77lin, isnul, stat)
		call f77upk (f77lin, xline, 1)
		if (strne (xline, "A"))
			next

# we have an alignment box ... get dSlitId
		call ftgcvj (lu, ic2, j, 1, 1, inul, SLITNO(bdat,ndx), isnul, stat)
# XXX ... and box size (col 6,8)?

# increment counter
		ndx = ndx + 1
	}
	nbox = ndx


# Get BluSlits table; search dSlitId, get average of xmm, ymm
# NB: the order in the table is important here; partially checked

	if (posn_hdu (lu, "BluSlits") != OK)
		call fatal (0, "")

	if (get_col_ndx (lu, "dSlitId", ic1) != OK)
		call fatal (0, "column must exist!")
	if (get_col_ndx (lu, "slitX1", ic2) != OK)
		call fatal (0, "column must exist!")
	if (get_col_ndx (lu, "slitY4", ic4) != OK)
		call fatal (0, "column must exist!")

# check format:
	if (ic4 != ic2+7 || nr != nslit) {
		call eprintf ("%s, nr=%d stat=%d\n")
		call pargstr (xname)
		call pargi (nr)
		call pargi (stat)
		call fatal (0, "Unexpected format in BluSlits table")
	}

	do i = 0, nbox-1 {
		ndx = SLITNO(bdat,i)
		do k = 0, nslit-1 {
			j = k + 1
			call ftgcvj (lu, ic1, j, 1, 1, inul, n, isnul, stat)
			if (n != ndx)
				next

			xsum = 0.
			ysum = 0.

			call ftgcve (lu, ic2,   j, 1, 1, rnul, x, isnul, stat)
			xsum = xsum + x
			call ftgcve (lu, ic2+1, j, 1, 1, rnul, y, isnul, stat)
			ysum = ysum + y

			call ftgcve (lu, ic2+2, j, 1, 1, rnul, x, isnul, stat)
			xsum = xsum + x
			call ftgcve (lu, ic2+3, j, 1, 1, rnul, y, isnul, stat)
			ysum = ysum + y

			call ftgcve (lu, ic2+4, j, 1, 1, rnul, x, isnul, stat)
			xsum = xsum + x
			call ftgcve (lu, ic2+5, j, 1, 1, rnul, y, isnul, stat)
			ysum = ysum + y

			call ftgcve (lu, ic2+6, j, 1, 1, rnul, x, isnul, stat)
			xsum = xsum + x
			call ftgcve (lu, ic2+7, j, 1, 1, rnul, y, isnul, stat)
			ysum = ysum + y

			XMM(bdat,i) = xsum / 4.
			YMM(bdat,i) = ysum / 4.
		}
	}
			


# Get SlitObjMap table; search on dSlitId to get ObjectId

	if (posn_hdu (lu, "SlitObjMap") != OK)
		call fatal (0, "")

	if (get_col_ndx (lu, "dSlitId", ic1) != OK)
		call fatal (0, "column must exist!")
	if (get_col_ndx (lu, "ObjectId", ic2) != OK)
		call fatal (0, "column must exist!")

	call ftgnrw (lu, nr, stat)
	nobj = nr
	do i = 0, nbox-1 {
		ndx = SLITNO(bdat,i)
		do k = 0, nobj-1 {
			j = k + 1
			call ftgcvj (lu, ic1, j, 1, 1, inul, n, isnul, stat)
			if (n != ndx)
				next

			call ftgcvj (lu, ic2, j, 1, 1, inul, OBJNO(bdat,i), isnul, stat)
			break
		}
	}



# Get ObjectCat table (1); search ObjectId (Col 1), get OBJECT (col 2), etc

	if (posn_hdu (lu, "ObjectCat") != OK)
		call fatal (0, "")

### WORK HERE
	if (get_col_ndx (lu, "ObjectId", ic1) != OK)
		call fatal (0, "column must exist!")
	if (get_col_ndx (lu, "OBJECT", ic2) != OK)
		call fatal (0, "column must exist!")
	if (get_col_ndx (lu, "mag", ic3) != OK)
		call fatal (0, "column must exist!")
	if (get_col_ndx (lu, "pBand", ic4) != OK)
		call fatal (0, "column must exist!")

	do i = 0, nbox-1 {
		ndx = OBJNO(bdat,i)
		do k = 0, nobj-1 {
			j = k + 1
			call ftgcvj (lu, ic1, j, 1, 1, inul, n, isnul, stat)
			if (n != ndx)
				next

			call ftgcvs (lu, ic2, j, 1, 1, f77nul, f77lin, isnul, stat)
			call f77upk (f77lin, OBJNAME(bdat,i), SZ_ID-1)
			call ftgcve (lu, ic3, j, 1, 1, -99., MAG(bdat,i), isnul, stat)
			call ftgcvs (lu, ic4, j, 1, 1, f77nul, f77lin, isnul, stat)
			call f77upk (f77lin, xline, 1)
			PBAND(bdat,i) = xline[1]
			if (xline[2] != EOS)
				call eprintf ("Warning; pband truncated\n")
			break
		}
	}


# Get MaskBlu table for guiname

	if (posn_hdu (lu, "MaskBlu") != OK)
		call fatal (0, "")

	if (get_col_ndx (lu, "guiname", ic1) != OK)
		call fatal (0, "column must exist!")

	call ftgcvs (lu, ic1, 1, 1, 1, f77nul, f77lin, isnul, stat)
	if (isnul)
		call strcpy ("", guiname, SZ_GUINAM)
	else
		call f77upk (f77lin, guiname, SZ_GUINAM)


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

#
# FITS_IZE: append appropriate suffix as needed for FITSIO routines
#

int	procedure fits_ize (name, maxlen)

char	name[ARB]
int	maxlen

char	xname[SZ_FNAME], kw[8]
int	len
pointer	fd

bool	streq()
int	access(), strlen(), fscan()
pointer	open()

begin
	len = strlen (name)
	if (len > maxlen-5)			# XXX could clean up
		call eprintf ("WARNING: file name length close to limit!!\n")

	call strcpy (name, xname, SZ_FNAME)

	
	if (access (xname, 0, 0) == NO) {
		call strcat (".fit", xname, SZ_FNAME)
		if (access (xname, 0, 0) == NO) {
			call strcat ("s", xname, SZ_FNAME)
			if (access (xname, 0, 0) == NO) {
				call eprintf ("WARNING: %s does not exist\n")
					call pargstr (name)
					return (ERR)
			}
		}
	}

	len = strlen (xname)
	if (len > maxlen)
		call fatal (0, "FITS_ized name too long!!")

	call strcpy (xname, name, maxlen)

# Check for SIMPLE keyword as test of FITS file:

	fd = open (name, READ_ONLY, TEXT_FILE)
	len = fscan (fd)
	call gargwrd (kw, 8)
	call close (fd)
	if (streq (kw, "SIMPLE")) {
		return (OK)
	} else {
		call eprintf ("ERROR: %s cannot be a FITS file!\n")
			call pargstr (xname)
		return (ERR)
	}
end


#
# POSN_HDU: position the FITSIO routines to a particular table
#

int	procedure posn_hdu (lu, extname)

int	lu			# unit number from FITSIO
char	extname[ARB]		# table name

int	stat
%	character*80 f77nam

begin
	call f77pak (extname, f77nam, SZ_FNAME)
	call ftmnhd (lu, ANY_HDU, f77nam, 0, stat)

	if (stat == BAD_HDU_NUM) {
		call eprintf ("Extension does not exist (%-s)\n")
			call pargstr (extname)
		return (ERR)
	} else {
		return (OK)
	}
end

#
# GET_COL_NDX: get a column index in a FITS table
#

int	procedure get_col_ndx (lu, name, ndx)

int	lu			# FITSIO unit number of file
char	name[ARB]		# column name
int	ndx			# returned index

int	stat
%	character*80 f77nam
begin
	stat = 0			# clear status flag

	call f77pak (name, f77nam, 80)
	call ftgcno (lu, true, f77nam, ndx, stat)

	if (stat == COL_NOT_FOUND) {
		call eprintf ("WARNING: column %s not found!\n")
			call pargstr (name)
		ndx = -1
		return (ERR)
	} else {
		return (OK)
	}
end

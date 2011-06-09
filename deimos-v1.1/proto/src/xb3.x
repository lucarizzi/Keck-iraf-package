# Check existence of file name.
# todo:  general routine to read in table names EXTNAME

######################### CFITSIO routines ##########################
## NB: local "fitsio.h" comes from all the defines in cfitsio.h (?) with
## appropriate characters replaced.  XXX NOTE IN MAINTENANCE DOC

include	<math.h>
include "instrument.h"
include "align.h"
include "fitsio.h"

procedure	t_itest ()

char	image[SZ_FNAME]
char	tblname[SZ_FNAME]

int	lu, ndx, stat

begin
	call clgstr ("image", image, SZ_FNAME)
	call clgstr ("table", tblname, SZ_FNAME)

	call get_tbl_ndx (image, tblname, YES, lu, ndx)

call eprintf ("%s = %s[%d]\n")
call pargstr (tblname)
call pargstr (image)
call pargi (ndx)

	call get_tbl_ndx (image, "DesiSlits", NO, lu, ndx)
call eprintf ("%s = %s[%d]\n")
call pargstr (tblname)
call pargstr (image)
call pargi (ndx)

	call get_tbl_ndx (image, "SlitObjMap", NO, lu, ndx)
call eprintf ("%s = %s[%d]\n")
call pargstr (tblname)
call pargstr (image)
call pargi (ndx)

	call ftclos (lu, stat)

	call get_tbl_ndx (image, "DesiSlits", YES, lu, ndx)
call eprintf ("%s = %s[%d]\n")
call pargstr (tblname)
call pargstr (image)
call pargi (ndx)

	call ftclos (lu, stat)

	call get_tbl_ndx (image, "SlitObjMap", YES, lu, ndx)
call eprintf ("%s = %s[%d]\n")
call pargstr (tblname)
call pargstr (image)
call pargi (ndx)

	call ftclos (lu, stat)

end

procedure	get_tbl_ndx (fname, desname, open_it, lu, ndx)

char	fname[ARB]		# name of image
char	desname[ARB]		# desired table name
int	open_it			# open the file?
int	lu			# returned (open) unit number for FITSIO
int	ndx			# returned 0-indexed number of HDU

int	stat

%	character*80 f77nam
%	character*80 extnam

int	access()
begin

call eprintf ("DEBUGA %d\n"); call pargi (stat)
	ndx = -1
	stat = 0	# need for FITSIO, which doesn't clean this out

	if (access (fname, 0, 0) == NO) {
		call eprintf ("File does not exist (%-s)\n")
			call pargstr (fname)
		return
	}

call eprintf ("DEBUGB %d\n"); call pargi (stat)
	call f77pak (fname, f77nam, SZ_FNAME)
	call f77pak (desname, extnam, SZ_FNAME)
	
	if (open_it == YES) {
		call ftgiou (lu, stat)
		call ftnopn (lu, f77nam, READONLY, stat)
	}

call eprintf ("DEBUGC %d\n"); call pargi (stat)
	call ftmnhd (lu, ANY_HDU, extnam, 0, stat)

	if (stat == BAD_HDU_NUM) {
		call eprintf ("Extension does not exist (%-s)\n")
			call pargstr (desname)
		return
	}

call eprintf ("DEBUGD\n")
	call ftghdn (lu, ndx)
	ndx = ndx - 1
end


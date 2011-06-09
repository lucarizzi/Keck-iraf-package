#
# WRITE_DESIGN: generate the FITS table design
#

include	<time.h>
include	"ftab.h"

procedure	write_design (indat, tdat, ntarg, sdat, nslit, fd, objfile)

pointer	indat			# pointer to instr/telesc data struct

pointer	tdat			# pointer to targets data struct
int	ntarg

pointer	sdat			# pointer to slits data struct
int	nslit

pointer	fd			# output file descriptor
char	objfile[ARB]		# name of input file

int	i

char	ident[SZ_ID]
int	stat
int	nselect


char	fmtstr[SZ_LINE]				# format for writing table row
int	fmlen
char	date[SZ_TIME]					# time/date of file
pointer	kdat				# FITS table data structure
# int	inul
# real	rnul
# double	dnul
char	snul[1]
#####
char	creatask[60]
char	person[30]
char	unknown[3]

char	acode

char	fmt_mm[4], mm[2], deg[3]	# define to avoid filling declaration space

char	cnul[6]
int	inul
real	rnul
double	dnul

int	sscan()
begin
	call strcpy ("F9.3", fmt_mm, 4)
	call strcpy ("mm", mm, 2)
	call strcpy ("deg", deg, 3)

	fmlen = SZ_LINE

	call strcpy ("INDEF", cnul, 6)
	inul = -9999
	rnul = -9999.0		# Careful! will write it out longer if it
	dnul = -9999.0		# doesn't fit into the field

# SET UP FITS FILE
	call strcpy ("", snul, 1)
	call strcpy ("Phillips <phillips@ucolick.org>", person, 35)
	call strcpy ("DSIMULATOR -- June2000", creatask, 35)
	call strcpy ("???", unknown, 3)

	call ft_date (date, SZ_TIME)

# Write the primary HDU
	call ftab_whdu0 (fd)

call eprintf ("    table 1 ")
# Only include selected objects
	nselect = 0
	do i = 0, ntarg-1 {
		if (SEL(tdat,i) == YES)
			nselect = nselect + 1
	}

########  WRITE ObjectCat TABLE
	call ftab_init (kdat, 18, fmtstr)

	call ftcol_defi (kdat, "ObjectId",  "I6", snul, inul,  6, fmtstr, fmlen)
	call ftcol_defc (kdat, "OBJECT",   "A68", snul, cnul, 68, fmtstr, fmlen)
	call ftcol_defd (kdat, "RA_OBJ", "F12.8", deg, dnul, 12, fmtstr, fmlen)
	call ftcol_defd (kdat, "DEC_OBJ","F12.8", deg, dnul, 12, fmtstr, fmlen)
	call ftcol_defc (kdat, "RADESYS",   "A8", snul, cnul, 8, fmtstr, fmlen)
	call ftcol_defd (kdat, "EQUINOX", "F8.3",   "a", dnul,  8, fmtstr, fmlen)
	call ftcol_defd (kdat, "MJD-OBS", "F11.3",  "d", dnul, 11, fmtstr, fmlen)
	call ftcol_defr (kdat, "mag",      "F7.3", snul, rnul,  7, fmtstr, fmlen)
	call ftcol_defc (kdat, "pBand",      "A6", snul, cnul,  6, fmtstr, fmlen)
	call ftcol_defr (kdat, "RadVel", "F10.3",snul, rnul, 10, fmtstr, fmlen)
	call ftcol_defr (kdat, "MajAxis", "F9.2", "arcsec", rnul, 9, fmtstr, fmlen)
	call ftcol_defr (kdat, "MajAxPA", "F8.2",   deg, rnul, 8, fmtstr, fmlen)
	call ftcol_defr (kdat, "MinAxis", "F9.2", "arcsec", rnul, 9, fmtstr, fmlen)
	call ftcol_defr (kdat, "PM_RA",  "F9.4", "arcsec/a", rnul, 9, fmtstr, fmlen)
	call ftcol_defr (kdat, "PM_Dec", "F9.4", "arcsec/a", rnul, 9, fmtstr, fmlen)
	call ftcol_defr (kdat, "Parallax", "F7.4", "arcsec", rnul, 7, fmtstr, fmlen)
	call ftcol_defc (kdat, "ObjClass",  "A20", snul, cnul, 20, fmtstr, fmlen)
	call ftcol_defi (kdat, "CatFilePK",  "I6", snul, inul,  6, fmtstr, fmlen)

	call ftab_whead (fd, kdat, nselect)

# Boilerplate keywords 1
	call bpstd (fd, kdat, "ObjectCat", 1, 1, "ObjectId", snul)
	call bpadd (fd, kdat, 1, "CatFilePK", "CatFilePK", "CatFiles")

	call pkwstr (fd, kdat, "CATNAME", objfile, snul)

# Write out the OBJECT CAT TABLE
		
	do i = 0, ntarg-1 {
		if (SEL(tdat,i) == NO)
			next
		stat = sscan (DATLINE(tdat,i))
		call gargwrd (ident, SZ_ID)

		call sprintf (BUFROW(kdat), TWUSED(kdat), fmtstr)
			call pargi (INDEX(tdat,i))
			call pargstr (ident)
			call pargd   (RADTODEG(RA(tdat,i)))
			call pargd   (RADTODEG(DEC(tdat,i)))
			call pargstr ("")	# RADESYS
			call pargd   (STD_EQX(indat))	# XXX
			call pargr   (rnul)	#
			call pargr   (rnul)	# XXX mag
			call pargstr (cnul)	#pband
			call pargr   (rnul)	#
			call pargr   (rnul)	# maj ax
			if (PA(tdat,i) == INDEF)	# XXX should be specific
				call pargr (rnul)
			else
				call pargr (RADTODEG(PA(tdat,i)))
			call pargr   (rnul)	# min_ax
			call pargr   (rnul)	# PM_ra
			call pargr   (rnul)	# PM_dec
			call pargr   (rnul)	# Parallax
			if (PCODE(tdat,i) == CODE_GS)
				call pargstr ("Guide_Star")
			else if (PCODE(tdat,i) == CODE_AS)
				call pargstr ("Alignment_Star")
			else
				call pargstr ("Program_Target")
			call pargi   (1)
			
		call ftab_wrow (fd, kdat)
	}
	call ftab_free (fd, kdat)

call eprintf ("OK\n    table 2 ")
##### WRITE CatFiles TABLE
	call ftab_init (kdat, 2, fmtstr)

	call ftcol_defi (kdat, "CatFilePK",     "I6", snul, inul,   6, fmtstr, fmlen)
	call ftcol_defc (kdat, "CatFileName", "A255", snul, cnul, 255, fmtstr, fmlen)
#	call ftcol_defc (kdat, "CatFileName", "A68", snul, cnul, 68, fmtstr, fmlen)

	call ftab_whead (fd, kdat, 1)

# Boilerplate 2
	call bpstd (fd, kdat, "CatFiles", 1, 0, "CatFilePK", snul)

# Write out the OBJECT CAT TABLE
	call sprintf (BUFROW(kdat), TWUSED(kdat), fmtstr)
		call pargi (1)
		call pargstr (cnul)
			
	call ftab_wrow (fd, kdat)

	call ftab_free (fd, kdat)


call eprintf ("OK\n    table 3 ")
########  WRITE MaskDesign TABLE
	call ftab_init (kdat, 12, fmtstr)

	call ftcol_defi (kdat, "DesId",     "I11", snul, inul, 11, fmtstr, fmlen)
	call ftcol_defc (kdat, "DesName",   "A68", snul, cnul, 68, fmtstr, fmlen)
	call ftcol_defc (kdat, "DesAuth",   "A68", snul, cnul, 68, fmtstr, fmlen)
	call ftcol_defc (kdat, "DesCreat",  "A68", snul, cnul, 68, fmtstr, fmlen)
	call ftcol_defc (kdat, "DesDate",   "A19", snul, cnul, 19, fmtstr, fmlen)
	call ftcol_defd (kdat, "RA_PNT",  "F12.8", deg, dnul, 12, fmtstr, fmlen)
	call ftcol_defd (kdat, "DEC_PNT", "F12.8", deg, dnul, 12, fmtstr, fmlen)
	call ftcol_defc (kdat, "RADEPNT",    "A8", snul, cnul,  8, fmtstr, fmlen)
	call ftcol_defd (kdat, "EQUINPNT","F13.6",   "a", dnul, 13, fmtstr, fmlen)
	call ftcol_defd (kdat, "PA_PNT",  "F12.7", deg, dnul, 12, fmtstr, fmlen)
	call ftcol_defc (kdat, "DATE_PNT",  "A19", snul, cnul, 19, fmtstr, fmlen)
	call ftcol_defd (kdat, "LST_PNT",   "F8.3", deg, dnul,  8, fmtstr, fmlen)

	call ftab_whead (fd, kdat, 1)

# Boilerplate 3
	call bpstd (fd, kdat, "MaskDesign", 1, 0, "DesId", snul)

	call sprintf (BUFROW(kdat), TWUSED(kdat), fmtstr)
		call pargi (1)
		call pargstr ("output")	# XXX
		call pargstr (person)
		call pargstr (creatask)
		call pargstr (date)
		call pargd (RADTODEG(RA_TEL(indat)))
		call pargd (RADTODEG(DEC_TEL(indat)))
		call pargstr ("")
		call pargd (STD_EQX(indat))
		call pargd (RADTODEG(PA_ROT(indat)))
		call pargstr (date)		# XXX
		call pargd (RADTODEG(HA_TEL(indat)))
			
	call ftab_wrow (fd, kdat)

	call ftab_free (fd, kdat)


call eprintf ("OK\n    table 4 ")
########  WRITE DesiSlits TABLE

	call ftab_init (kdat, 9, fmtstr)

	call ftcol_defi (kdat, "dSlitId", "I11",  snul, inul, 11, fmtstr, fmlen)
	call ftcol_defi (kdat, "DesId",   "I11",  snul, inul, 11, fmtstr, fmlen)
	call ftcol_defr (kdat, "slitRA", "F12.8", deg, rnul, 12, fmtstr, fmlen)
	call ftcol_defd (kdat, "slitDec", "F12.8", deg, dnul, 12, fmtstr, fmlen)
	call ftcol_defc (kdat, "slitTyp",  "A1",  snul, snul,  1, fmtstr, fmlen)
	call ftcol_defr (kdat, "slitLen", "F11.3", "arcsec", rnul, 11, fmtstr, fmlen)
	call ftcol_defr (kdat, "slitLPA", "F8.3", deg, rnul,  8, fmtstr, fmlen)
	call ftcol_defr (kdat, "slitWid", "F11.3", "arcsec", rnul, 11, fmtstr, fmlen)
	call ftcol_defr (kdat, "slitWPA", "F8.3", deg, rnul,  8, fmtstr, fmlen)

	call ftab_whead (fd, kdat, nslit)

# Boilerplate 4
	call bpstd (fd, kdat, "DesiSlits", 1, 1, "dSlitID", "")
	call bpadd (fd, kdat, 1, "DesId", "DesId", "MaskDesign")

	do i = 0, nslit-1 {

		if (PCODE(sdat,i) == CODE_GS || PCODE(sdat,i) == CODE_RF)
			acode = 'G'
		else if (PCODE(sdat,i) == CODE_AS)
			acode = 'A'
		else
			acode = 'P'
		
		call sprintf (BUFROW(kdat), TWUSED(kdat), fmtstr)
			call pargi (INDEX(sdat,i))
			call pargi (1)
			call pargd (RADTODEG(RA(sdat,i)))
			call pargd (RADTODEG(DEC(sdat,i)))
			call pargc (acode)
			call pargr (LEN1(sdat,i)+LEN2(sdat,i))
			if (PA(sdat,i) == INDEF)	# XXX should be specific
				call pargr (0.)
			else
				call pargr (RADTODEG(PA(sdat,i)))
			call pargr (SLWID(sdat,i))
			call pargd (RADTODEG(PA_ROT(indat))+90.)	# XXX
			
		call ftab_wrow (fd, kdat)
	}

	call ftab_free (fd, kdat)


call eprintf ("OK\n    table 5 ")
########  WRITE SlitObjMap TABLE

	call ftab_init (kdat, 3, fmtstr)

	call ftcol_defi (kdat, "DesId",    "I11", snul, inul, 11, fmtstr, fmlen)
	call ftcol_defi (kdat, "ObjectId", "I11", snul, inul, 11, fmtstr, fmlen)
	call ftcol_defi (kdat, "dSlitId",  "I11", snul, inul, 11, fmtstr, fmlen)

	call ftab_whead (fd, kdat, nselect)

# Boilerplate 5
	call bpstd (fd, kdat, "SlitObjMap", 2, 3, "DesId", "ObjectId")
	call bpadd (fd, kdat, 1, "DesId", "DesId", "MaskDesign")
	call bpadd (fd, kdat, 2, "dSlitId", "dSlitId", "DesiSlits")
	call bpadd (fd, kdat, 3, "ObjectId", "ObjectId", "ObjectCat")


	do i = 0, ntarg-1 {
		if (SEL(tdat,i) == NO)
			next
	    call sprintf (BUFROW(kdat), TWUSED(kdat), fmtstr)
		call pargi (1)
		call pargi (i)
		call pargi (SLNDX(tdat,i))
			
	    call ftab_wrow (fd, kdat)
	}

	call ftab_free (fd, kdat)


call eprintf ("OK\n    table 6 ")
########  WRITE MaskBlu TABLE
	call ftab_init (kdat, 16, fmtstr)

	call ftcol_defi (kdat, "BluId",     "I11", snul, inul, 11, fmtstr, fmlen)
	call ftcol_defi (kdat, "DesId",     "I11", snul, inul, 11, fmtstr, fmlen)
	call ftcol_defc (kdat, "BluName",   "A68", snul, cnul, 68, fmtstr, fmlen)
	call ftcol_defc (kdat, "BluObsvr",  "A68", snul, cnul, 68, fmtstr, fmlen)
	call ftcol_defc (kdat, "BluCreat",  "A68", snul, cnul, 68, fmtstr, fmlen)
	call ftcol_defc (kdat, "BluDate",   "A19", snul, cnul, 19, fmtstr, fmlen)
	call ftcol_defd (kdat, "LST_Use",  "F8.3",deg, dnul,  8, fmtstr, fmlen)
	call ftcol_defc (kdat, "Date_Use",  "A19", snul, cnul, 19, fmtstr, fmlen)
	call ftcol_defc (kdat, "TELESCOP",  "A68", snul, cnul, 68, fmtstr, fmlen)
	call ftcol_defc (kdat, "RefrAlg",   "A68", snul, cnul, 68, fmtstr, fmlen)
	call ftcol_defd (kdat, "AtmTempC", "F5.1","degC", dnul, 5, fmtstr, fmlen)
	call ftcol_defd (kdat, "AtmPres",  "F6.1","mmHg", dnul, 6, fmtstr, fmlen)
	call ftcol_defd (kdat, "AtmHumid", "F5.3", "K/m", dnul, 5, fmtstr, fmlen)
	call ftcol_defd (kdat, "AtmTTLap", "F7.5", snul, dnul,  7, fmtstr, fmlen)
	call ftcol_defd (kdat, "RefWave",  "F7.2", "nm", dnul,  7, fmtstr, fmlen)
	call ftcol_defc (kdat, "DistMeth",  "A68", snul, cnul, 68, fmtstr, fmlen)

	call ftab_whead (fd, kdat, 1)

# Boilerplate 6
	call bpstd (fd, kdat, "MaskBlu", 1, 1, "BluId", snul)
	call bpadd (fd, kdat, 1, "DesId", "DesId", "MaskDesign")

	call sprintf (BUFROW(kdat), TWUSED(kdat), fmtstr)
		call pargi (1)
		call pargi (1)
		call pargstr ("output")		# XXX
		call pargstr (person)
		call pargstr (creatask)
		call pargstr (date)
		call pargd   (RADTODEG(HA_TEL(indat)+RA_TEL(indat)))
		call pargstr (date)		# Date_use XXX
		call pargstr ("Keck II")
		call pargstr (cnul)		# RefrAlg
		call pargd   (TEMP(indat))
		call pargd   (PRES(indat))
		call pargd   (0.)		# TMP: no clue  AtmHumid XXX
		call pargd   (0.)		# TMP: no clue  AtmTTLap XXX
		call pargd   (WAVER(indat))
		call pargstr (cnul)		# DistAlg
			
	call ftab_wrow (fd, kdat)

	call ftab_free (fd, kdat)


call eprintf ("OK\n    table 7 ")
##### write the BluSlits table
	call ftab_init (kdat, 11, fmtstr)

	call ftcol_def (kdat, "bSlitId", "I11", snul, inul, 11, fmtstr, fmlen)
	call ftcol_def (kdat, "BluId",   "I11", snul, inul, 11, fmtstr, fmlen)
	call ftcol_def (kdat, "dSlitId", "I11", snul, inul, 11, fmtstr, fmlen)
	call ftcol_def (kdat, "slitX1", fmt_mm, mm, dnul, 9, fmtstr, fmlen)
	call ftcol_def (kdat, "slitY1", fmt_mm, mm, dnul, 9, fmtstr, fmlen)
	call ftcol_def (kdat, "slitX2", fmt_mm, mm, dnul, 9, fmtstr, fmlen)
	call ftcol_def (kdat, "slitY2", fmt_mm, mm, dnul, 9, fmtstr, fmlen)
	call ftcol_def (kdat, "slitX3", fmt_mm, mm, dnul, 9, fmtstr, fmlen)
	call ftcol_def (kdat, "slitY3", fmt_mm, mm, dnul, 9, fmtstr, fmlen)
	call ftcol_def (kdat, "slitX4", fmt_mm, mm, dnul, 9, fmtstr, fmlen)
	call ftcol_def (kdat, "slitY4", fmt_mm, mm, dnul, 9, fmtstr, fmlen)

	call ftab_whead (fd, kdat, nslit)

# Boilerplate 7
	call bpstd (fd, kdat, "BluSlits", 1, 2, "bSlitId", "")
	call bpadd (fd, kdat, 1, "BluId", "BluId", "MaskBlu")
	call bpadd (fd, kdat, 2, "dSlitId", "dSlitId", "DesiSlits")

	do i = 0, nslit-1 {
		call sprintf (BUFROW(kdat), TWUSED(kdat), fmtstr)
			call pargi (i)			# XXX
			call pargi (1)			# XXX
			call pargi (i)			# XXX
			call pargd (XMM1(sdat,i))
			call pargd (YMM1(sdat,i))
			call pargd (XMM2(sdat,i))
			call pargd (YMM2(sdat,i))
			call pargd (XMM3(sdat,i))
			call pargd (YMM3(sdat,i))
			call pargd (XMM4(sdat,i))
			call pargd (YMM4(sdat,i))
		call ftab_wrow (fd, kdat)
	}

	call ftab_free (fd, kdat)


call eprintf ("OK\n    table 8 ")
### append the RDBmap:
	call rdb_map (fd)


call eprintf ("OK\n")

end

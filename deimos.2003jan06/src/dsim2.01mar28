# TBD:  get various params into proper routines

include	<math.h>
# include	<gset.h>
# include	<gim.h>
include	<time.h>
include	"dsimulator.h"
include	"deimos.h"

#
# int procedure	line_count (fd)
# procedure	data_init (indat)
# procedure	targ_init (fd, tdat, ntarg, indat)
# procedure	tel_coords (tdat, ntarg, indat)
# int procedure	chk_stat (x, y)

#
# LINE_COUNT: (I know I have another version of this around somewhere ...)
#

int procedure	line_count (fd)

pointer	fd			# file descriptor of open file

char	tchar
int	ndx

int	fscan(), nscan()
begin

# Count the entries
	ndx = 0
	while (fscan (fd) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#' || nscan() == 0) {
			next
		}
		ndx = ndx + 1
	}
	call seek (fd, BOF)

	return (ndx)
end


#
# DATA_INIT: initialize data structure for telescope/background data
# 

procedure	data_init (indat)

pointer	indat

double	clgetd ()
real	clgetr()

begin

# Allocate the data structure and vectors
	call malloc (indat, NINDAT, TY_STRUCT)
	call malloc (PTTELDAT(indat), NTELPAR, TY_DOUBLE)
	call malloc (PTDEFDAT(indat), NDEFPAR, TY_REAL)

# Read in params
	RA_FLD(indat)  = DEGTORAD(15. * clgetd ("ra0"))
	DEC_FLD(indat) = DEGTORAD(clgetd ("dec0"))
	PA_ROT(indat) = DEGTORAD(clgetd ("PA0"))
	HA_TEL(indat) = DEGTORAD(clgetr ("ha"))
	TEMP(indat) = clgetr ("temp")
	PRES(indat) = clgetr ("pressure")
	WAVER(indat) = clgetr ("lambda_cen") * 1.e-4		# in microns
	SLIT_GAP(indat) = clgetr ("sep_slit")		# sep. bet. slits, asec
	DEF_HLEN(indat) = 0.5 * clgetr ("min_slit")	# min. length in arcsec
	DEF_BOXR(indat) = 0.5 * clgetr ("box_sz")	# box 1/2-length in arcs
	DEF_SLWID(indat) = clgetr ("slit_width")	# slit width in arcsec
	STD_EQX(indat) = clgetr ("equinox")
#	epoch = clgetr ("epoch")


# TMP XXX TEST ONLY
	PAR_ANG(indat) = DEGTORAD (60.)

end


#
# TARG_INIT: initialize data structure for targets; fill
# 

procedure	targ_init (fd, tdat, ntarg, indat)

pointer	fd
pointer	tdat
int	ntarg
pointer	indat

real	eqnx_std		# TMP, needs input

char	tchar
char	idstr[SZ_ID]
int	prior, nlist, selcode
real	equinox, pangle, l1, l2
double	alpha, delta

int	ndx, n

int	fscan(), nscan()
int	line_count()

begin

# Count the entries
	ntarg = line_count (fd)

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


# Read in data
	ndx = 0
	while (fscan (fd) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#' || nscan() == 0) {
			next
		}
		call reset_scan ()

		call gargwrd (idstr, SZ_ID)
		call gargd (alpha)
		call gargd (delta)
		call gargr (equinox)
		call gargi (prior)

		if (nscan() < 5) {
			call eprintf ("Bad format on input line -- skipped\n`%s'\n")
			next
		}

		call gargi (nlist)
		call gargi (selcode)

		call gargr (pangle)
		call gargr (l1)
		call gargr (l2)

# Supply default arguments as needed:
		n = nscan()

# Are lengths present (either both or neither)
		if (n < 10) {
			l1 = DEF_HLEN(indat)
			l2 = DEF_HLEN(indat)
		}
		if (n < 8)
			pangle = INDEF
		if (n < 7)
			selcode = NO
		if (n < 6)
			nlist = PRIMARY

# Check for some special cases:
# XXX should check that nlist = reasonable value

		if (prior == CODE_AS) {
			l1 = DEF_BOXR(indat)
			l2 = DEF_BOXR(indat)
		}

		call reset_scan ()
		call gargstr (DATLINE(tdat,ndx), SZ_LINE-1)

# ASSUME no proper motion updates needed, for now, but add here in future
# XXX  will need epoch, epoch_std

# put on standard equinox, if needed
#	eqnx_std = STD_EQX(indat)
# XXX	if (equinox != eqnx_std)
# XXX		call precess (RA(tdat,ndx), DEC(tdat,ndx), equinox, eqnx_std)

# refract coordinates (ASSUME PA does not need refraction); to first order
# this is true since PA's are relative to rotator PA, similarly affected.
# It is unclear to me that this is necessary _now_ ...
#	ha = HA_TEL(indat)
# XXX	call refract (RA(tdat,ndx), DEC(tdat,ndx), ha)


# Assign values, always in radians for angles; arcsec for short lengths (???)
		INDEX(tdat,ndx) = ndx
		RA(tdat,ndx) = DEGTORAD(alpha*15.)
		DEC(tdat,ndx) = DEGTORAD(delta)
		if (pangle != INDEF)
			PA(tdat,ndx) = DEGTORAD(pangle)
		else
			PA(tdat,ndx) = INDEF
		PCODE(tdat,ndx) = prior
		LEN1(tdat,ndx) = l1
		LEN2(tdat,ndx) = l2
		SAMPL(tdat,ndx) = nlist
		SEL(tdat,ndx) = selcode
			
# XXX Assign slit-width
		if (PCODE(tdat,ndx) == CODE_AS) {
			PA(tdat,ndx) = INDEF
			RELPA(tdat,ndx) = INDEF
			SLWID(tdat,ndx) = 2.*DEF_BOXR(indat)
		} else {
			SLWID(tdat,ndx) = DEF_SLWID(indat)
		}

		ndx = ndx + 1
	}
	ntarg = ndx

# refract coordinates of telescope axis
# XXX	call refract (RA_TEL(indet), DEC_TEL(indat), ha)

end


#
# TEL_COORDS: Convert (refracted) alpha,dec into offsets from telescope center.
#

procedure	tel_coords (tdat, ntarg, indat)

pointer	tdat
int	ntarg
pointer	indat

int	i
double	r			# radius of object from tel-axis
double	p			# PA on sky from tel-axis
double	dec_obj, del_ra
double	cosr, sinp, cosp
double	ra0, dec0, pa0		# RA, Dec and PA on axis

real	x, y

int	chk_stat()
begin
	ra0  = RA_TEL(indat)
	dec0 = DEC_TEL(indat)
	pa0  = PA_ROT(indat)

	do i = 0, ntarg-1 {
		dec_obj = DEC(tdat,i)
		del_ra = RA(tdat,i) - ra0
		cosr = sin (dec_obj) * sin (dec0) +
			cos (dec_obj) * cos (dec0) * cos (del_ra)
		r = acos (cosr)

		sinp = cos (dec_obj) * sin (del_ra) / sqrt (1. - cosr*cosr)
		cosp = sqrt (max ((1. - sinp*sinp), 0.))
		if (dec_obj < dec0)
			cosp = -cosp
		p = atan2 (sinp, cosp)

# For now, convert radii to arcsec XXX
# XXX NB: I am not sure this is correct!  We should still be on SPHERICAL surf.
# 	but these are EUCLIDEAN relations.  Options: work in spherical coord
#	OR work in tan projection.
# More:  The difference at 10 arcmin between the tan and angle is < 0.002 arcsec
# If we convert "r" to tan(r) we should have the tan projection
		r = tan(r) * 206264.8
#		r = RADTODEG(r) * 3600.
		XARCS(tdat,i) = r * cos (p - pa0)
		YARCS(tdat,i) = r * sin (p - pa0)
		if (PA(tdat,i) == INDEF)
			RELPA(tdat,i) = INDEF
		else
			RELPA(tdat,i) = PA(tdat,i) - pa0

# Calc STAT
		x = XARCS(tdat,i)	# unclear if XYARCS will be real or dbl
		y = YARCS(tdat,i)
		STAT(tdat,i) = chk_stat (x, y)
	}
end

#
## CHK_STAT: is object within the slitmask? REAL INPUTS; INSTRUMENT SPECIFIC
## This could also be replaced by a bunch of limiting curves
#

int procedure	chk_stat (x, y)

pointer	tdat
int	ntarg
real	x, y

real	r

begin
	r = sqrt (x*x + y*y)

# Is object within 10 arcmin radius?
	if (r > 600.)
		return (NO)

# inner edge of mask
	if (y < 3*60.)			# TMP HARDCODES, test only
		return (NO)

# outer edge of mask
	if (y > 8*60.)			# TMP HARDCODES, test only
		return (NO)

# outer edge of mask
	if (x > XUPP_LIM || x < XLOW_LIM)
		return (NO)

# within radius of camera obscuration/vignetting?
	if (x*x+(y-YCAMCEN)**2 < RADVIGN**2)
		return (NO)
	

# near gaps in mosaic?
# XXX needs definition;  things like this show be contained in defines and
# limits (how close) in parameters

	if (abs (x+250.) < 4. || abs (x) < 4. || abs(x-250.) < 4.)
		return (NO)

# appears OK...
	return (YES)
end

	


#
# GEN_SLITS: initialize data structure for slits; fill
# 

procedure	gen_slits (tdat, ntarg, sdat, nslit, indat)

pointer	tdat
int	ntarg
pointer	sdat
int	nslit
pointer	indat

int	ndx, i
real	tana

begin

# Count the selected targets
	nslit = 0
	do i = 0, ntarg-1 {
		if (SEL(tdat,i) == YES)
			nslit = nslit + 1
	}

# Allocate the vectors: note that many quantities match the tdat structure
	call malloc (sdat, SDATLEN, TY_STRUCT)
	call malloc (PTINDEX(sdat), nslit, TY_INT)
	call malloc (PTRA(sdat), nslit, TY_DOUBLE)
	call malloc (PTDEC(sdat), nslit, TY_DOUBLE)
	call malloc (PTPA(sdat), nslit, TY_REAL)
	call malloc (PTLEN1(sdat), nslit, TY_REAL)
	call malloc (PTLEN2(sdat), nslit, TY_REAL)
	call malloc (PTWID(sdat), nslit, TY_REAL)
	call malloc (PTPCODE(sdat), nslit, TY_INT)
	call malloc (PTXARCS(sdat), nslit, TY_REAL)
	call malloc (PTYARCS(sdat), nslit, TY_REAL)
	call malloc (PTRELPA(sdat), nslit, TY_REAL)
	call malloc (PTSTAT(sdat), nslit, TY_INT)

	call malloc (PTSCOOR(sdat), nslit*NSCOOR, TY_DOUBLE)

# Set up slits for selected targets:
	ndx = 0
	do i = 0, ntarg-1 {
		if (SEL(tdat,i) == YES)	{	# or != 0

			INDEX(sdat,ndx) = ndx 
			PA(sdat,ndx) = PA(tdat,i)
			RELPA(sdat,ndx) = RELPA(tdat,i)
			PCODE(sdat,ndx) = PCODE(tdat,i)

			if (RELPA(tdat,i) == INDEF)
				tana = 0.
			else
				tana = tan (RELPA(tdat,i))
## XXX actually, we want to work out the actual PA on the metal if we just
# want it to be ortho on the slitmask.
			X1(sdat,ndx) = XARCS(tdat,i) - LEN1(tdat,i)
			Y1(sdat,ndx) = YARCS(tdat,i) - LEN1(tdat,i) * tana
			X2(sdat,ndx) = XARCS(tdat,i) + LEN2(tdat,i)
			Y2(sdat,ndx) = YARCS(tdat,i) + LEN2(tdat,i) * tana
# XXX NB: until the final sky_coords are calc'd, want X/YARCS to repr. objects
			XARCS(sdat,ndx) = XARCS(tdat,i)
			YARCS(sdat,ndx) = YARCS(tdat,i)
# XXX cuidado!  I am not sure that the tan-projection of the rel PA is the
# same as the rel PA -- MUST CHECK!

			SLWID(sdat,ndx) = SLWID(tdat,i)

# This is where we also assign slit index to object
			SLNDX(tdat,i) = ndx

			ndx = ndx + 1
		}
	}

# OK, we have assigned a slit to each object -- now ... XXX

# Sort  in x
# report conflicts
# display and allow binding ...
			
end


#
# BPSTD: Standard Boilerplate
# BPADD: Additional Boilerplate
#

procedure	bpstd (fd, kdat, tabname, npk, nfk, nam1, nam2, nam3)

pointer	fd		# file descriptor
pointer	kdat		# 
char	tabname[ARB]	# tablename
int	npk		# no. of primary keys
int	nfk		# no. of foreign keys
char	nam1[ARB], nam2[ARB], nam3[ARB]

# char	pnam1[ARB]	# first primary name
# char	pnam2[ARB]	# second primary name

char	date[SZ_TIME]					# time/date of file

int	nkey

char	snul
char	creatask[60]
char	person[30]

begin
# XXX
	call strcpy ("", snul, 1)
	call strcpy ("Phillips <phillips@ucolick.org>", person, 35)
	call strcpy ("DSIMULATOR -- June2000", creatask, 35)

	call ft_date (date, SZ_TIME)
# XXX

	call pkwstr (fd, kdat, "EXTNAME", tabname, snul)
	call pkwi   (fd, kdat, "EXTVER", 0, snul)
	call pkwstr (fd, kdat, "DATE", date, snul)
	call pkwstr (fd, kdat, "AUTHOR", person, snul)
	call pkwstr (fd, kdat, "CREATOR", creatask, snul)

	call pkwi   (fd, kdat, "NPRIKEY", npk, snul)

	if (nfk > 0) {
		call pkwi   (fd, kdat, "NFORKEY", nfk, snul)
		call pkwi   (fd, kdat, "NFORHDU", nfk, snul)
	}

	call pkwstr (fd, kdat, "PKTYP1", nam1, snul)
	if (npk == 2)
		call pkwstr (fd, kdat, "PKTYP2", nam2, snul)

	return

entry	bpadd (fd, kdat, nkey, nam1, nam2, nam3)

	if (nkey == 1) {
		call pkwstr (fd, kdat, "FKTYP1", nam1, snul)
		call pkwstr (fd, kdat, "FKFYP1", nam2, snul)
		call pkwi   (fd, kdat, "FKHDU1", nkey, snul)
		call pkwstr (fd, kdat, "HDLOC1", snul, snul)
		call pkwstr (fd, kdat, "HDXTN1", "TABLE", snul)
		call pkwstr (fd, kdat, "HDNAM1", nam3, snul)
	} else if (nkey == 2) {
		call pkwstr (fd, kdat, "FKTYP2", nam1, snul)
		call pkwstr (fd, kdat, "FKFYP2", nam2, snul)
		call pkwi   (fd, kdat, "FKHDU2", nkey, snul)
		call pkwstr (fd, kdat, "HDLOC2", snul, snul)
		call pkwstr (fd, kdat, "HDXTN2", "TABLE", snul)
		call pkwstr (fd, kdat, "HDNAM2", nam3, snul)
	} else if (nkey == 3) {
		call pkwstr (fd, kdat, "FKTYP3", nam1, snul)
		call pkwstr (fd, kdat, "FKFYP3", nam2, snul)
		call pkwi   (fd, kdat, "FKHDU3", nkey, snul)
		call pkwstr (fd, kdat, "HDLOC3", snul, snul)
		call pkwstr (fd, kdat, "HDXTN3", "TABLE", snul)
		call pkwstr (fd, kdat, "HDNAM3", nam3, snul)
	}

	return
end


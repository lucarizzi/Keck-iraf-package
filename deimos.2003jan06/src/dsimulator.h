####### Definitions for DSIMULATOR -- DEIMOS slitmask design

# May want to replace STAT with "OK"; then can use YES/NO unambiguously

define	SZ_ID		16	## TMP? reconcile elsewhere
define	PRIMARY		 1	# TMP? resolve how to name lists

# Define the struct for the target data
define	TDATLEN		16
define	PTINDEX		Memi[$1]
define	PTRA		Memi[$1+1]
define	PTDEC		Memi[$1+2]
define	PTPA		Memi[$1+3]
define	PTLEN1		Memi[$1+4]
define	PTLEN2		Memi[$1+5]
define	PTWID		Memi[$1+6]
define	PTPCODE		Memi[$1+7]
define	PTXARCS		Memi[$1+8]
define	PTYARCS		Memi[$1+9]
define	PTRELPA		Memi[$1+10]
define	PTSTAT		Memi[$1+11]
define	PTSAMPL		Memi[$1+12]
define	PTSEL		Memi[$1+13]
define	PTSLNDX		Memi[$1+14]
define	PTLINE		Memi[$1+15]

# Note that the number of targets can also be described here.

define	INDEX		Memi[PTINDEX($1)+$2]		# index
define	RA		Memd[PTRA($1)+$2]		# RA (rad)
define	DEC		Memd[PTDEC($1)+$2]		# Dec (rad)
define	PA		Memr[PTPA($1)+$2]		# PA on sky (rad)
define	LEN1		Memr[PTLEN1($1)+$2]		# length 1
define	LEN2		Memr[PTLEN2($1)+$2]		# length 2
define	SLWID		Memr[PTWID($1)+$2]		# width
define	PCODE		Memi[PTPCODE($1)+$2]		# pcode
define	XARCS		Memr[PTXARCS($1)+$2]		# X in arcsec
define	YARCS		Memr[PTYARCS($1)+$2]		# Y in arcsec
define	RELPA		Memr[PTRELPA($1)+$2]		# Rel PA wrt tel (rad)
define	STAT		Memi[PTSTAT($1)+$2]		# OK for selection?
define	SAMPL		Memi[PTSAMPL($1)+$2]		# sample code
define	SEL		Memi[PTSEL($1)+$2]		# Selected?
define	SLNDX		Memi[PTSLNDX($1)+$2]		# index to slit
define	DATLINE		Memc[PTLINE($1)+$2*SZ_LINE]	# full data line


# We want to bundle other items here:
# ra_tel, dec_tel, eqnx_std, pa_rot, ha0, temp, pressure

define	NINDAT		2
define	PTTELDAT	Memi[$1]
define	PTDEFDAT	Memi[$1+1]

define	NTELPAR		11
define	RA_TEL		Memd[PTTELDAT($1)]	# RA of tel. axis (rad)
define	DEC_TEL		Memd[PTTELDAT($1)+1]	# Dec of tel. axis (rad)
define	HA_TEL		Memd[PTTELDAT($1)+2]	# HA of tel. axis (rad)
define	PA_ROT		Memd[PTTELDAT($1)+3]	# PA of rotator (rad)
define	RA_FLD		Memd[PTTELDAT($1)+4]	# RA of fld (rad)
define	DEC_FLD		Memd[PTTELDAT($1)+5]	# dec of fld (rad)
define	PAR_ANG		Memd[PTTELDAT($1)+6]	# parallactic angle (rad)
define	STD_EQX		Memd[PTTELDAT($1)+7]	# standard equinox (yr)
define	TEMP		Memd[PTTELDAT($1)+8]	# temp (C)
define	PRES		Memd[PTTELDAT($1)+9]	# atm. press (mm Hg)
define	WAVER		Memd[PTTELDAT($1)+10]	# wavel. for refract. (microns)

define	NDEFPAR		4
define	SLIT_GAP	Memr[PTDEFDAT($1)]	# slit separation (asec)
define	DEF_HLEN	Memr[PTDEFDAT($1)+1]	# min slit 1/2-length (asec)
define	DEF_BOXR	Memr[PTDEFDAT($1)+2]	# box radius (asec)
define	DEF_SLWID	Memr[PTDEFDAT($1)+3]	# default slit width (asec)

# NB We ASSUME that PA_ROT is the on-axis PA of the rotator, || +X axis

# Definitions for slits: Note that many are same as for targets and use th
# same definition

define	SDATLEN		13
#define	PTINDEX		Memi[$1]	# defined above with tdat
#define	PTRA		Memi[$1+1]	# defined above with tdat
#define	PTDEC		Memi[$1+2]	# defined above with tdat
#define	PTPA		Memi[$1+3]	# defined above with tdat
#define	PTLEN1		Memi[$1+4]	# defined above with tdat
#define	PTLEN2		Memi[$1+5]	# defined above with tdat
#define	PTWID		Memi[$1+6]	# defined above with tdat
#define	PTPCODE		Memi[$1+7]	# defined above with tdat
#define	PTXARCS		Memi[$1+8]	# defined above with tdat
#define	PTYARCS		Memi[$1+9]	# defined above with tdat
#define	PTRELPA		Memi[$1+10]	# defined above with tdat
#define	PTSTAT		Memi[$1+11]	# defined above with tdat
define	PTSCOOR		Memi[$1+12]

define	NSCOOR		12			# number of bundled elements

define	X1		Memd[PTSCOOR($1)+$2*NSCOOR]	# double? # In asec
define	Y1		Memd[PTSCOOR($1)+$2*NSCOOR+1]	# double? # In asec
define	X2		Memd[PTSCOOR($1)+$2*NSCOOR+2]	# double? # In asec
define	Y2		Memd[PTSCOOR($1)+$2*NSCOOR+3]	# double? # In asec
define	XMM1		Memd[PTSCOOR($1)+$2*NSCOOR+4]	# double? # In mm
define	YMM1		Memd[PTSCOOR($1)+$2*NSCOOR+5]	# double? # In mm
define	XMM2		Memd[PTSCOOR($1)+$2*NSCOOR+6]	# double? # In mm
define	YMM2		Memd[PTSCOOR($1)+$2*NSCOOR+7]	# double? # In mm
define	XMM3		Memd[PTSCOOR($1)+$2*NSCOOR+8]	# double? # In mm
define	YMM3		Memd[PTSCOOR($1)+$2*NSCOOR+9]	# double? # In mm
define	XMM4		Memd[PTSCOOR($1)+$2*NSCOOR+10]	# double? # In mm
define	YMM4		Memd[PTSCOOR($1)+$2*NSCOOR+11]	# double? # In mm



define	KEYSFILE	"deimos$lib/dsimulator.keys"

define	RDB_MAP_FILE	"deimos$lib/RDBmap.hdu"


#############################################################################
# FOCAL PLANE OUTLINES
# This is pure drawing only, and is described by FP_FILE below.
# Need file of (x,y,pen) where pen 0=move,1=solid,2=dashed,3=dotted
# In principle, by simply changing the file we should end up with
# a new instrument description


define	MDATLEN		4
define	NFPDES		Memi[$1]
define	PTFPLX		Memi[$1+1]
define	PTFPLY		Memi[$1+2]
define	PTFPLZ		Memi[$1+3]		# Z is the pen action

define	FPLX		Memr[PTFPLX($1)+$2]
define	FPLY		Memr[PTFPLY($1)+$2]
define	FPLZ		Memi[PTFPLZ($1)+$2]	# Z is the pen action

#############################################################################

# INSTRUMENT SPECIFIC DEFINITIONS
# Adapting to new instrument includes not only changes here but in check_stat()
# several items needed here for general case:
# - inversion? Attempt to do with SENSX,SENSY
# - axis of dispersion?
# -- we could simply adopt that X will always be slit direction, y-dispersion

# This is the focal-plane outline in arcsec on sky. NB follows same convention
# of N=x-axis, E=y-axis

define	FP_FILE		"deimos$src/foc_plane.dat"

define	SENSX		-1.		# either 1 or -1
define	SENSY		 1.		# either 1 or -1
define	FLDCEN_X	0.		# Offset (asec) between tel & fld XXX
define	FLDCEN_Y	280.		# Offset (asec) between tel & fld XXX
define	XLOW_LIM	-500.		# lowest X value for slits
define	XUPP_LIM	500.		# highest X value for slits
define	YCAMCEN		700.		# camera center, in arcsec (11.66')
define	RADVIGN		302.		# radius of cam vign., asec (8.62")

#############################################################################

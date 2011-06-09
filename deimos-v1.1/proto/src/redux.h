# Define the struct for the REDUX data
define	RDATLEN		18
define	PTSLTNO		Memi[$1]
define	PTOBJNO		Memi[$1+1]
define	PTMAG		Memi[$1+2]
define	PTPBAND		Memi[$1+3]
define	PTXMM		Memi[$1+4]
define	PTYMM		Memi[$1+5]
define	PTX1MM		Memi[$1+6]
define	PTY1MM		Memi[$1+7]
define	PTX2MM		Memi[$1+8]
define	PTY2MM		Memi[$1+9]
define	PTXPX		Memi[$1+10]
define	PTYPX		Memi[$1+11]
define	PTX1PX		Memi[$1+12]
define	PTY1PX		Memi[$1+13]
define	PTX2PX		Memi[$1+14]
define	PTY2PX		Memi[$1+15]
define	PTCHIP		Memi[$1+16]
define	PTNAME		Memi[$1+17]

define	SLITNO		Memi[PTSLTNO($1)+$2]		# dSlitId
define	OBJNO		Memi[PTOBJNO($1)+$2]		# ObjectId
define	MAG		Memr[PTMAG($1)+$2]		# Magnitude
define	PBAND		Memc[PTPBAND($1)+$2]		# Passband for magnitude
define	XMM		Memr[PTXMM($1)+$2]		# X slitmask (mm)
define	YMM		Memr[PTYMM($1)+$2]		# Y slitmask (mm)
define	X1MM		Memr[PTX1MM($1)+$2]		# X slitmask (mm)
define	Y1MM		Memr[PTY1MM($1)+$2]		# Y slitmask (mm)
define	X2MM		Memr[PTX2MM($1)+$2]		# X slitmask (mm)
define	Y2MM		Memr[PTY2MM($1)+$2]		# Y slitmask (mm)
define	XPX		Memr[PTXPX($1)+$2]		# X in pixels (object)
define	YPX		Memr[PTYPX($1)+$2]		# Y in pixels (object)
define	X1PX		Memr[PTX1PX($1)+$2]		# X in pixels
define	Y1PX		Memr[PTY1PX($1)+$2]		# Y in pixels
define	X2PX		Memr[PTX2PX($1)+$2]		# X in pixels
define	Y2PX		Memr[PTY2PX($1)+$2]		# Y in pixels
define	OBJCHIP		Memi[PTCHIP($1)+$2]		# Chip number for object
define	OBJNAME		Memc[PTNAME($1)+$2*SZ_ID]	# full data line

define	SZ_ID	41
define	SZ_GUINAM	8	# Max char in GUI Name of slitmask



# Define the struct for the CHIP data
define	CDATLEN		8
#define	PTSLTNO		Memi[$1]
define	PTOVERL		Memi[$1+1]
define	PTCV1		Memi[$1+2]
define	PTCV2		Memi[$1+3]
define	PTCVO		Memi[$1+4]
define	PTCVS		Memi[$1+5]
define	PTCVA		Memi[$1+6]
define	PTCVG		Memi[$1+7]

#define	SLITNO		Memi[PTSLTNO($1)+$2]		# dSlitId
define	OVERLP		Memi[PTOVERL($1)+$2]		# overlap flag
define	CV1		Memi[PTCV1($1)+$2]		# CV ptr for slit edge1
define	CV2		Memi[PTCV2($1)+$2]		# CV ptr for slit edge2
define	CVO		Memi[PTCVO($1)+$2]		# CV ptr for overlap div
define	CVS		Memi[PTCVS($1)+$2]		# CV ptr for sky fit
define	CVALP		Memi[PTCVA($1)+$2]		# CV ptr for inp. alpha
define	CVGAM		Memi[PTCVG($1)+$2]		# CV ptr for inp. gamma

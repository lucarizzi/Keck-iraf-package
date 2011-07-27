# Define the struct for the box data
define	BDATLEN		9
define	PTSLTNO		Memi[$1]
define	PTOBJNO		Memi[$1+1]
define	PTMAG		Memi[$1+2]
define	PTPBAND		Memi[$1+3]
define	PTXMM		Memi[$1+4]
define	PTYMM		Memi[$1+5]
define	PTXPX		Memi[$1+6]
define	PTYPX		Memi[$1+7]
define	PTNAME		Memi[$1+8]

define	SLITNO		Memi[PTSLTNO($1)+$2]		# dSlitId
define	OBJNO		Memi[PTOBJNO($1)+$2]		# ObjectId
define	MAG		Memr[PTMAG($1)+$2]		# Magnitude
define	PBAND		Memc[PTPBAND($1)+$2]		# Passband for magnitude
define	XMM		Memr[PTXMM($1)+$2]		# X slitmask (mm)
define	YMM		Memr[PTYMM($1)+$2]		# Y slitmask (mm)
define	XPX		Memr[PTXPX($1)+$2]		# X in pixels
define	YPX		Memr[PTYPX($1)+$2]		# Y in pixels
define	OBJNAME		Memc[PTNAME($1)+$2*SZ_ID]	# full data line

define	SZ_ID	41
define	SZ_GUINAM	8	# Max char in GUI Name of slitmask

define	DYNA_DIR "/net/polo/local/kroot/data/deiccd/dyna/"

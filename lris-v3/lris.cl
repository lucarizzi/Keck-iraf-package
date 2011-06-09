#{ LRIS
# 
# W. M. Keck Observatory Low-Resolution Imaging Spectrometer 
# Observing Package
#
# Author:
#   Andrew C. Phillips [phillips@ucolick.org]
#   Gregory D. Wirth [wirth@keck.hawaii.edu]
#
#-----------------------------------------------------------------------

# load required packages...
noao
imred
ccdred
fitsutil	# fxheader needed by multi2simple
obsutil		# need nmisc.starfocus for lstarfocus
tables		# need tables.tstat for lspecfocus
mscred		# needed for lrisdisplay
ctio		# needed for multi2simple
artdata		# need mkheader for multi2simple
lred		# new slitmask alignment package by ACP
utilities	# 'translit' needed for multi2simple
;

package lris, bin = keckbin$

#	task	lstarfocus	= "lris$src/lstarfocus.cl"
	task	lspecfocus	= "lris$src/lspecfocus.cl"
	task	check_boxes	= "lris$src/check_boxes.cl"
	task	qbox		= "lris$src/qbox.cl"
#	task	lccdproc	= "lris$src/lccdproc.cl"
#	task	lrisbias	= "lris$src/lrisbias.cl"
#	task	xdistcor	= "lris$src/xdistcor.cl"
	task	get_seeing	= "lris$src/get_seeing.cl"
	task	do_check_boxes	= "lris$src/do_check_boxes.cl"
	task	do_lbox		= "lris$src/do_lbox.cl"
	task	$cdata		= "lris$src/cdata.cl"
	task	$lastimageb	= "$rsh lrisserver lastimageb"
	task	$lastimager	= "$rsh lrisserver lastimager"
	task	$wfib		= "$rsh lrisserver wfib"
	task	$wfir		= "$rsh lrisserver wfir"
	task	lrisdisplay	= "lris$src/lrisdisplay.cl"
	task	multi2simple	= "lris$src/multi2simple.cl"
	task	count_hdu	= "lris$src/count_hdu.cl"
	task	$sed		= "$foreign"
	task	lris_add_keywords = "lris$src/lris_add_keywords.cl"
	task	test		= "lris$src/test.cl"

	cache sections

	hidetask lastimageb,lastimager,wfib,wfir,count_hdu,qbox,sed

	if ( !defvar("tmp"))
		set tmp = "/tmp"

	type "lris$motd"

clbye()

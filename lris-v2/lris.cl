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
obsutil		# need nmisc.starfocus for lstarfocus
tables		# need tables.tstat for lspecfocus
mscred		# needed for lrisdisplay
;

package lris, bin = keckbin$

task	maskalign,
	mboxfind,
	mshift,
	mshiftb,
	xbox,
	xboxb,
	l2process = "lris$src/x_lris.e"

	task	lstarfocus	= "lris$src/lstarfocus.cl"
	task	lspecfocus	= "lris$src/lspecfocus.cl"
	task	check_boxes	= "lris$src/check_boxes.cl"
	task	check_boxesb	= "lris$src/check_boxesb.cl"
	task	check_boxesb2	= "lris$src/check_boxesb2.cl"
	task	qbox		= "lris$src/qbox.cl"
	task	lccdproc	= "lris$src/lccdproc.cl"
	task	lrisbias	= "lris$src/lrisbias.cl"
	task	xdistcor	= "lris$src/xdistcor.cl"
	task	get_seeing	= "lris$src/get_seeing.cl"
	task	$do_check_boxesb= "lris$src/do_check_boxesb.cl"
	task	$do_xboxb	= "lris$src/do_xboxb.cl"
	task	$cdata		= "lris$src/cdata.cl"
	task	$lastimageb	= "$rsh lrisserver lastimageb"
	task	$lastimager	= "$rsh lrisserver lastimager"
	task	$wfib		= "$rsh lrisserver wfib"
	task	lrisdisplay	= "lris$src/lrisdisplay.cl"

	cache sections

	hidetask lastimageb,lastimager,wfib

	type "lris$motd"

clbye()

#{ LRIS2
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
obsutil   # need nmisc.starfocus for lstarfocus
;

package lris2, bin = keckbin$

task	maskalign,
	mboxfind,
	mshift,
	mshiftb,
	xbox,
	xboxb,
	l2process = "lris2$src/x_lris.e"

	task	lstarfocus	= "lris2$src/lstarfocus.cl"
	task	check_boxes	= "lris2$src/check_boxes.cl"
	task	check_boxesb	= "lris2$src/check_boxesb.cl"
	task	qbox		= "lris2$src/qbox.cl"
	task	lccdproc	= "lris2$src/lccdproc.cl"
	task	lrisbias	= "lris2$src/lrisbias.cl"
	task	xdistcor	= "lris2$src/xdistcor.cl"

	cache sections

	type "lris2$motd"

clbye()

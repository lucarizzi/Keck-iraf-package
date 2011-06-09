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
obsutil   # need nmisc.starfocus for lstarfocus
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
	task	check_boxes	= "lris$src/check_boxes.cl"
	task	check_boxesb	= "lris$src/check_boxesb.cl"
	task	qbox		= "lris$src/qbox.cl"
	task	lccdproc	= "lris$src/lccdproc.cl"
	task	lrisbias	= "lris$src/lrisbias.cl"
	task	xdistcor	= "lris$src/xdistcor.cl"

	cache sections

	type "lris$motd"

clbye()

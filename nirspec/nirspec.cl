# NIRSPEC
# 
# W. M. Keck Observatory 
# Near-Infrared Spectrometer (NIRSPEC) Reduction Package
#
# Type:
#   IRAF command language (CL) script
# 
# Author:
#   Gregory D. Wirth [wirth@keck.hawaii.edu]
#
# Modification history:
#   2001-Jan-10     GDW     Original version
#-----------------------------------------------------------------------

# load required packages...
noao
onedspec
twodspec
longslit
artdata		# mk1dspec	(skyplot)
apextract	# apall		(donspec)
tables motd-
ttools		# tsel		(skyplot)

package nirspec
	task xdistcor		= "nirspec$xdistcor.cl"
	task ydistcor		= "nirspec$ydistcor.cl"
	task skyinterp		= "nirspec$skyinterp.cl"
	task plotcuts		= "nirspec$plotcuts.cl"
	task mktracer		= "nirspec$mktracer.cl"
	task spec_s2n		= "nirspec$spec_s2n.cl"
	task skyplot		= "nirspec$skyplot.cl"
	task donspec		= "nirspec$donspec.cl"

	type "nirspec$motd"
clbye()

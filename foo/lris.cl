#{ LRIS
# 
# W. M. Keck Observatory 
# Low-Resolution Imaging Spectrometer (LRIS) Reduction Package
#
# Author:
#   Gregory D. Wirth [wirth@keck.hawaii.edu]
#
#-----------------------------------------------------------------------

# load required packages...
noao
imred
ccdred
;

package lris
    # subpackage...
    task $lrisobs.pkg		= "lrisobs$lrisobs.cl"

    # tasks...
    task lccdproc		= "lris$src/lccdproc.cl"
    task lrisbias		= "lris$src/lrisbias.cl"
    task xdistcor		= "lris$src/xdistcor.cl"

    cache sections

clbye()

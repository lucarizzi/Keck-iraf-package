#{ LWS
# 
# W. M. Keck Observatory Low-Resolution Imaging Spectrometer 
# Reduction Package
#
# Author:
#   Gregory D. Wirth [wirth@keck.hawaii.edu]
#
#-----------------------------------------------------------------------

package lws, bin = keckbin$
    task lwscoadd		= "lws$x_lws.e"
    type "lws$motd"
clbye()

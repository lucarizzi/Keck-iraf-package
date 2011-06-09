#-----------------------------------------------------------------------
procedure mktracer( input, output)
#-----------------------------------------------------------------------
# nirspec$mktracer
# 
# Purpose:
#   From a NIRSPEC low-D flatfield exposure, create an image which
#   can be used to remove y distortions.
# 
# Author:
#	Gregory D. Wirth, W. M. Keck Observatory
# 
# Modification history:
#	2001-Jan-03		GDW		Original version
#
# See accompanying file mktracer.html for full documentation.
#-----------------------------------------------------------------------
string  input   { prompt="Name of input image" }
string  output  { prompt="Name of output image" }

begin
    string  i_input     # internal copy of input
    string  i_output    # internal copy of output
    string  tmpprefix   # temporary file prefix
    string  smoothflat  # smoothed version of flat
    string  gradflat    # gradient version of smoothflat
    string  absflat     # absolute value of gradflat
    real    maxval      # peak value in image

    # prompt...
    i_input = input
    i_output = output

    # temp files..
    tmpprefix = "mktracer"
    smoothflat = mktemp( tmpprefix) + ".fits"
    gradflat = mktemp( tmpprefix) + ".fits"
    absflat = mktemp( tmpprefix) + ".fits"

    # smooth the flat...
    boxcar( i_input, smoothflat, 11, 11)

    # compute the gradient to find the edges of the illuminated region...
    gradient( smoothflat, gradflat, gradient="0", boundary="nearest")
    imdelete( smoothflat, ver-)

    # get the maximum value...
    imstat( gradflat, field="max", format-) | scan( maxval)

    # normalize the image to prevent subsequent overflow...
    imarith( gradflat, "/", maxval, gradflat)

    # square the image to make negative peaks positive and 
    # to tighten up the peaks...
    imexpr( "a*a", absflat, gradflat, dims="auto",
        intype="auto", outtype="auto", verb-)
    imdelete( gradflat, ver-)

    # smooth with a median filter to iron out low spots...
    median( absflat, absflat, 11, 1, verb-)

    # and smooth with a Gaussian to merge double peaks...
    gauss( absflat, i_output, sigma=10, ratio=0.01, theta=0., nsigma=3.)
    imdelete( absflat, ver-)

end

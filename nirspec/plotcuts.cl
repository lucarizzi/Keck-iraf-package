#-----------------------------------------------------------------------
procedure plotcuts( image)
#-----------------------------------------------------------------------
# nirspec$plotcuts
# 
# Purpose:
#	Overplot cuts from multiple places along a spectrum
# 
# Author:
#	Gregory D. Wirth, W. M. Keck Observatory
# 
# Modification history:
#	2001-Jan-03		GDW		Original version
#
# See accompanying file plotcuts.html for full documentation.
#-----------------------------------------------------------------------
string  image   { prompt="Image(s) to plot" }
int     y1      { prompt="First row to plot", min=1 }
int     y2      { prompt="Last row to plot", min=1 }
int     nsum=5  { prompt="Number of columns to coadd", min=1 }
struct  *ilist

begin
    string  tmpprefix   # prefix for temp files
    string  i_image     # image name
    string  inlist      # list of input files
    string  in          # current input image
    string  apertlist   # list of apertures
    string  sxn         # image section
    real    maxval      # maximum value in image section
    int     nx          # number of columns in image

    # prompt...
    i_image = image

    # tempfiles...  
    tmpprefix = "PC"
    inlist = mktemp( tmpprefix)

    # expand file template...
    sections( i_image, >inlist)

    # loop over input images....
    ilist = inlist
    while( fscan( ilist, in)!=EOF){

        # get image size...
        hselect( in, "naxis1", "yes") | scan( nx)

        # generate aperture list...
        apertlist = str(nx/2)
        apertlist = apertlist + "," + str(nsum)
        apertlist = apertlist + "," + str(nx-nsum)

        # determine peak pixel...
        sxn = "[*," + str(y1) + ":" + str(y2) + "]"
        imstat( in//sxn, fields="max", format-) | scan( maxval)
        maxval = maxval * nsum

        # generate the plot...
        onedspec.nsum = nsum
        specplot( in, apert=apertlist, xmin=y1, xmax=y2, autolay-, 
            ysca+, ymin=0., ymax=maxval, title=in)
    }

    # clean up...
    delete( inlist, ver-)
    ilist = ""

end

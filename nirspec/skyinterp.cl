#-----------------------------------------------------------------------
procedure skyinterp( inimage, outimage, skyimage)
#-----------------------------------------------------------------------
# nirspec$skyinterp
# 
# Purpose:
#   Given a 1-D spectrum, produce a "pretty" version in which all
#   regions dominated by night sky emission lines are replaced by a 
#   smoothed version of the underlying spectrum.
#
# Author:
#   Gregory D. Wirth, W. M. Keck Observatory
# 
# Modification history:
#   2000-Dec-28     GDW     Original version
#
# See accompanying file skyinterp.html for full documentation.
#-----------------------------------------------------------------------

string  inimage     { prompt="input image(s)" }
string  outimage    { prompt="output image(s)" }
string  skyimage    { prompt="image containing sky spectrum" }
real    thresh      { prompt="sky value above which to excise" }
bool    verbose=no  { prompt="print operations?" }
struct  *ilist
struct  *olist

begin
    string  tmpprefix   # prefix for temp files
    string  iinimage    # internal copy of inimage
    string  ioutimage   # internal copy of ioutimage
    string  iskyimage   # internal copy of skyimage
    string  sxn         # image section
    string  in          # current input image
    string  out         # current output image
    string  inlist      # list of input images
    string  outlist     # list of output images
    string  medianimage # image containing median data
    string  cleanimage  # image containing cleaned data
    int     nin         # number of input images
    int     nout        # number of output images
    int     istat       # status of scan
    int     n1=29       # width of median filter for first pass
    int     n2=9        # width of median filter for second pass
    int     naxis2      # number of apertures in image section
    int     naxis3      # number of bands in image section

    # queries...
    iinimage = inimage
    ioutimage = outimage
    iskyimage = skyimage

    # check for non-existent sky image...
    if( ! imaccess( iskyimage))
        error( 1, "  WARNING: sky image "//iskyimage//" does not exist!")

    # temp files...
    tmpprefix = "skyinterp"
    inlist = mktemp( tmpprefix)
    outlist = mktemp( tmpprefix)
    medianimage = mktemp( tmpprefix)
    cleanimage = mktemp( tmpprefix)

    # caches...
    cache( "sections")

    # expand input lists...
    sections( iinimage, >inlist)
    nin = sections.nimages
    sections( ioutimage, >outlist)
    nout = sections.nimages

    # check lists for errors...
    if( nin!=nout)
        error( 1, "different numbers of input and output images!")

    # check that the sky image is a section...
    hselect( iskyimage, "naxis2 naxis3", "yes") | scan( naxis2, naxis3)
    if( nscan() > 0)    
        error( 1, 'please specify an image section to use for sky')

    # loop over input/output pairs...
    ilist = inlist
    olist = outlist
    while( fscan( ilist, in)!=EOF && fscan( olist, out)!=EOF){

        # check for non-existent image...
        if( ! imaccess( in)){
            print( "  WARNING: "//in//" does not exist! -- skipping")
            next
        }

        # check for clobber...
        if( imaccess( out) && out!=in){
            print( "  WARNING: "//out//" already exists! -- skipping")
            next
        }

        # generate a cleaned image...
        median( in, medianimage, n1, 1, verbose=verbose)

        # expunge night sky lines...
        imexpr( "c<d ? a : b", cleanimage,
            in, medianimage, iskyimage, thresh, verbose=verbose)
        imdelete( medianimage, ver-)

        # generate an even cleaner image...
        median( cleanimage, medianimage, n2, 1, verbose=verbose)
        imdelete( cleanimage, ver-)
        
        # expunge night sky lines...
        imexpr( "c<d ? a : b", out,
            in, medianimage, iskyimage, thresh, verbose=verbose)
        imdelete( medianimage, ver-)
    }

    # clean up...
    delete( inlist, ver-)
    delete( outlist, ver-)
    ilist = ""
    olist = ""

end

#-----------------------------------------------------------------------
procedure spec_s2n( input)
#-----------------------------------------------------------------------
# nirspec$spec_s2n
# 
# Purpose:
#   Create simulated sky spectrum and corresponding line list
#
# Type:
#   IRAF command language (CL) script
# 
# Author:
#   Gregory D. Wirth, W. M. Keck Observatory
# 
# Note:
#   See accompanying file spec_s2n.html for full documentation
# 
# Modification history:
#   2001-Jan-09     GDW     Original version
#-----------------------------------------------------------------------

string  input   { prompt="Multispec images to measure"}
struct  *ilist

begin
    string  tmpprefix   # prefix for temp files
    string  i_input     # internal copy of input
    string  inlist      # file listing input images
    string  in          # current input image
    string  s2nimage    # temporary image containing s/n per pixel
    real    midpt       # median S/N per pixel
    real    sum=0       # sum of squares
    int     naxis2      # number of apertures in image
    int     naxis3      # number of bands in image
    int     i           # counter
    int     nin         # number of images

    # queries...
    i_input = input

    # temp files...
    tmpprefix = "S2N"
    inlist = mktemp( tmpprefix)
    s2nimage = mktemp( tmpprefix)

    # expand input template...
    sections( i_input, >inlist)
    nin = sections.nimages

    # loop over input images...
    ilist = inlist
    while( fscan( ilist, in)!=EOF){

        # check for non-existent image...
        if( ! imaccess( in)){
            print( "  WARNING: "//in//" does not exist; skipped")
            next
        }

        # get image sizes...
        hselect( in, "naxis2 naxis3", "yes") | scan( naxis2, naxis3)
        if( nscan() != 2){
            print( "  WARNING: Could not get image dimensions for image "//in//"; skipped")
            next
        }   

        # check for 4 bands...
        if( naxis3 != 4){
            print( "  WARNING: Image "//in//" does not have 4 bands;  skipped")
            next
        }

        # loop over apertures...
        for( i=1 ; i<=naxis2 ; i+=1 ){
            imarith( in//"[*,"//i//",1]", "/", in//"[*,"//i//",4]", s2nimage)
            imstat( s2nimage, field="midpt", format-) | scan( midpt)
            sum += midpt**2
            printf( "  Image: %s", in)
            if( naxis2 > 1)
                printf( "  Aperture: %d", i)
            printf( "  S/N: %.1f per pixel (median)", midpt)
            printf( "\n")
            imdelete( s2nimage, ver-)
        }

    }

    if( nin>1 && naxis2==1 ){
        sum = sqrt(sum)
        printf( "    Quadrature sum: %.1f per pixel\n", sum)
    }

    # clean up...
    delete( inlist, ver-)
    ilist = ""

end

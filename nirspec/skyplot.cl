#-----------------------------------------------------------------------
procedure skyplot( inlines, outlines, outimage)
#-----------------------------------------------------------------------
# nirspec$skyplot
# 
# Purpose:
#   Given a list of sky emission features and their relative
#   intensities, generate a spectrum of these lines over a specified
#   wavelength interval.  If desired, create a plot and generate a list
#   of the strongest lines in the spectrum.
# 
# Type:
#   IRAF command language (CL) script
# 
# Author:
#   Gregory D. Wirth, W. M. Keck Observatory
# 
# Note:
#   See accompanying file skyplot.html for full documentation
# 
# Modification history:
#   2001-Jan-09     GDW     Original version
#-----------------------------------------------------------------------

string  inlines="nirspec$lowd_ir_ohlines.dat" \
            {prompt="File listing sky lines"}
string  outlines    { prompt="Output line list [optional]" }
string  outimage    { prompt="Name of image to create [optional]" }
real    wavelength=10000. \
            {prompt="Central wavelength [Angstroms]", min=1000 }
real    dispersion=2.077    {prompt="Dispersion [Angstroms/px]", min=0.}
real    fwhm=5.     {prompt="Width of sky lines [px]", min=0.}
int     ncols=1024  {prompt="number of pixels in spectrum" }
bool    identify=yes    { prompt="Identify lines?" }
int     nlines=50   { prompt="Number of lines to identify" }
string  graphics="stdgraph" { prompt="Graphics outimage device" }


begin
    string  tmpprefix       # prefix for temp files
    string  i_outimage      # name of outimage image
    string  i_inlines       # internal copy of inlines
    string  i_outlines      # internal copy of outlines
    string  cursor          # input file for identify
    string  temp_table      # temporary table
    real    wstart          # starting wavelength
    real    wend            # ending wavelength
    bool    saveimage       # save the output image?
    bool    savelines       # save the output image?

    # prompt...
    i_inlines = inlines 
    i_outlines = outlines   
    i_outimage = outimage

    # temp files...
    tmpprefix = "SKYPLOT"
    cursor = mktemp( tmpprefix)
    temp_table = mktemp( tmpprefix)

    if( i_outimage == ""){
        i_outimage = mktemp( tmpprefix)
        saveimage = no
    } else {
        saveimage = yes
    }

    if( i_outlines == ""){
        i_outlines = mktemp( tmpprefix)
        savelines = no
    } else {
        savelines = yes
    }

    # check for data file...
    if( !access( i_inlines))
        error( 1, "cannot access linelist file "//i_inlines)

    # check whether to clobber output line list...
    if( i_outlines != "" && access( i_outlines)){
        delete( i_outlines, ver+, default+)
        if( access( i_outlines))
            error( 1, "Operation would clobber existing file "//i_outlines)
    }

    # check whether to clobber output image...
    if( i_outimage != "" && imaccess( i_outimage)){
        imdelete( i_outimage, ver+, default+)
        if( imaccess( i_outimage))
            error( 1, "Operation would clobber existing image "//i_outimage)
    }

    # generate simulated sky image...
    wstart = wavelength - 0.5*ncols*dispersion
    wend = wstart + ncols*dispersion
    print( "Creating...")
    mk1dspec( i_outimage, ap=1, rv=0., z=no, ncols=ncols, naps=1, 
        header="artdata$stdheader.dat", wstart=wstart, wend=wend, 
        continuum=0., slope=0., temperature=0., fnu-, 
        lines=i_inlines, nlines=0, profile="gaussian", peak=1, 
        gfwhm=fwhm, seed=1, comments+)

    # generate linelist...
    if( identify || savelines){

        # extract data from the stated wavelength range...
        printf( "c1>%d && c1<%d\n", wstart, wend) | scan( line)
        tselect( i_inlines, temp_table, line)

        # sort these by intensity...
        tsort( temp_table, "c2", asc-)

        # extract the first nlines...
        printf( "row() <= %d\n", nlines) | scan( line)
        tselect( temp_table, i_outlines, line)
        delete( temp_table)

        # sort new linelist by wavelength...
        tsort( i_outlines, "c1", asc+)

    }

    # plot it...
    if( identify){
        print( "y", >>cursor)
        print( ":label coo", >>cursor)
        print( "no") | \
            onedspec.identify( i_outimage, section="middle line", 
                coordlist=i_outlines, units="Angstroms", 
                nsum="1", match=-3., maxfeatures=nlines, zwidth=100., 
                ftype="emission", fwidth=fwhm, cradius=fwhm, threshold=0., 
                minsep=2., autowrite-, graphics=graphics, cursor=cursor,
                > "dev$null" )
        delete( cursor, ver-)
    }

    # clean up...
    if( ! saveimage)
        imdelete( i_outimage, ver-)
    if( ! savelines)
        delete( i_outlines, ver-)

end

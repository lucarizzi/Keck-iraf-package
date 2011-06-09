#-----------------------------------------------------------------------
procedure donspec( input, root, docheck)
#-----------------------------------------------------------------------
# nirspec$donspec
# 
# Purpose:
#   Coordinate all steps in the reduction of NIRSPEC low-dispersion data.
# 
# Author:
#	Gregory D. Wirth, W. M. Keck Observatory
# 
# Modification history:
#	2001-Jan-10		GDW		Original version
#
# See accompanying file donspec.html for full documentation.
#-----------------------------------------------------------------------
string  input       { prompt="List of 2-D images to reduce" }
string  root        { prompt="Rootname of output 1-D spectra" }
string  xref        { prompt="Name of flat image to use for xdistcor" }
string  yref=""     { prompt="Name of line image to use for ydistcor" }
string  coordlist   { prompt="File listing calibration lines" }
int     x1          { min=1, max=1024, prompt="Starting column to extract" }
int     x2          { min=1, max=1024, prompt="Ending column to extract" }
real    dy=0.       { prompt="Feature shift in y from x1 to x2" }
real    gain=5      { prompt="Detector gain", min=0. }
real    rdnoise=30  { prompt="Detector read noise", min=0. }
real    thresh=-1   { prompt="Remove pixels with sky above this level" }
string  logfile="logfile"   { prompt="Logfile" }
bool    verbose=yes         { prompt="Give feedback?" }
bool    docheck=yes         { prompt="Check results?" }
struct  *ilist

begin
    string  i_input     # internal copy of input
    string  i_root      # internal copy of root
    string  tmpprefix   # prefix for temp files
    string  masterlist  # master list of input files
    string  inlist      # list of input files
    string  outlist     # list of output files
    string  in          # current input image
    string  out         # current output image
    string  cursor      # cursor input for graphics tasks
    string  new_trace   # transformed version of trace image
    string  trace_image     # name of trace image to create
    string  combined_image  # name of combined image to create
    string  sky_image       # name of sky image to create
    string  dispcor_image   # name of dispersion-corrected image
    string  skyinterp_image # name of sky interpolated image
    string  ref         # reference image
    bool    yesno=yes   # prompt
    int     nin         # number of input images
    int     nout        # number of output images
    int     x           # column to plot
    int     y           # input value
    int     istat       # status of scan
    int     y1          # first row for plotcut zoom
    int     y2          # last row for plotcut zoom

    # prompt...
    i_input = input
    i_root = root

    # set up temp files...
    tmpprefix = "DONSPEC"
    masterlist = mktemp( tmpprefix)
    inlist = mktemp( tmpprefix)
    outlist = mktemp( tmpprefix)
    cursor = mktemp( tmpprefix)
    new_trace = mktemp( tmpprefix) // ".fits"

    # set up new image names...
    trace_image = i_root//"_trace"
    sky_image = i_root//"_sky"
    combined_image = i_root//"_avsigclip"
    dispcor_image = i_root//"_dispcor"
    skyinterp_image = i_root//"_skyinterp"

    # generate list of input spectra...
    files( i_input, >masterlist)
    count( masterlist) | scan( nin)
    if( nin < 1)
        error( 1, "no input spectra")

    # verify that all are acessible...
    ilist = masterlist
    while( fscan( ilist, in)!=EOF){
        if( ! imaccess( in))
            error( 1, "image "//in//" not accessible")
    }

    # verify reference images...
    if( ! imaccess( xref))
        error( 1, "image "//xref//" not accessible")
    if( yref!="" && ! imaccess( yref))
        error( 1, "image "//yref//" not accessible")

    # remove images which are in the way...
    if( imaccess( trace_image))
        imdelete( trace_image, ver+, def+)
    if( imaccess( trace_image))
        error( 1, "trace image "//trace_image//" already exists")

    if( imaccess( sky_image))
        imdelete( sky_image, ver+, def+)
    if( imaccess( sky_image))
        error( 1, "sky image "//sky_image//" already exists")

    if( imaccess( combined_image))
        imdelete( combined_image, ver+, def+)
    if( imaccess( combined_image))
        error( 1, "combined image "//combined_image//" already exists")

    if( imaccess( dispcor_image))
        imdelete( dispcor_image, ver+, def+)
    if( imaccess( dispcor_image))
        error( 1, "dispcor image "//dispcor_image//" already exists")

    if( imaccess( skyinterp_image))
        imdelete( skyinterp_image, ver+, def+)
    if( imaccess( skyinterp_image))
        error( 1, "skyinterp image "//skyinterp_image//" already exists")

    if( ! access( coordlist))
        error( 1, "coordlist file "//coordlist//" not found")

    # modify image headers...
    hedit( "@"//masterlist, "dispaxis", "2", add+, verify-, update+)

    # generate trace image...
    if( verbose)
        print( "  Generating trace image...")
    mktracer( xref, trace_image)

    # remove x distortion...
    if( verbose)
        print( "  Removing x distortion...")
    copy( masterlist, inlist)
    files( "@"//masterlist//"//%.fits%x%", >outlist)
    imdelete( "@"//outlist, ver+, def+)
    xdistcor( "@"//inlist, "@"//outlist, ref=trace_image, verbose=verbose,
        useold-)

    # allow checking of results...
    if( docheck){

        # transform trace image...
        hedit( trace_image, "dispaxis", "2", add+, verify-, update+)
        xdistcor( trace_image, new_trace, ref=trace_image, verbose=verbose,
            useold+)

        # plot results...
        print( ":a10", >>cursor)
        print( ":l512", >>cursor)
        print( ":o", >>cursor)
        print( ":l10", >>cursor)
        print( ":o", >>cursor)
        print( ":l1010", >>cursor)
        implot( new_trace, line=512, coords=cursor)

        # clean up...
        delete( cursor, ver-)
        imdelete( new_trace, ver-)
        print( "Check the above plot to verify alignment of features.")
        printf( "Continue processing? (yes): ")
        istat = scan( yesno)
        if( ! yesno)
            goto abort
    }

    # remove y distortion...
    if( verbose)
        print( "  Removing y distortion...")
    delete( inlist, ver-)
    copy( outlist, inlist)
    delete( outlist, ver-)
    files( "@"//masterlist//"//%.fits%y%", >outlist)
    imdelete( "@"//outlist, ver+, def+)

    # If yref is defined, use that image.  
    # If undefined, grab first image from input file list...
    if( yref == "")
        files( "@"//inlist) | scan( ref)
    else
        ref = yref

    # perform correction...
    ydistcor( "@"//inlist, "@"//outlist, ref=ref, verbose=verbose,
        x1=x1, x2=x2, dy=dy)
    imdelete( "@"//inlist, ver-)

    # allow checking of results...
    if( docheck){

        # generate a plot...
        files( "@"//outlist) | scan( ref)
        print( "q", >>cursor)
        x = 0.5*(x2-x1)
        splot( ref, x, cursor=cursor)

        # get plotting limits...
        print( "Please define one corner of zoom window by moving cursor and pressing any key")
        istat = fscan( gcur, y1)
        print( "Please define other corner of zoom window by moving cursor and pressing any key")
        istat = fscan( gcur, y2)

        # make plot...
        plotcuts( ref, y1=y1, y2=y2, nsum=5)

        # clean up...   
        delete( cursor, ver-)
        print( "Please verify alignment of emission features as plotted here.")
        printf( "Continue processing? (yes): ")
        istat = scan( yesno)
        if( ! yesno)
            goto abort
    }

    # extract...
    if( verbose)
        print( "  Extracting spectra...")
    delete( inlist, ver-)
    copy( outlist, inlist)
    delete( outlist, ver-)
    files( "@"//inlist//"//%.fits%.ms%", >outlist)
    imdelete( "@"//outlist, ver+, def+)
    apall( input="@"//inlist, nfind=1, output="@"//outlist, apertures="",
        format="multispec", references="", profiles="", interactive+, find+,
        recenter-, resize-, edit+, trace+, fittrace+, extract+, extras+,
        review-, line=INDEF, nsum=-25, lower=-10., upper=10., apidtable="",
        b_function="chebyshev", b_order=1, b_sample="-30:-6,6:30",
        b_naverage=-100, b_niterate=0, b_low_reject=3., b_high_reject=3.,
        b_grow=0., width=5., radius=10., threshold=0., minsep=5.,
        maxsep=1000., order="increasing", aprecenter="", npeaks=INDEF, shift+,
        llimit=INDEF, ulimit=INDEF, ylevel=0.1, peak+, bkg+, r_grow=0.,
        avglimits-, t_nsum=25, t_step=25, t_nlost=100, t_function="chebyshev",
        t_order=2, t_sample="*", t_naverage=1, t_niterate=10, t_low_reject=3.,
        t_high_reject=3., t_grow=0., background="fit", skybox=1,
        weights="variance", pfit="fit1d", clean+, saturation=INDEF,
        readnoise=rdnoise, gain=gain, lsigma=3., usigma=3., nsubaps=1)

    # get S/N...
    if( verbose)
        print( "  Checking S/N...")
    delete( inlist, ver-)
    copy( outlist, inlist)
    spec_s2n( "@"//inlist)
    if( logfile!="")
        spec_s2n( "@"//inlist, >>logfile)

    # combine spectra...
    if( verbose)
        print( "  Combining spectra...")
    scombine( "@"//inlist, combined_image, noutput="", logfile=logfile,
        apertures="", group="apertures", combine="average", 
        reject="avsigclip", first-, w1=INDEF, w2=INDEF, dw=INDEF,
        nw=INDEF, log-, scale="median", zero="none", weight="median", 
        sample="",
        lthreshold=INDEF, hthreshold=INDEF, nlow=1, nhigh=1, nkeep=1, mclip+,
        lsigma=3., hsigma=3., rdnoise=rdnoise, gain=gain, snoise="0.",
        sigscale=0.1, pclip=-0.5, grow=0, blank=0.)

    # generate sky image...
    if( verbose)
        print( "  Generating sky image...")
    files( "@"//inlist) | scan( ref)
    imcopy( ref//"[*,1,3]", sky_image, verbose-)

    # solve wavelength function...
    if( verbose)
        print( "  Identifying lines...")
    onedspec.identify( sky_image, section="", database="database",
        coordlist=coordlist, units="", nsum="1", match=-3., 
        maxfeatures=50, zwidth=100., ftype="emission", fwidth=4.,
        cradius=5., threshold=0., minsep=2., function="spline3", order=1,
        sample="*", niterate=0, low_reject=3., high_reject=3., grow=0., ,
        graphics="stdgraph", cursor="", aidpars="")

    # apply solution...
    if( verbose)
        print( "  Applying wavelength solution...")
    hedit( combined_image, "refspec1", sky_image, add+, upd+, ver-)
    dispcor( combined_image, dispcor_image, linearize-, database="database",
        table="", w1=INDEF, w2=INDEF, dw=INDEF, nw=INDEF, log-,
        flux+, samedisp-, global-, ignoreaps-, confirm-,
        listonly-, verbose+)

    # check sky level, as needed...
    if( thresh<0.){
        print( "q", >>cursor)
        splot( sky_image, x, cursor=cursor)

        # get plotting limits...
        print( "Defining sky threshold level...")
        istat = fscan( gcur, x, y)
    } else {
        y = thresh
    }

    # replace sky...
    if( verbose)
        print( "  Removing sky lines...")
    skyinterp( dispcor_image, skyinterp_image, sky_image, thresh=y,
        verbose=verbose)

    # print message...
    printf( "\tOutput trace image     = %s\n", trace_image)
    printf( "\tOutput sky image       = %s\n", sky_image)
    printf( "\tOutput combined image  = %s\n", combined_image)
    printf( "\tOutput dispcor image   = %s\n", dispcor_image)
    printf( "\tOutput skyinterp image = %s\n", skyinterp_image)
    printf( "Processing completed.\n")
    beep()

    abort:
    delete( inlist, ver-)
    delete( outlist, ver-)
    ilist = ""
    
end

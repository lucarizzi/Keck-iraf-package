#-----------------------------------------------------------------------
procedure xdistcor( input, output, ref, useold)
#-----------------------------------------------------------------------
# nirspec$xdistcor
# 
# Purpose:
#   Remove x distortion from NIRSPEC image
# 
# Author:
#   Gregory D. Wirth, W. M. Keck Observatory
# 
# Note:
#	Refer to accompanying file xdistcor.html for full documentation.
# 
# Modification history:
#   2001-Jan-04     GDW     Original version
#-----------------------------------------------------------------------

string  input   { prompt="Name of input image(s)" }
string  output  { prompt="Name of output image(s)" }
string  ref     { prompt="Name of reference (trace) image" }
bool    verbose { yes, prompt="Give feedback?" }
bool    useold  { yes, prompt="Use existing solution?" }

begin
    string  i_input     # internal copy of input
    string  i_output    # internal copy of output
    string  i_ref       # internal copy of ref
    string  tmpprefix   # temporary file prefix
    string  inlist      # file listing input images
    string  outlist     # file listing output images
    string  logfile     # output file for tasks
    string  database    # name of directory storing solutions
    string  solution    # name of solution file
    string  cursor      # file for graphics cursor input to identify
    real    shift       # feature shift bewteen consecutive fits
    real    dx=0        # Column shift of features from bottom to top
    int     nrows       # number of rows in image
    int     nstep=50    # number of rows to step between fits
    int     nsum=5      # number of rows to sum
    bool    skipfit=no  # prompt whether to skip fitting

    # prompt...
    i_input = input
    i_output = output
    i_ref = ref

    # temp files..
    tmpprefix = "XDC"
    inlist = mktemp( tmpprefix)
    outlist = mktemp( tmpprefix)
    cursor = mktemp( tmpprefix)

    # expand input lists...
    sections( i_input, >inlist)
    sections( i_output, >outlist)

    # define disposition of output...
    if( verbose)
        logfile = "STDOUT"
    else
        logfile = ""

    # check for existing solution. If found, prompt whether to use it...
    database = "database"
    solution =  database + "/fc" + i_ref
    if( access( solution))
        skipfit = useold

    if( ! skipfit){

        # set up cursor input file for ident...
        print( "y", >>cursor)
        print( "f", >>cursor)
        print( "q", >>cursor)
        print( "q", >>cursor)

        # identify slit edges...
        identify( i_ref, section="middle line", database=database,
            coordlist="", units="", nsum=nsum, match=-5., maxfeatures=2, 
            zwidth=100., ftype="emission", fwidth=25., cradius=25., 
            threshold=0., minsep=100., function="chebyshev", 
            order=2, sample="*", niterate=0, low_reject=3., high_reject=3., 
            grow=0., autowrite=yes, graphics="stdgraph", cursor=cursor)
        delete( cursor, ver-)

        # reidentify features across spectrum...
        reidentify( reference=i_ref, images=i_ref, 
            interactive-, section="middle line", newaps=yes, override=yes,
            refit=yes, trace=no, step=nstep, nsum=nsum, shift=0., 
            search=INDEF, nlost=0, cradius=15., threshold=0.,
            addfeatures-, match=-10., maxfeatures=2, minsep=6., 
            database=database, logfiles=logfile, 
            plotfile="", verbose=verbose, graphics="stdgraph")

        # derive the fit...
        fitcoords( images=i_ref, fitname="", interactive+, 
            combine-, database=database, function="chebyshev", 
            xorder=2, yorder=3, logfiles=logfile, plotfile="", 
            graphics="stdgraph", deletions="")
    }

    # transform science images...
    transform( input="@"//inlist, output="@"//outlist, fitnames=i_ref, 
        database=database, interptype="linear", 
        x1=INDEF, x2=INDEF, dx=INDEF, nx=INDEF, 
        xlog=no, y1=INDEF, y2=INDEF, dy=INDEF, ny=INDEF, 
        ylog=no, flux=yes, logfiles=logfile)
    
    # clean up...
    delete( inlist, ver-)
    delete( outlist, ver-)
    
end

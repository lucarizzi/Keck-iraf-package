#-----------------------------------------------------------------------
procedure lrisbias( images)
#-----------------------------------------------------------------------
# wmkolris.lrisbias
#
# Authors:
#   Gregory D. Wirth [wirth@keck.hawaii.edu]
#   Chris Fassnacht  [fassnacht@physics.ucdavis.edu] (minor modifications only)
#
# Purpose:
#     Subtract overscan from LRIS images based on a fit to the overscan region, 
#   plus add header keywords useful for processing with CCDPROC.
# 
# Modification history:
#  05 Jun 1998	G. Wirth	Original version
#  24 Feb 1999	G. Wirth	Rewrote using CCDPROC calls for bias removal
#  20 Apr 1999	G. Wirth	Fixed Attempt to access undefined local 
#				variable `upp1' bug with numamps=1
#  06 Jun 2000	G. Wirth	Revised treatment of 2-amp mode to better deal
#				with the case of all data on RIGHT amp;
#				also fixed treatment of binned images
#  13 Jun 2000	G. Wirth	Changed access() to imaccess()
#  13 Jul 2001	G. Wirth	Specified 'output=""' on ccdproc() call
#  12 Feb 2002	G. Wirth	Removed 'biassec' from being added to header
#  16 Feb 2002	G. Wirth	Added "padline=30" to prevent "hedit"
#				task from overwriting data section of images
#  28 Apr 2004	C. Fassnacht	Added a numamps=4 option to deal with new
#				LRIS-blue setup
#  13 Nov 2007	G. Wirth	Added "addonly" and "del" params to all hedit calls
#-----------------------------------------------------------------------

string  images       { prompt="Images to correct"}
bool    interactive  { no, prompt="Fit overscan interactively?"}
string  function     { "legendre", prompt="Fitting function" }
int     order        { 1, prompt="Number of polynomial terms or spline pieces", min=1 }
string  sample       { "*", prompt="Sample points to fit" }
int     naverage     { 1, prompt="Number of sample points to combine" }
int     niterate     { 1, prompt="Number of rejection iterations", min=0 }
real    low_reject   { 3., prompt="Low sigma rejection factor", min=0. }
real    high_reject  { 3., prompt="High sigma rejection factor", min=0. }
real    grow         { 0., prompt="Rejection growing radius", min=0. }
real    rdnoise      { INDEF, prompt="Readout Noise [e-]"}
real    gain         { INDEF, prompt="Inverse CCD gain [e-/ADU]"}
int     dispaxis     { -1, prompt="Dispersion Axis (if > 0)"}
string  observat     { "keck", prompt="Observatory database key"}
bool    noproc       { no, prompt="List processing steps only?"}
string	version		 {"1.0"}

# file structure pointers for reading lists
struct  *ilist       { prompt="Just leave this blank"}

begin
     string     i_images               # internal copy of images
     #-----
     string     tmpprefix              # prefix for temp files
     string     ifile, ofile           # temporary files containing expanded lists
     string     in                     # loop variable for stepping thru file lists
     string     ccdsec                 # readout window expressed a la IRAF
     string     ccdsum                 # the on chip summation
     string     biasstr                # bias section (all amps)
     string     biassxn1               # bias section (amp 1)
     string     biassxn2               # bias section (amp 2)
     string     biassxn3               # bias section (amp 3)
     string     biassxn4               # bias section (amp 4)
     string     datastr                # data section (all amps)
     string     datastr1               # data section (amp 1)
     string     datastr2               # data section (amp 2)
     string     datastr3               # data section (amp 3)
     string     datastr4               # data section (amp 4)
     string     overstr                # value of the overscan card on output
     string     timestr                # label for output header cards
     string     s1,s2,s3,s4,s5,s6      # pieces of the date
     string     buf                    # string buffer
     string     temp                   # working image
     int        cdelt1, cdelt2         # binning or on chip summation
     int        crval1, crval2         # offset of readout window from chip origin
     int        prepix                 # number of prescan pixels per amplifier
     int        postpix                # number of postscan pixels per amplifier
     int        numamps                # number of amplifiers used
     int        winchip                # chip number to which subraster applies
     int        winxstart              # starting column of subraster [0 indexed]
     int        winystart              # starting row of subraster [0 indexed]
     int        winxlen                # subraster extent along columns
     int        winylen                # subraster extent along rows
     int        x1                     # first column with data [logical c.s.]
     int        x2                     # last column with data (1 or 2 amp) [logical c.s.]
     int        ncols=2048             # number of CCD columns
     int        x_mid                  # last column in amp 1 readout
     int        x_mid2                 # last column in amp 2 readout
     int        x_mid3                 # last column in amp 3 readout
     int        y1                     # first line with data             
     int        y2                     # last line with data              
     int        l1                     # first line of bias section       
     int        l2                     # last line of bias section        
     int        bx11                   # first col of bias section 1
     int        bx12                   # last col of bias section 1
     int        bx21                   # first col of bias section 1
     int        bx22                   # last col of bias section 1
     int        bx31                   # first col of bias section 1
     int        bx32                   # last col of bias section 1
     int        bx41                   # first col of bias section 1
     int        bx42                   # last col of bias section 1
     int        nxlow                  # last column of lower-half of CCD 
     int        x11                    # first column in amp 1 data section
     int        x12                    # last column in amp 1 data section
     int        x21                    # first column in amp 2 data section
     int        x22                    # last column in amp 2 data section
     int        x31                    # first column in amp 3 data section
     int        x32                    # last column in amp 3 data section
     int        x41                    # first column in amp 4 data section
     int        x42                    # last column in amp 4 data section
     bool       debug=no               # print out diagnostic info?

     # define temp files...
     tmpprefix = "lrisbias"
     ifile = mktemp ( tmpprefix)
     temp = mktemp( tmpprefix) + ".fits"

     # query...
     i_images = images

     # expand image template...
     sections( i_images, >ifile)
	 if (sections.nimages == 0)
		error( 1, "no images in input list")

     # set file pointer to look at our input list, and loop...
     ilist = ifile
     while( fscan( ilist, in) != EOF && imaccess( in)) {

          # be sure that the baseline correction has not been done before now
          s1 = ""
          hselect( in, "overscan", yes) | scan( s1)
          if( s1 != "") {
              print( in//":\n\tAlready baselined --- no action taken.")
              next
          }

          # parse the header cards to get the number of bias samples 
          # both before and after the data...
          hselect( in, "prepix,postpix,numamps", "yes") \
          | scan( prepix, postpix, numamps)

          # legal values for numamps are 1, 2, and 4...
          if( numamps!=1 && numamps!=2 && numamps!=4)
               error( 1, "illegal value for NUMAMPS in image "//in)

          # parse the WINDOW header card to obtain image size and offset...
          hselect( in, "window", "yes") \
          | translit( "STDIN", ",", " ", delete-, collapse-) \
          | scan( winchip, winxstart, winystart, winxlen, winylen)
          if( nscan()<5){
               printf( "  Warning: image %s skipped --- can't read WINDOW header card\n", in)
               next
          }

          # Add ugly kludge here to fix winxlen for LRIS blue 4-amp case
          #  At least for data taken in April 2004, the WINDOW header card for 
          #   LRIS-B images is set to 1,0,0,2048,4096.  In fact, the window 
          #   x size appears to be 4096 and not 2048. (CDF)
          if (numamps==4) {
             winxlen = 4096
          }

          # read the header cards to get the binning (or on chip summation)...
          hselect( in, "binning", "yes") \
          | translit( "STDIN", "\",", " ", delete-, collapse-) \
          | scan( cdelt1, cdelt2)

          # compute the window offset from the CCD origin...
          crval1 = winxstart
          crval2 = winystart

          #--------------------------------------------------
          # Determine parameters for ccdproc...
          #--------------------------------------------------

          # define the data section (logical coordinates)...
          x1 = 1 + prepix * numamps
          x2 = x1 + winxlen - 1
          y1 = 1
          y2 = winylen
          printf( "[%d:%d,%d:%d]\n", x1, x2, y1, y2) | scan( datastr)

          # define the ROWS to sample for bias determination...
          l1 = y1
          l2 = y2

          # define columns for which bias is taken from first amp...
          bx11 = x2 + 1
          bx12 = x2 + postpix
          printf( "[%d:%d,%d:%d]\n", bx11, bx12, l1, l2) | scan( biassxn1)

          # define columns for which bias is taken from second amp...
          if( numamps>1){
               bx21 = x2 + postpix + 1
               bx22 = x2 + postpix + postpix
               printf( "[%d:%d,%d:%d]\n", bx21, bx22, l1, l2) | scan( biassxn2)
          }

          # define columns for which bias is taken from third and fourth amp...
          if( numamps==4){
               bx31 = x2 + postpix + postpix + 1
               bx32 = x2 + postpix + postpix + postpix
               printf( "[%d:%d,%d:%d]\n", bx31, bx32, l1, l2) | scan( biassxn3)
               bx41 = x2 + postpix + postpix + postpix + 1
               bx42 = x2 + postpix + postpix + postpix + postpix
               printf( "[%d:%d,%d:%d]\n", bx41, bx42, l1, l2) | scan( biassxn4)
          }

          # define the data sections...
          if( numamps==1) {

               datastr1 = datastr

          } else {

               # determine the columns at which the amps change,
               # accounting for on-chip binning via 
               # use of CDELT1...
               # First handle the first two amps.
               x_mid = x1 + ncols/2/cdelt1 - winxstart - 1
               x_mid2 = x1 + ncols/cdelt1 - winxstart - 1
               if( numamps == 4) {
                  x_mid3 = x1 + 3*ncols/2/cdelt1 - winxstart - 1
               }
               if( debug){ printf("x1=%d winxstart=%d x_mid=%d\n", 
                    x1, winxstart, x_mid) }

               # case 1: readout window lies on both amps...
               if( x1<=x_mid && x2>x_mid) {

                    x11 = x1
                    x12 = x_mid
                    printf( "[%d:%d,%d:%d]\n", x11, x12, y1, y2) \
                         | scan( datastr1)

                    x21 = x_mid+1
                    x22 = x_mid2
                    printf( "[%d:%d,%d:%d]\n", x21, x22, y1, y2) \
                         | scan( datastr2)

               # case 2: readout window is entirely on one amp...
               } else {

                    numamps = 1
                    printf( "[%d:%d,%d:%d]\n", x1, x2, y1, y2) | scan( datastr1)

                    # if using only RIGHT amp, then use the correct bias sxn...
                    if( x1>x_mid){ biassxn1 = biassxn2 }

               }
               # Now add section for 4-amp case
               if( numamps == 4) {
                 x31 = x_mid2+1
                 x32 = x_mid3
                 printf( "[%d:%d,%d:%d]\n", x31, x32, y1, y2) \
                         | scan( datastr3)
                 x41 = x_mid3+1
                 x42 = x2
                 printf( "[%d:%d,%d:%d]\n", x41, x42, y1, y2) \
                         | scan( datastr4)

               }
          }

          # WARNING:  It is not obvious what CCDSEC should be when binned.
          # We adopt here the assumption that a binned detector is an entirely
          # different detector which has fewer pixels. ( This seems to be most
          # consistent with the IRAF way of viewing the world.)
          # If in the Lick header(  CRVALn % CDELTn == 0 ), i.e. there
          # is an integral number of binned pixels preceding the readout
          # window, then this definition of CCDSEC makes sense.  If the
          # above is not true, then all bets are off because the result
          # must depend on which detector, how many amplifiers, and other
          # factors which have never been specified for Lick CCD readouts.
           printf( "[%d:%d,%d:%d]\n",
               1+(crval1/cdelt1),
               0+winxlen+(crval1/cdelt1),
                1+(crval2/cdelt2),
               0+winylen+(crval2/cdelt2)) \
          | scan( ccdsec)

          # Define binning or on chip summation parameters...
          # [Note that as far as CCDPROC is concerned, a chip using binning
          #  is a different "instrument" than the same chip not using binning]
          ccdsum=str(cdelt1)+' '+str(cdelt2)

          #--------------------------------------------------
          # Give diagnostic output...
          #--------------------------------------------------

          if( noproc){
               print( in, ":")
               if( numamps==1){
                    print( "  [TO BE DONE] Overscan section is ", biassxn1)
               } else {
                    print( "  [TO BE DONE] Overscan section #1 is ", biassxn1)
                    print( "  [TO BE DONE] Overscan section #2 is ", biassxn2)
                    if( numamps==4) {
                       print( "  [TO BE DONE] Overscan section #3 is ", biassxn3)
                       print( "  [TO BE DONE] Overscan section #4 is ", biassxn4)
                    }
               }
               next
          }

          #--------------------------------------------------
          # Invoke ccdproc for bias removal...
          #--------------------------------------------------

          if( debug) { printf( "datasec1=%s\n", datastr1) }
          hedit( in, "ccdsec,datasec", add-, delete+, verify-, show-, update+, addonly-)
          hedit( in, "datasec", datastr1, add+, delete-, verify-, show-, update+, addonly-)
          ccdproc( in, output="", ccdtype="", noproc=noproc, fixpix-,
               overscan+, trim-,
               zerocor-, darkcor-, flatcor-, illumcor-, fringecor-, readcor-, 
               scancor-, readaxis="line", biassec=biassxn1,
               interactive=interactive, function=function, order=order, 
               sample=sample, naverage=naverage, niterate=niterate, 
               low_reject=low_reject, high_reject=high_reject, grow=grow)

          if( numamps>1) {
               if( debug) { printf( "datasec2=%s\n", datastr2) }
               hedit( in, "ccdsec,datasec,overscan", add-, delete+, verify-, 
                    show-, update+, addonly-)
               hedit( in, "datasec", datastr2, add+, delete-, verify-, 
                    show-, update+, addonly-)
               ccdproc( in, output="", ccdtype="", noproc=noproc, fixpix-,
                    overscan+, trim-,
                    zerocor-, darkcor-, flatcor-, illumcor-, fringecor-, readcor-, 
                    scancor-, readaxis="line", biassec=biassxn2,
                    interactive=interactive, function=function, order=order, 
                    sample=sample, naverage=naverage, niterate=niterate, 
                    low_reject=low_reject, high_reject=high_reject, grow=grow)
          }

          if( numamps==4) {
               if( debug) { printf( "datasec3=%s\n", datastr3) }
               if( debug) { printf( "datasec4=%s\n", datastr4) }
               hedit( in, "ccdsec,datasec,overscan", add-, delete+, verify-, 
                    show-, update+, addonly-)
               hedit( in, "datasec", datastr3, add+, delete-, verify-, 
                    show-, update+, addonly-)
               ccdproc( in, output="", ccdtype="", noproc=noproc, fixpix-,
                    overscan+, trim-,
                    zerocor-, darkcor-, flatcor-, illumcor-, fringecor-, readcor-, 
                    scancor-, readaxis="line", biassec=biassxn3,
                    interactive=interactive, function=function, order=order, 
                    sample=sample, naverage=naverage, niterate=niterate, 
                    low_reject=low_reject, high_reject=high_reject, grow=grow)
               if( debug) { printf( "datasec4=%s\n", datastr4) }
               hedit( in, "ccdsec,datasec,overscan", add-, delete+, verify-, 
                    show-, update+, addonly-)
               hedit( in, "datasec", datastr4, add+, delete-, verify-, 
                    show-, update+, addonly-)
               ccdproc( in, output="", ccdtype="", noproc=noproc, fixpix-,
                    overscan+, trim-,
                    zerocor-, darkcor-, flatcor-, illumcor-, fringecor-, readcor-, 
                    scancor-, readaxis="line", biassec=biassxn4,
                    interactive=interactive, function=function, order=order, 
                    sample=sample, naverage=naverage, niterate=niterate, 
                    low_reject=low_reject, high_reject=high_reject, grow=grow)
          }

          #--------------------------------------------------
          # Modify the image header...
          #--------------------------------------------------

          # define a card to indicate who did the processing...
          date | scan(s1,s2,s3,s4,s5,s6)
          printf( "LRISBIAS processing completed %s %s %s %s %s %s\n", s1, s6, s2, s3, s4, s5) \
          | scan( line)
          timestr = line

          # set the OVERSCAN card....
          if( numamps==1){
               printf( "Overscan section is %s\n", biassxn1) \
               | scan( line)
               biasstr = biassxn1
          } else {
               if( numamps==2) {
                  printf( "Overscan sections are %s and %s\n", biassxn1, biassxn2) \
                  | scan( line)
                  biasstr = biassxn1//","//biassxn2
               } else {
                  printf( "Overscan sections are %s, %s, %s, and %s\n", biassxn1, biassxn2,biassxn3,biassxn4) \
                  | scan( line)
                  biasstr = biassxn1//","//biassxn2//","//biassxn3//","//biassxn4
               }
          }
          overstr = line

          #----------------------------------------
          # Create working version of input image
          # while allocating sufficient space to 
          # add header keywords
          #----------------------------------------
          imcopy( in, temp//"[padline=30]", verb-)

          # Add in relevant header cards for the sake of ccdproc.
          # Note that since the data are already overscan subtracted, 
          # we do NOT add the card for BIASSEC.
          # Doing so would cause CCDPROC to choke...
          hedit( temp, 'ccdsec'  , ccdsec  , add+, update+, show-, ver-, del-, addonly-)
          hedit( temp, 'ccdsum'  , ccdsum  , add+, update+, show-, ver-, del-, addonly-)
          hedit( temp, 'overscan', overstr , add+, update+, show-, ver-, del-, addonly-)
          hedit( temp, 'datasec' , datastr , add+, update+, show-, ver-, del-, addonly-)
          hedit( temp, 'trimsec' , datastr , add+, update+, show-, ver-, del-, addonly-)
          hedit( temp, 'lrisbias', timestr , add+, update+, show-, ver-, del-, addonly-)

          # define observatory and dispaxis keywords...
          if( observat != "" )
               hedit( temp, 'observat', observat, add+, update+, show-, ver-, del-, addonly-)
          if( dispaxis > 0 )
               hedit( temp, 'dispaxis', dispaxis, add+, update+, show-, ver-, del-, addonly-)

          # store gain and read noise in header...
          if( rdnoise >= 0. )
              hedit( temp, 'rdnoise' , rdnoise , add+, update+, show-, ver-, del-, addonly-)
          if(  gain >0. )
              hedit( temp, 'gain'    , gain    , add+, update+, show-, ver-, del-, addonly-)

          # replace input with new image...
          imdelete( in, verify-)
          imrename( temp, in, verb-)

    }

    # clear struct...
    ilist = ""

    # remove temp file...
    delete( ifile, ver-)

end

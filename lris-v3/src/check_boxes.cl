#-----------------------------------------------------------------------
procedure check_boxes( image, boxfile)
#-----------------------------------------------------------------------
# Name:
#   check_boxes
# 
# Purpose:
#   define and verify alignment box coordinates
# 
# Author:
#   Gregory D. Wirth, W. M. Keck Observatory
# 
# Notes:
#   Coordinate systems --- all internal coordinates are maintained in
#   the coordinate system of the image.  When coordinates are written back
#   to the lbox input file they are put back into the lbox (physical CCD)
#   coordinate system.
#
# Modification history:
#   2001-Jun-02     GDW     original version
#   2001-Sep-09     GDW     Accounted for NUMAMPS and PREPIX as in lbox
#   2002-Jun-01     AES     Calls lbox, displays image "fill+"
#   2008-May-07     GDW     add "minmax" option, remove "fill+"
#   2009-Apr-18     GDW     add default values for image and boxfile
#	2009-Jun-02		GDW		Add support for multi-HDU images
#	2009-Jun-21		GDW		Ad printout of boxfile name
#-----------------------------------------------------------------------

string  image="last"   	{ prompt="Name of input image ('last' for most recent)" }
string  boxfile="default" { prompt="Name of input file listing box centers" }
bool    dolbox  		{ yes, prompt="Launch lbox when done?" }
bool	minmax			{ no, prompt="Scale display using min/max?" }
struct  *ilist
struct  iline

begin
    string  tmpprefix   # prefix for temp files
    string  i_image     # internal copy of image
    string  i_boxfile	# internal copy of boxfile
    string  buf         # string buffer
    string  file1       # coords file (PHYSICAL coord system)
    string  file2       # new file
    string  name0       # designation for current box
    string  slitname	# name of slitmask
    string  dir			# directory with image
    char    key         # key pressed by user
    real    d2          # squared distance between two points
    real    d2_min      # minimum value of d2 in sample
    real    minval=INDEF      # minimum value in image display
    real    maxval=INDEF      # minimum value in image display
	real	pscale=0.135	# pixel scale
    int     x,y         # cursor position
    int     dx,dy       # shift in rows and columns
    int     x0,y0       # column and row coordinates of box
    int     x1,y1       # column and row coordinates of box
    int     xtmp,ytmp   # coordinates read from file
    int     xoff,yoff   # lbox offsets
    int     frame=1     # frame in which to display
    int     i           # counter
    int     i_min       # counter of target at minimum distance d2
    int     n           # number of items scanned
    int     box_size    # box length for tvmark
    int     offset      # x/y offset for tvmark labels
	int		nhdu		# number of header data units in input image
    bool    display=yes # flag whether to redisplay image
    bool    mark=yes    # flag whether to remark image
    bool    change=no   # flag whether list of boxes has changed since original
    bool    swap=no     # flag whether list of boxes has changed this iteration
    bool    do_label=no # flag whether to label boxes with names
    bool    do_number=yes   # flag whether to label boxes with numbers

    # translate last image...
    i_image = image
    if( i_image == "last" ) {
		print("scanning for last image...")
		lastimageb() | scan( i_image)
    }

    # verify image...
    if ( ! imaccess( i_image))
        error( 1, "requested input image "//i_image//" not found")

    # generate boxfile name as needed...
    i_boxfile = boxfile
    if( i_boxfile == "default" ) {

        hselect( i_image//"[0]", "slitname", "yes") | scan( slitname)
        if ( slitname == "direct" ){
            print( "I can only determine the box locations if you use a slitmask!")
            error( 1, "")
        }

        i_boxfile = slitname + '.box'
    }
    print( "boxfile is ", i_boxfile)

    # define temp files...
    tmpprefix = "CHECK_BOXES"
    file1 = mktemp( tmpprefix)
    file2 = mktemp( tmpprefix)

    # define params for tvmark...
    box_size = min(lbox.xsz,lbox.ysz)/pscale
    offset = 3+box_size/2

    # make local copy of box list...
    if( access( i_boxfile)){

        printf( "Input box list found; got the following boxes:\n\n")
        printf( "  %1s %12s   %11s   %11s\n", "", "", "Physical", "Screen")
        printf( "  %1s %12s   %11s   %11s\n", "", "", "Coordinates", "Coordinates")
        printf( "  %1s %-12s   %5s %5s   %5s %5s\n", "N", "Name", "x", "y", "x", "y")
        printf( "  %1s %-12s   %5s %5s   %5s %5s\n", 
            "-", "------------", "-----", "-----", "-----", "-----")

        # transform coordinates...
        ilist = i_boxfile
        i = 0
        while( fscan( ilist, iline)!=EOF){

            # default name is blank...
            name0 = ""

            # attempt to parse line...
            n = fscan( iline, x0, y0, name0)
    
            # if x0,y0 were not read, assume the line is a 
            # comment and pass it to the output file...
            if( n<2){
                print( iline, >>file1)
                next
            }

            # if name was scanned, use names for labels...
            if ( n>2){
                do_label = yes
                do_number = no
            }

            # print input coords...
            i += 1
            printf( "  %d %-12s   %5d %5d", i, name0, x0, y0)

#            # correct coordinates for the offset...
#            x0 += xoff
#            y0 += yoff

            # print output coords...
            printf( "   %5d %5d\n", x0, y0)

            # send corrected coordinates to output file...
            print( x0, y0, name0, >>file1)
        }
        printf( "\n")

    } else {

        # generate empty file...
        printf( "No input box list found --- defining new file.\n")
        type( "", >file1)

    }

#    # if using minmax scaling, derive scaling limits...
#    if ( minmax ) {
#		printf( "Computing image min/max values.\n")
#		minmax( i_image, force+, update-, verb-)
#		minval = minmax.minval
#		maxval = minmax.maxval
#    }

    # start loop...
    while( yes ){

        # display image...
        if( display){
            printf( "Displaying image %s.\n", i_image)
	    if ( minmax )
           	lrisdisplay( i_image, frame, zsc-, zra+, ztr="log")
#           	lrisdisplay( i_image, frame, zsc-, zra+, z1=minval, z2=maxval, 
#				fill-, ztr="log")
	    else
                lrisdisplay( i_image, frame, zsc+, zra-)
        }

        # mark boxes...
        if( mark){
            printf( "Marking box locations on image.\n")
            tvmark( frame=1, coords=file1, logfile="", autolog-,
                outimage="", deletions="", commands="", mark="circle",
                lengths=str(box_size), radii=str(nint(0.5*box_size+0.5)), 
                font="raster", color=204, 
                label=do_label, number=do_number,
                nxoffset=offset, nyoffset=offset, pointsize=2,
                txsize=6, tolerance=1.5, interactive-)
        }

        # define actions...
        display = no
        mark = no
        swap = no

        # print instructions...
        printf( "\nCommand options:\n" )
        printf( "  ")
        printf( "%-25s", "m = mark a new box" )
        printf( "%-25s", "q = quit and save" )
        printf( "%-25s", "I = abort (no save)" )
        printf( "\n")
        printf( "  ")
        printf( "%-25s", "d = delete a box" )
        printf( "%-25s", "D = delete all boxes" )
        printf( "%-25s", "- = shift all by -21 px" )
        printf( "\n")
        printf( "  ")
        printf( "%-25s", "s = shift one box" )
        printf( "%-25s", "S = shift all boxes" )
        printf( "%-25s", "+ = shift all by +21 px" )
        printf( "\n")
        printf( "Please place cursor on image display and press a key: ")

        # get user input...
        n = fscan( imcur, x, y, buf, key)
		if ( nscan() != 4){
			print("WARNING: invalid keystroke!")
			next
		}
		print( key)

        # parse input...
        if( key == "m"){

            #----------------------------------------
            # m = mark new box
            #----------------------------------------

            # add new box...
            printf( "Adding new box at image coordinates %d,%d.\n",
                x, y)
            printf( "%d %d\n", x, y, >>file1)
            change = yes
            display = no
            mark = yes

        } else if( key == "d") {

            #----------------------------------------
            # d = delete box
            #----------------------------------------

            # locate closest box to cursor...
            ilist = file1
            d2_min = 1.e32
            i_min = INDEF
            i = 0
            while( fscan( ilist, x0, y0)!=EOF){
                n = nscan()
                if( n==2){
                    i += 1
                    d2 = (x-x0)**2 + (y-y0)**2
                    if( d2 < d2_min){
                        d2_min = d2
                        i_min = i
                    }
                }
            }

            # check to ensure that a box was found...
            if( i_min == INDEF){
                print( "ERROR -- no boxes available to delete!")
                beep()
                next
            }

            # output the list minus the offending box...
            ilist = file1
            i = 0
            while( fscan( ilist, iline)!=EOF){
                n = fscan( iline, x0, y0)
                if( n==2){ i += 1 }
                if( n==2 && i==i_min){
                    printf( "Deleting box at image coordinates %d,%d.\n",
                        x0, y0)
                } else {
                    print( iline, >>file2)
                }
            }

            swap = yes
            display = yes
            mark = yes

        } else if( key == "D") {

            #----------------------------------------
            # D = delete all boxes
            #----------------------------------------

            # output the list minus all boxes...
            ilist = file1
            i = 0
            while( fscan( ilist, iline)!=EOF){
                n = fscan( iline, x0, y0)
                if( n==2){
                    printf( "Deleting box at image coordinates %d,%d.\n",
                        x0, y0)
                } else {
                    print( iline, >>file2)
                }
            }

            swap = yes
            display = yes
            mark = yes

        } else if ( key == "s") {

            #----------------------------------------
            # s = shift one box
            #----------------------------------------

            # get starting and ending positions...
            print( "Preparing to shift ONE box." )
            print( "Step 1: center cursor on box to be shifted and press key")
            n = fscan( imcur, x, y, buf, key)
            if( n < 4){
                printf( "ERROR reading cursor input --- action aborted")
                beep()
                next
            }

            # locate closest box to cursor...
            ilist = file1
            d2_min = 1.e32
            i = 0
            i_min = INDEF
            while( fscan( ilist, xtmp, ytmp)!=EOF){
                n = nscan()
                if( n==2){
                    i += 1
                    d2 = (x-xtmp)**2 + (y-ytmp)**2
                    if( d2 < d2_min){
                        d2_min = d2
                        i_min = i
                        x0 = xtmp
                        y0 = ytmp
                    }
                }
            }
            if( d2_min == INDEF){
                printf( "ERROR finding nearest box --- action aborted")
                beep()
                next
            }

            print( "Step 2: center cursor on new box position and press key")
            n = fscan( imcur, x1, y1, buf, key)
            if( n < 4)
                next
            dx = x1 - x0
            dy = y1 - y0
            printf( "Shifting box %d by %d px in x and %d px in y...\n",
                i_min, dx, dy)

            # make changes to coordinates of the i-th box...
            ilist = file1
            i = 0
            while( fscan( ilist, iline)!=EOF){
                name0 = ""
                n = fscan( iline, x0, y0, name0)
                if( n<2){
                    print( iline, >>file2)
                    next
                }
                i += 1
                if( i==i_min){
                    x0 += dx
                    y0 += dy
                }
                print( x0, y0, name0, >>file2)
            }
            swap = yes
            change = yes
            display = yes
            mark = yes

        } else if ( key == "S") {

            #----------------------------------------
            # S = shift all boxes
            #----------------------------------------

            # get starting and ending positions...
            printf( "Preparing to shift ALL boxes.\n" )
            printf( "Step 1: center cursor on any box and press key: ")
            n = fscan( imcur, x, y, buf, key)
            if( n < 4){
                printf( "ERROR reading cursor input --- action aborted")
                beep()
                next
            }

            # locate closest box to cursor...
            ilist = file1
            d2_min = 1.e32
            i = 0
            i_min = INDEF
            while( fscan( ilist, xtmp, ytmp)!=EOF){
                n = nscan()
                if( n==2){
                    i += 1
                    d2 = (x-xtmp)**2 + (y-ytmp)**2
                    if( d2 < d2_min){
                        d2_min = d2
                        i_min = i
                        x0 = xtmp
                        y0 = ytmp
                    }
                }
            }

            if( x0 == INDEF){
                printf( "\nERROR finding nearest box --- action aborted")
                beep()
                next
            }

            printf( "got box at %d,%d.\n", x0, y0)

            printf( "Step 2: center cursor on new box position and press key: ")
            n = fscan( imcur, x1, y1, buf, key)
            if( n < 4)
                next
            dx = x1 - x0
            dy = y1 - y0
            printf( "done.\nShifting ALL boxes by %d px in x and %d px in y.\n",
                dx, dy)

            # make changes to coordinates of ALL boxes...
            ilist = file1
            while( fscan( ilist, iline)!=EOF){

                name0 = ""
                n = fscan( iline, x0, y0, name0)
                if( n<2){
                    print( iline, >>file2)
                    next
                }
                x0 += dx
                y0 += dy
                print( x0, y0, name0, >>file2)
            }
            swap = yes
            change = yes
            display = yes
            mark = yes

        } else if ( key == "-") {

            #----------------------------------------
            # - = shift all by -21 px
            #----------------------------------------

            dx = -21
            dy = 0
            printf( "Shifting ALL boxes by %d px in x and %d px in y.\n",
                dx, dy)

            # make changes to coordinates of ALL boxes...
            ilist = file1
            while( fscan( ilist, iline)!=EOF){

                name0 = ""
                n = fscan( iline, x0, y0, name0)
                if( n<2){
                    print( iline, >>file2)
                    next
                }
                x0 += dx
                y0 += dy
                print( x0, y0, name0, >>file2)
            }
            swap = yes
            change = yes
            display = yes
            mark = yes

        } else if ( key == "+") {

            #----------------------------------------
            # + = shift all by +21 px
            #----------------------------------------

            dx = 21
            dy = 0
            printf( "Shifting ALL boxes by %d px in x and %d px in y.\n",
                dx, dy)

            # make changes to coordinates of ALL boxes...
            ilist = file1
            while( fscan( ilist, iline)!=EOF){

                name0 = ""
                n = fscan( iline, x0, y0, name0)
                if( n<2){
                    print( iline, >>file2)
                    next
                }
                x0 += dx
                y0 += dy
                print( x0, y0, name0, >>file2)
            }
            swap = yes
            change = yes
            display = yes
            mark = yes

        } else if ( key == "q" ) {

            #----------------------------------------
            # q = quit and save
            #----------------------------------------

            # replace old file with new...
            if( change){

                # transform coordinates...
                ilist = file1
                while( fscan( ilist, iline)!=EOF){
                    name0 = ""
                    n = fscan( iline, x0, y0, name0)
                    if( n<2){
                        print( iline, >>file2)
                        next
                    }
#                    x0 -= xoff
#                    y0 -= yoff
                    print( x0, y0, name0, >>file2)
                }
                if( access( i_boxfile))
                    delete( i_boxfile, ver-)
                copy( file2, i_boxfile, ver-)
                printf( "Normal exit --- saved revised version of box file '%s'\n",
                    i_boxfile)
            } else {
                printf( "Normal exit --- no changes made to box file '%s'\n",
                    i_boxfile)
            }
            break

        } else if ( key == "I" ){

            #----------------------------------------
            # I = interrupt
            #----------------------------------------

            if( access( file1))
                delete( file1, ver-)
            if( access( file2))
                delete( file2, ver-)
            error( 1, "Aborting!")

        } else {

            printf( "ERROR: '%s' is not a valid choice.  Please try again.\n",
                key)
            beep()
            next

        }

        # swap files...
        if( swap){              
            delete( file1, ver-)
            file1 = file2
            file2 = mktemp( tmpprefix)

        }

        # ensure that file1 exists...
        if( ! access( file1))
            type( "", >file1)
    }

    # clean up...
    ilist = ""
    if( access( file1))
        delete( file1, ver-)
    if( access( file2))
        delete( file2, ver-)

    # run lbox using the default lbox parameters...
    if( dolbox)
        lbox( i_image, i_boxfile, xoff=0, yoff=0, practice-)
end

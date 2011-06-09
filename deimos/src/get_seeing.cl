#-----------------------------------------------------------------------
procedure get_seeing( image)
#-----------------------------------------------------------------------
string	image="last"	{ prompt="Name of DEIMOS alignment image ('last' for most recent)" }
string	boxfile=""		{ prompt="File listing box positions (blank=>auto)" }
string	logfile="get_seeing.log"	{ prompt="Name of logfile" }
bool	check=no		{ prompt="Run check_boxes for verification?" }
struct	*ilist 

begin
	string	i_image
    string  box_file    # box file
    string  slmsknam	# name of slitmask in use
	string	cursorfile	# file containing cursor input
	string	plotfile	# file containing graphics output input
	string	datafile	# file containing data
	string	parfile		# file containing parameters
	string	tmpprefix	# prefix for temp files
	string	host="deimos"
	string	subimage	# image extension to display
	string	ampmode
	real	xc,yc		# centroids
	real	fwhm_px
	real	fwhm_asec
	real	crval2
	real		ellipticity
	int		xsize
	int		frame=1
	int		n
	int		m
	int		k
	int		amp[8]
	int		flip[8]
	int		istat
	int		nx,ny
	bool	debug=no

	# prompt...
	i_image = image

	# declare temp files...
	tmpprefix = "tmp$get_seeing"
	cursorfile = mktemp( tmpprefix)
	plotfile = mktemp( tmpprefix)
	datafile = mktemp( tmpprefix)
	parfile = mktemp( tmpprefix)

	# translate last image...
	if( i_image == "last" ) {
		lastimage() | scan( i_image)
		i_image = "/s" // i_image
	}

	# verify that the image exists...
	if( ! imaccess( i_image))
		error( 1, "Image "//i_image//" not found")
	print( "Processing image ", i_image)

	if ( boxfile == "" ) {
		# read mask name from header...
		hselect( i_image//"[0]", "SLMSKNAM", "yes") | scan( slmsknam)
		if ( nscan() != 1 )
		   error( 1, "Can't read SLMSKNAM from image header")

		box_file = "box." // slmsknam
		if( ! access( box_file))
		    get_boxes( i_image)
		    
	} else {
		if( access( boxfile)){
			box_file = boxfile
		} else {
			error( 1, "Box file '"//boxfile// "' not found -- abort!")
		}	
	}

	# determine readout mode...
	hselect( i_image//"[0]", "AMPMODE", "yes") | scan( ampmode)

	# configure amp array...
	if( ampmode == "SINGLE:A" ){
		amp[1] = 1
		amp[2] = 2
		amp[3] = 3
		amp[4] = 4

		flip[1] = 1
		flip[2] = 1
		flip[3] = 1
		flip[4] = 1

		xsize = 2048
	} else if ( ampmode == "SINGLE:B" ){
		amp[1] = 1
		amp[2] = 2
		amp[3] = 3
		amp[4] = 4

		flip[1] = 0
		flip[2] = 0
		flip[3] = 0
		flip[4] = 0

		xsize = 2048
	} else if ( ampmode == "DUAL:A+B" ){
		amp[1] = 2
		amp[2] = 1
		amp[3] = 4
		amp[4] = 3
		amp[5] = 6
		amp[6] = 5
		amp[7] = 8
		amp[8] = 7

		flip[1] = 1
		flip[2] = 0
		flip[3] = 1
		flip[4] = 0
		flip[5] = 1
		flip[6] = 0
		flip[7] = 1
		flip[8] = 0

		xsize = 1024
	} else {
	  error( 1, "unsupported ampmode")
	}

	# execute check_boxes for alignment...
	if ( check)
		check_boxes( i_image, box_file, doxbox-)

	# store params...
	dpar( "rimexam", >parfile)
	rimexam.banner = yes
	rimexam.title = ""
	rimexam.xlabel = "Radius"
	rimexam.ylabel = "Pixel Value"
	rimexam.fitplot = yes
	rimexam.fittype = "gaussian"
	rimexam.center = yes
	rimexam.background = yes
	rimexam.radius = 12.
	rimexam.buffer = 0.
	rimexam.width = 3.
	rimexam.iterations = 1
	rimexam.xorder = 0
	rimexam.yorder = 0
	rimexam.magzero = 25.
	rimexam.beta = INDEF
	rimexam.rplot = 15.
	rimexam.x1 = INDEF
	rimexam.x2 = INDEF
	rimexam.y1 = INDEF
	rimexam.y2 = INDEF
	rimexam.pointmode = yes
	rimexam.marker = "plus"
	rimexam.szmarker = 1.
	rimexam.logx = no
	rimexam.logy = no
	rimexam.box = yes
	rimexam.ticklabels = yes
	rimexam.majrx = 5
	rimexam.minrx = 5
	rimexam.majry = 5
	rimexam.minry = 5
	rimexam.round = no

	# parse the file...
	ilist = box_file
	n = 0
	while( fscan( ilist, xc, yc) != EOF){
		n = n + 1

		# convert position in whole image to position in subimage...
		m = 1 + int(xc/xsize)
		k = amp[m]
		xc = mod(xc,xsize)

		# correct y position for CCD windowing...
		subimage = i_image//"["//str(k)//"]"
		imgets( subimage, "CRVAL2")
		istat = fscan( imgets.value, crval2)
		yc = yc - crval2

		# correct x position for image flip...
		if ( flip[k] == 1 ) {
			xc = xsize - xc
		}
		xc += 16

		print( xc, yc, " 100 r", >cursorfile)

		# examine file...
		if( debug )
			printf( "Measuring box at (%d,%d) in image %s\n", xc, yc, subimage)
		display( subimage, frame, zsc-, zra-, z1=1000, z2=2500, erase-, >"dev$null")

		if( debug ){
		    tvmark( frame=frame, coords=cursorfile, logfile="", autolog-,
			outimage="", deletions="", commands="", mark="circle",
			lengths="10", radii="10", 
			font="raster", color=204, 
			label-, number+,
			nxoffset=0, nyoffset=0, pointsize=2,
			txsize=6, tolerance=1.5, interactive-)
		    imexamine()
		}

		## printf( "running imexamine...\n")
		## type( cursorfile)
		imexamine( input="", frame=frame, image="", defkey="r",
			imagecur=cursorfile, keeplog-, >>Gplotfile, >>datafile)
		delete( cursorfile, ver-)
	}


	# generate plot...
	if ( n <= 4 ) {
		nx = 2
		ny = 2
	} else if ( n <= 9 ) {
		nx = 3
		ny = 3
	} else {
		nx = 4
		ny = 4
	}
	gkiextract( plotfile, "1-99") \
		| gkimosaic( device="stdgraph", output="", nx=nx, ny=ny, 
			rotate=no, fill=no, interactive=no)

	# print results...
	printf( "\n Star    FWHM    FWHM  Ellip   Flux\n")
	printf(   "         [px]   [asec]         [DN]\n")
	printf(   "------ ------- ------- ----- --------\n")
	tcalc( datafile, "c10", "c9*0.119")
	tprint( datafile, col="c9,c10,c6,c3", showhdr-)
	tstat( datafile, "c6", outtab="STDOUT", >"dev$null")
	ellipticity = tstat.median
	tstat( datafile, "c9", outtab="STDOUT", >"dev$null")
	fwhm_px = tstat.median
	tstat( datafile, "c10", outtab="STDOUT", >"dev$null")
	fwhm_asec = tstat.median
	printf(   "------ ------- ------- ----- --------\n")
	printf( "%-6s %7.2f %7.2f %5.2f\n", "Median", fwhm_px, fwhm_asec, 
		ellipticity)

	# print to logfile...
	if( logfile != "")
		printf( "%-40s FWHM = %5.2f px = %5.2f asec  Ellipticity = %5.2f\n", 
			i_image, fwhm_px, fwhm_asec, ellipticity, >>logfile)

	# clean up...
	cl( <parfile)
	if( access( cursorfile))
		delete( cursorfile, ver-)
	if( access( plotfile))
		delete( plotfile, ver-)
	if( access( datafile))
		delete( datafile, ver-)
	if( access( parfile))
		delete( parfile, ver-)
	ilist=""
	
end

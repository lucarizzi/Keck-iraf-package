#-----------------------------------------------------------------------
procedure seeing( image)
#-----------------------------------------------------------------------
string	image="last"	{ prompt="Name of DEIMOS alignment image ('last' for most recent)" }
string	logfile="get_seeing.log"	{ prompt="Name of logfile" }
bool	check=no	{ prompt="Run check_boxes for verification?" }
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
	real	xc,yc		# centroids
	real	fwhm_px
	real	fwhm_asec
	int		frame=1
	int		n

	# prompt...
	i_image = image

	# declare temp files...
	tmpprefix = "tmp$get_seeing"
	cursorfile = mktemp( tmpprefix)
	plotfile = mktemp( tmpprefix)
	datafile = mktemp( tmpprefix)
	parfile = mktemp( tmpprefix)

	# translate last image...
	if( i_image == "last" )
		lastimage("-tail") | scan( i_image)

	# verify that the image exists...
	if( ! imaccess( i_image))
		error( 1, "Image "//i_image//" not found")
	print( "Processing image ", i_image)

	# read mask name from header...
    hselect( i_image//"[0]", "SLMSKNAM", "yes") | scan( slmsknam)
    if ( nscan() != 1 )
		error( 1, "Can't read SLMSKNAM from image header")

	box_file = "box." // slmsknam
	if( ! access( box_file))
	    get_boxes( i_image)

	# if box file not found, try using the most recent one...
	if( ! access( box_file)){
			ls("-1t","box.*") | scan( box_file)
		if( access( box_file)){
			beep()
			printf( "WARNING: using alternate box file %s\n", box_file)
		} else {
		    error( 1, "Could not generate box file "//box_file)
		}
	}

	# execute check_boxes for alignment...
	if ( check)
		check_boxes( i_image, box_file, doxbox-)
	else 
		dmosdisplay( i_image, frame, zscale-, zra-, z1=1000, z2=10000)

	# parse the file...
	ilist = box_file
	n = 0
	while( fscan( ilist, xc, yc) != EOF){
		n = n + 1
		if ( n == 1 )
			print( xc, yc, "a", >>cursorfile)
		print( xc, yc, "r", >>cursorfile)
	}

	# store params...
	dpar( "rimexam2", >parfile)
	rimexam2.banner = yes
	rimexam2.title = ""
	rimexam2.xlabel = "Radius"
	rimexam2.ylabel = "Pixel Value"
	rimexam2.fitplot = yes
	rimexam2.fittype = "gaussian"
	rimexam2.center = yes
	rimexam2.background = yes
	rimexam2.radius = 5.
	rimexam2.buffer = 5.
	rimexam2.width = 5.
	rimexam2.iterations = 1
	rimexam2.xorder = 0
	rimexam2.yorder = 0
	rimexam2.magzero = 25.
	rimexam2.beta = INDEF
	rimexam2.rplot = 8.
	rimexam2.x1 = INDEF
	rimexam2.x2 = INDEF
	rimexam2.y1 = INDEF
	rimexam2.y2 = INDEF
	rimexam2.pointmode = yes
	rimexam2.marker = "plus"
	rimexam2.szmarker = 1.
	rimexam2.logx = no
	rimexam2.logy = no
	rimexam2.box = yes
	rimexam2.ticklabels = yes
	rimexam2.majrx = 5
	rimexam2.minrx = 5
	rimexam2.majry = 5
	rimexam2.minry = 5
	rimexam2.round = no

	# examine file...
	mscexamine( frame=frame, defkey="r", imagecur=cursorfile,
		>Gplotfile, >datafile)

	# generate plot...
	gkiextract( plotfile, "1-99") \
		| gkimosaic( device="stdgraph", output="", nx=3, ny=2, 
			rotate=no, fill=yes, interactive=no)

	# print results...
	printf( "\nStar FHMW [px] FWHM [asec]\n")
	printf( "\n---- --------- -----------\n")
	tcalc( datafile, "c10", "c9*0.119")
	tprint( datafile, col="c9,c10", showhdr-)
	tstat( datafile, "c9", outtab="STDOUT", >"dev$null")
	fwhm_px = tstat.median
	tstat( datafile, "c10", outtab="STDOUT", >"dev$null")
	fwhm_asec = tstat.median
	printf( "%-8s %5.2f %5.2f\n", "Median", fwhm_px, fwhm_asec)

	# print to logfile...
	if( logfile != "")
		printf( "%-40s FWHM = %5.2f px = %5.2f asec\n", 
			i_image, fwhm_px, fwhm_asec, >>logfile)

	# clean up...
	cl( <parfile)
	delete( cursorfile, ver-)
	delete( plotfile, ver-)
	delete( datafile, ver-)
	delete( parfile, ver-)
	ilist=""
	
end

#-----------------------------------------------------------------------
procedure get_seeing( image)
#-----------------------------------------------------------------------
# Name:
#	get_seeing
#
# Purpose:
#	Measure the seeing in an LRIS slitmask alignment image by fitting
#	Gaussian profiles to the stars in the alignment boxes.
#
# Usage:
#	get_seeing image boxfile
# 
# Parameters:
#
#	image = name of alignment image to process.  If the value is the
#		special string "last", then the program will obtain the 
#		name of the last LRIS blue-side image acquired.
#	boxfile = name of file listing alignment box coordinates 
#		(same format as for xbox)
#	logfile = name of a logfile in which to record results
#
# Restrictions:
#	- only works with LRIS blue-side images
#
# Modification history:
#	2008-Oct-31	GDW	Original version
#-----------------------------------------------------------------------

string	image="last"	{ prompt="Name of LRIS alignment image ('last' for most recent)" }
string	logfile="get_seeing.log"	{ prompt="Name of logfile" }
struct	*ilist 

begin
	string	i_image
	string  boxfile="default"
	string	i_boxfile
	string	cursorfile	# file containing cursor input
	string	plotfile	# file containing graphics output input
	string	datafile	# file containing data
	string	parfile		# file containing parameters
	string	tmpprefix	# prefix for temp files
	string	slitname	# slitmask name
	string	dir		# name of directory with image
	real	xc,yc		# centroids
	real	fwhm_px
	real	fwhm_asec
	int		frame=1
	int		n
	int		m
	int		k
	int		nx,ny
	int	debug=0

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
		lastimageb() | scan( i_image)
	}

	# verify that the image exists...
	if( ! imaccess( i_image))
		error( 1, "Image "//i_image//" not found")
	print( "Processing image ", i_image)

	# generate boxfile name as needed...
	i_boxfile = boxfile
	if( i_boxfile == "default" ) {

	    hselect( i_image, "slitname", "yes") | scan( slitname)
	    if ( slitname == "direct" ){
		print( "I can only determine the box locations if you use a slitmask!")
		error( 1, "")
	    }

	    breakname( i_image) | scan( dir)
	    i_boxfile = dir + slitname + '.box'
	    print( "boxfile is ", i_boxfile)
	}

	# verify that box file exists...
	if( ! access( i_boxfile))
	    error( 1, "I can't find the expected boxfile ", i_boxfile)

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
	rimexam.radius = 8.
	rimexam.buffer = 4.
	rimexam.width = 3.
	rimexam.iterations = 1
	rimexam.xorder = 0
	rimexam.yorder = 0
	rimexam.magzero = 25.
	rimexam.beta = INDEF
	rimexam.rplot = 8.
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

	# display the image...
	display( i_image, frame, zsc-, zra-, z1=1000, z2=2500)

	# parse the file...
	ilist = i_boxfile
	n = 0
	while( fscan( ilist, xc, yc) != EOF){

		n+=1
		print( xc, yc, " 100 r", >cursorfile)

		# examine file...
		if ( debug == 1) {
		    printf( "Measuring box at (%d,%d) in image %s\n", xc, yc, i_image)
		    printf( "running imexamine...\n")
		    type( cursorfile)
		}
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
	printf( "\nStar FHMW [px] FWHM [asec]")
	printf( "\n---- --------- -----------\n")
	tcalc( datafile, "c10", "c9*0.135")
	tprint( datafile, col="c9,c10", showhdr-)
	tstat( datafile, "c9", outtab="STDOUT", >"dev$null")
	fwhm_px = tstat.median
	tstat( datafile, "c10", outtab="STDOUT", >"dev$null")
	fwhm_asec = tstat.median
	printf( "%-8s %5.2f  %5.2f\n", "Median", fwhm_px, fwhm_asec)

	# print to logfile...
	if( logfile != "")
		printf( "%-40s FWHM = %5.2f px = %5.2f asec\n", 
			i_image, fwhm_px, fwhm_asec, >>logfile)

	# clean up...
	cl( <parfile)
	delete( plotfile, ver-)
	delete( datafile, ver-)
	delete( parfile, ver-)
	ilist=""
	
end

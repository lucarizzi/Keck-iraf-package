#-----------------------------------------------------------------------
procedure lspecfocus( images)
#-----------------------------------------------------------------------
# keck.lris.lspecfocus
#
# Author:
#	Gregory D. Wirth [wirth@keck.hawaii.edu]
# 
# Purpose:
#	Derive best focus from a series of LRIS arc exposures.
#
# Modification history:
#	2006-Jan-26	GDW		Original version
#	2008-Jan-04	GDW	Fix case of missing INSTRUME in header
#	2008-Jan-05	GDW	Account for change in image filename prefix
#-----------------------------------------------------------------------

string	images="last 7 red"	{ prompt="List of CCD images to correct" }
int		start1=INDEF	{ prompt="Lower edge slit1", min=1}
int		end1=INDEF		{ prompt="Upper edge slit1", min=1}
int		start2=INDEF	{ prompt="Lower edge slit2", min=1}
int		end2=INDEF		{ prompt="Upper edge slit2", min=1}
bool	interactive=yes	{ prompt="Interactively measure slit edges?" }
real	bestfocus		{ prompt="Measured best focus" }

begin
	string	i_images		# internal copy of images
	string	tmpprefix		# prefix for temp files
	string	imfile			# file listing input images
	string	inst			# instrument keyword
	string	inst1=""		# instrument name
	string	im				# image name
	string	side			# side of LRIS to analyze (red/blue)
	string	keyword			# keyword name corresponding to focus
	string	buf				# string buffer
	string	type			# row or column
	string	key				# keypress
	string	word1,word2,word3	# words
	string	test			# test to apply in hselect
	string	tmpfile			# temp file for filenames
	string	tmpfile2		# temp file for filenames
	string	tmpfile3		# temp file for filenames
	int		nimages			# number of input images
	int		nx,ny			# number of columns and rows in imagea
	int		dispaxis		# dispersion axis
	int		ndisp			# number of dispersion samples
	int		prepix			# number of prescan pixels per amp
	int		numamps			# number of amplifiers
	int		xstart			# starting col in readout window
	int		ystart			# starting col in readout window
	int		nxdata			# number of data cols
	int		nydata			# number of data rows
	int		col1,col2		# range of columns to plot
	int		row1,row2		# range of rows to plot
	int		istat			# status
	int		max_regions=2	# maximum number of regions
	int		n_regions		# counter for number of regions
	int		startpix[2]		# starting row/col in regions
	int		endpix[2]		# ending row/col in regions
	int		i,j				# counters
	int		maxpix=150		# maximum pixels per spectral sample
	int		nspectra		# number of spectrak samples for specfocus
	real	x,y				# graphics coordinates
	real	focus			# best focus for region n
	real	sum=0			# sum of focus values
	bool	debug=yes		# enable debug mode?

	# prompt...
	i_images = images

	# lock parameters in cache...
	cache ("sections")
	cache ("tstat")

	# define temp files...
	if ( !defvar("tmp"))
		set tmp = "/tmp"
	tmpprefix = "tmp$lspecfocus"
	imfile = mktemp( tmpprefix)
	tmpfile = mktemp( tmpprefix)
	tmpfile2 = mktemp( tmpprefix)
	tmpfile3 = mktemp( tmpprefix)

	# build input image list...
	istat = fscan( i_images, buf, nimages, side)
	if ( istat==3 && buf=="last") {

		# get last "nimages" images from side "side"...
		if ( side=="red"){
			test = "LRIS$"
		} else if ( side == "blue"){
			test = "LRISBLUE$"
		} else {
			error( 1, "invalid syntax in 'images' parameter; try 'last n red|blue'")
		}

		# create sorted list of images...
		ls( "-1tr", "*.fits", >tmpfile3)

		# select images for this side...
		hselect( "@" // tmpfile3, "$I,INSTRUME", "yes", >tmpfile2)
		match( test, tmpfile2, stop-, print-, meta+ ) \
		| fields( "STDIN", "1", lines="1-", quit-, print-, >tmpfile)

		# select n most recent...
		tail( tmpfile, nlines=nimages, >imfile)

		# count images...
		printf( "Input image list:\n")
		sections( "@"//imfile)

		# remove temp files...
		delete( tmpfile, ver-)
		delete( tmpfile2, ver-)
		delete( tmpfile3, ver-)
	} else {
		sections( i_images, opt="fullname", >imfile)
	}

	# count images...
	nimages = sections.nimages
	if ( nimages < 3 )
		error( 1, "must have at least 3 images to run focus check!")

	# check for all images acquire with same side of LRIS...
	list = imfile
	while( fscan( list, im) != EOF) {

		# verify image access...
		if( ! imaccess( im))
			error( 1, "can't access input image "//im)

		# grab instrument from header...
		hselect( im, "instrume", "yes") | scan( inst)

		# save setting if not defined...
		if( inst1 == "" )
			inst1 = inst

		# stop if not all instrume values are alike...
		if( inst != inst1 )
			error( 1, "all images must be acquired with the same side")
	}

	# check for red or blue side...
	if ( inst1 == "LRIS"){
		if ( debug )
			print( "Using RED side.")
		side = "red"
		keyword = "REDFOCUS"
		dispaxis = 1
		ndisp = 1
		type = "column"
	} else if ( inst1 == "LRISBLUE"){
		if ( debug )
			print( "Using BLUE side.")
		side = "blue"
		keyword = "BLUFOCUS"
		dispaxis = 2
		ndisp = 1
		type = "row"
	} else {
		error( 1, "illegal value for INSTRUME: "//inst1)
	}	

	# check whether we to get slit ends interactively...
	if( interactive){

		# get first image name...
		head( imfile, nlines=1) | scan( im)

		# get image size...
		hselect( im, "naxis1,naxis2", "yes") | scan( nx, ny)
		hselect( im, "prepix", "yes") | scan( prepix)
		hselect( im, "numamps", "yes") | scan( numamps)
		hselect( im, "window", "yes") | \
			translit( "STDIN", ",", " ", del-, coll-) | \
			scan( buf, xstart, ystart, nxdata, nydata)

		# plot rows/columns...
		print( "\nPlotting cut through image.  Please use the graphics window to mark the\nilluminated region(s).")
		if ( side == "red" ) {
			col1 = 2*prepix + 1
			col2 = col1 + nxdata - 1
			pcols( im, col1, col2, wcs="logical", wx1=0., wx2=0.,
				wy1=0., wy2=0., pointmode=no, marker="box",
				szmarker=0.005, logx=no, logy=no, xlabel="wcslabel",
				ylabel="", title="imtitle", xformat="wcsformat", vx1=0.,
				vx2=0., vy1=0., vy2=0., majrx=5, minrx=5, majry=5,
				minry=5, round=no, fill=yes, append=no)
		} else {
			row1 = 1
			row2 = ny
			prows( im, row1, row2, wcs="logical", wx1=0., wx2=0.,
				wy1=0., wy2=0., pointmode=no, marker="box",
				szmarker=0.005, logx=no, logy=no, xlabel="wcslabel",
				ylabel="", title="imtitle", xformat="wcsformat", vx1=0.,
				vx2=0., vy1=0., vy2=0., majrx=5, minrx=5, majry=5,
				minry=5, round=no, fill=yes, append=no)
		}

		# get regions...
		n_regions = 0
		while ( n_regions < max_regions ) {

			# get regions...
			printf( "Mark starting %s for region %d, or 'q' to quit...", 
				type, n_regions)
			istat = fscan( gcur, x, y, buf, key)
			print( "OK")
			if ( key == "q")
				break
			n_regions=n_regions+1
			startpix[n_regions] = nint(x)

			printf( "Mark ending %s, ...", type)
			istat = fscan( gcur, x, y, buf, key)
			print( "OK")
			endpix[n_regions] = nint(x)

		}
	} else {
		n_regions = 0

		# capture values from task params...
		if ( start1 != INDEF && end1 != INDEF ){
			n_regions += 1
			startpix[n_regions] = start1
			endpix[n_regions] = end1
		}

		if ( start2 != INDEF && end2 != INDEF ){
			n_regions += 1
			startpix[n_regions] = start2
			endpix[n_regions] = end2
		}
	}

	# verify that we got valid regions...
	if ( n_regions < 1 )
		error( 1, "no valid regions found")

	# loop over regions...
	for ( i=1 ; i<=n_regions ; i=i+1){

		# swap endpoints?
		if ( endpix[i] < startpix[i]) {
			j = endpix[i]
			endpix[i] = startpix[i]
			startpix[i] = j
		}

		# analyze data...
		nspectra = int((endpix[i]-startpix[i])/maxpix) + 1
		print( "Analyzing data in region ", i, " of ", n_regions)
		specfocus( "@"//imfile, focus=keyword, corwidth=20,
			level=0.5, shifts=yes, dispaxis=dispaxis, nspectra=1,
			ndisp=ndisp, slit1=startpix[i], slit2=endpix[i], 
			logfile="STDOUT") \
			| match( "Best", "STDIN", stop-, print-, meta-) \
			| scan( buf, buf, buf, buf, focus)
		printf( "Region %d: range = %d:%d focus = %d\n", 
			i, startpix[i], endpix[i], focus)
		sum += focus
	}

	print( "Images analyzed:")
	type( imfile)

	bestfocus = sum/n_regions
	printf( "Overall best focus = %d\n", bestfocus)

	# clean up...
	delete( imfile, ver-)
end

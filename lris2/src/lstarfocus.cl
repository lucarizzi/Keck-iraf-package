#-----------------------------------------------------------------------
procedure lstarfocus( side)
#-----------------------------------------------------------------------
# Name:
#	gdwlrispkg.lstarfocus
# 
# Purpose:
#	Determine secondary focus using stellar images in LRIS direct mode.
# 
# Author:
#	Gregory D. Wirth, W. M. Keck Observatory
# 
# Modification history:
#	2001-Jun-04		GDW		Original version
#-----------------------------------------------------------------------

string	side	{ "red", prompt="LRIS camera to focus", enum="red|blue" }
int		n		{ 7, prompt="Number of focus steps", min=5, max=15}
real	center_focus	{ INDEF, prompt="Center focus (default=current)"}
real	incr	{ 0.1, prompt="Focus stepsize [microns]", min=0., max=1.}

begin
	string	i_side				# internal copy of side
	string	cmd=""				# remote command
	string	machine="punaluu"	# summit host
	string	image				# full name of the image acquired
	string	new_image			# root name of the image acquired
	string	root				# rootname of image
	string	keyword="TELFOC"	# keyword giving telescope focus
	string	buf					# string buffer
	string	logfile="lstarfocus.log"	# text output
	int		sbuffer				# width of buffer between annuli [px]
	int		frame=1				# frame for image display
	real	telfoc				# central secondary focus value [mm]
	real	telfoc_new			# derived best value of telfocus
	real	focus1				# starting focus value [mm]
	real	delta=10.			# number of arcsec between images
	real	step				# number of pixels between images
	real	pscale=0.215		# arcsec per pixel
	real	radius				# max radius for photometry
	bool	debug=no			# flag whether to disable datataking for tests
	bool	answer=yes			# user response

	# prompt...
	i_side = side

	# determine central focus...
	if ( center_focus == INDEF ){
		cmd = "telfoc"
		rsh( machine, cmd) | scan( telfoc)
		if( nscan() != 1)
			error( 1, "error reading keyword telfoc on "//machine)
	} else {
		telfoc = center_focus
	}
	printf( "telfoc is %f\n", telfoc)

	# acquire data...
	if( debug){
		cmd = "echo "
	} else {
		cmd = ""
	}
	cmd = cmd + "starfocusloop " + i_side + " " + str(n)
	cmd = cmd + " " + str(telfoc) + " " + str(incr)
	printf( "Taking starfocus image on %s using command:\n\t%s\n",
		machine, cmd)
	rsh( machine, cmd)
	print( "Data acquisition completed.")

	# pause to prevent 'lastimage' from returning wrong value...
	sleep( 5)

	# obtain the name of the file created...
	printf( "Getting image name...")
	if ( i_side == "red" ){
		cmd = "lastimager"
	} else {
		cmd = "lastimageb"
	}
	rsh( machine, cmd) | scan( image)

	# prepend /s to filename if required...
	if( substr( image, 1, 3) != "/s/"){
		image = "/s" + image
	}
	print( image)

	if ( ! imaccess(image))
		error( 1, "Can't access output image.")

	# process the image to ensure that the bias level difference
	# between the two amplifiers will not be a problem...
	print( "Copying image to local directory")
	imcopy( image, ".") | scan( buf, buf, new_image)
	if( nscan() != 3 )
		error( 1, "error copying file to local working directory")
	print( "Processing to subtract overscan")
	lccdproc( new_image, ccdtype="", max_cache=0, noproc-, fixpix-,
		overscan+, trim+, zerocor-, darkcor-, flatcor-, illumcor-,
		fringecor-, readcor-, scancor-, readaxis="line", fixfile="",
		zero="", dark="", flat="", illum="", fringe="", minreplace=1.,
		scantype="shortscan", nscan=1, interactive-, function="leg",
		order=1, sample="*", naverage=1, niterate=1, low_reject=3.,
		high_reject=3., grow=0., rdnoise=INDEF, gain=INDEF)

	# analyze the image...
	focus1 = telfoc - (n-1)/2*incr
	step = delta / pscale
	radius = 0.4*step
	sbuffer = nint(0.1*step)

	if( debug)
		new_image = "stars"

	# analyze data...
	printf( "Processing done --- displaying image in buffer %d\n", frame)
	printf( "Follow these steps to complete analysis:\n")
	printf( "\t-Center the image cursor on the LEFTMOST star image.\n")
	printf( "\t-Press the 'm' key to measure the star images.\n")
	printf( "\t-Wait for the image cursor to return.\n")
	printf( "\t-Press the 'q' key to quit marking and view plot.\n")
	printf( "\t-Delete any bad data points on the plot using the 'd' key.\n")
	printf( "\t-Press the 'q' key to print final results.\n")
	printf( "The best focus will be printed below when done.\n\n")
	beep()
	starfocus( new_image, focus=focus1, fstep=incr, nexposures=n,
		step=step, direction="+column", gap="none", coords="mark1",
		wcs="logical", display+, frame=frame, level=0.5, size="FWHM",
		beta=INDEF, scale=0.21, radius=10., iterations=2,
		sbuffer=5., swidth=5., saturation=60000.,
		ignore_sat+, xcenter=INDEF, ycenter=INDEF, logfile=logfile,
		imagecur="", graphcur="")

	# remove the local copy of the image...
	if( ! debug)
		imdelete( new_image, ver-)

	# get results...
	tail( logfile, n=1) | scan( buf, buf, buf, telfoc_new)
	if( nscan() != 4)
		error( 1, "error reading focus value from file")
	printf( "\nDo you want to set the telescope focus to %.3f mm now? (y/n) [%b]: ",
		telfoc_new, answer)
	beep()
	scan( answer)
	if( answer){
		printf( "Setting telescope focus to %.3f\n", telfoc_new)
		if( debug){
			cmd = "echo "
		} else {
			cmd = ""
		}
		cmd = cmd + "telfoc " + str(telfoc_new)
		rsh( machine, cmd)
	}
	print( "Done.")	
	beep()
end

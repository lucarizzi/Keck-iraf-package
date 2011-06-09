#-----------------------------------------------------------------------
procedure do_check_boxes()
#-----------------------------------------------------------------------
# do_check_boxes / GDW / 2008-Nov-03
# 
# Purpose:
#	Automate the startup of check_boxes by waiting for the image to be
#	written and capturing the image name automatically.
#
# Usage:
#	do_check_boxes
# 
# Arguments:
#	None
# 
# Restrictions:
#	- Should not be started until after the alignment exposure has 
#	started, or else it will analyze the previous image.
# 
#	- box file must be named SLITNAME.box, where SLITNAME is the name of
#	the mask as represented by the SLITNAME keyword, and must reside in 
#	the same directory with the images.
# 
# Modification history:
#	2008-Nov-03	GDW	Original version
#-----------------------------------------------------------------------

string	side="blue"	{ prompt="Side of LRIS acquiring image", enum="red|blue" }

begin
	string	image
	string	slitname
	string	boxfile
	string	dir

	# get name of last blue image...
	print( "Waiting for image to appear...")
	if ( side == "blue" ) {
		wfib( "-n")
		lastimageb | scan( image)
	} else {
		wfir( "-n")
		lastimager | scan( image)
	}

	print( "Image is ", image)
	beep()

	# get name of box file...
	hselect( image//"[0]", "slitname", "yes") | scan( slitname)
	if ( slitname == "direct" ){
		print( "I can only determine the correct box file when you're using a slitmask!")
		print( "Please use the 'check_boxes' command when in direct mode.")
		error( 1, "")
	}
	print( "Slitname is ", slitname)

	# check for box file...
	boxfile = slitname + '.box'
	print( "Boxfile is ", boxfile)
	if ( ! access(boxfile))
		error( 1, "Unable to open box coord file "//boxfile//" -- abort!")

	# run xbox...
	print( "Running check_boxes on image ", image)
	check_boxes( image, boxfile, dolbox+)	
end

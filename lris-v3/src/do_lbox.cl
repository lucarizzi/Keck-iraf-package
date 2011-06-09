#-----------------------------------------------------------------------
procedure do_lbox()
#-----------------------------------------------------------------------
# do_lbox / GDW / 2008-Nov-03
# 
# Purpose:
#	Automate the startup of lbox by waiting for the image to be
#	written and capturing the image name automatically.
#
# Usage:
#	do_lbox
# 
# Arguments:
#	None
# 
# Restrictions:
#	- Should not be started until after the alignment exposure has 
#	started, or else it will analyze the previous image.
# 
#	- box file must be named SLITNAME.box, where SLITNAME is the name of
#	the mask as represented by the SLITNAME keyword.
# 
# Modification history:
#	2008-Nov-03	GDW	Original version
#	2009-Jun-06	GDW	Adapt for lbox
#-----------------------------------------------------------------------

string	side="blue"	{ prompt="Side of LRIS acquiring image", enum="red|blue" }

begin
	string	image
	string	slitname
	string	boxfile
	string	dir
	string	tmpprefix = "DO_LBOX"
	string	null = ""

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
	print( "I can only determine the correct box file when you using a slitmask!")
	print( "Please use the 'lbox' command when in direct mode.")
	error( 1, "")
	}
	print ( "Slitname is ", slitname)

	# check for box file...
##	breakname( image) | scan( dir)
##	boxfile = dir + slitname + '.box'
	boxfile = slitname + '.box'
	if ( ! access(boxfile))
	error( 1, "Unable to open box coord file "//boxfile//" -- abort!")
	print( "Boxfile is ", boxfile)

	# run xbox...
	print( "Running lbox on image ", image)
	lbox( image, boxfile, practice-)	
end

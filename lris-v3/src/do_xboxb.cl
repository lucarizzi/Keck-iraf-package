#-----------------------------------------------------------------------
procedure do_xboxb()
#-----------------------------------------------------------------------
# do_xboxb / GDW / 2008-Nov-03
# 
# Purpose:
#	Automate the startup of xboxb by waiting for the image to be
#	written and capturing the image name automatically.
#
# Usage:
#	do_xboxb
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
#	2009-Jun-02	GDW	Add support for multi-HDU images
#-----------------------------------------------------------------------

bool	overscan	{ no, prompt="Subtract overscan?" }

begin
	string	iimage
	string	slitname
	string	boxfile
	string	dir
	string	tmpprefix = "DO_XBOXB"
	string	null = ""
	string	simpleimage

	# initialize simpleimage...
	simpleimage = null

	# get name of last blue image...
	print( "Waiting for image to appear...")
	wfib( "-n") 
	beep()

	lastimageb | scan( iimage)
	print( "Image is ", iimage)

	# get name of box file...
	hselect( iimage, "slitname", "yes") | scan( slitname)
	if ( slitname == "direct" ){
	print( "I can only determine the correct box file when you using a slitmask!")
	print( "Please use the 'xboxb' command when in direct mode.")
	error( 1, "")
	}
	print ( "Slitname is ", slitname)

	# check for box file...
	breakname( iimage) | scan( dir)
	boxfile = dir + slitname + '.box'
	if ( ! access(boxfile))
	error( 1, "Unable to open box coord file "//boxfile//" -- abort!")
	print( "Boxfile is ", boxfile)

	# convert multi-HDU to simple file if needed...
	count_hdu( iimage) | scan( nhdu)
	if ( nhdu < 0 ) {
		error( 1, "Unable to count HDUs")
	} else if ( nhdu > 0 )
		print( "Creating temporary simple image from ", iimage)
		simpleimage = mktemp(tmpprefix) // ".fits"
		multi2simple( iimage, simpleimage, overscan=overscan)
		iimage = simpleimage
	}

	# run xbox...
	print( "Running xbox_b on image ", iimage)
	xboxb( iimage, boxfile)	

	# clean up...
	if ( simpleimage != null )
		imdelete( simpleimage, ver-)
end

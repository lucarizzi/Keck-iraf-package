#-----------------------------------------------------------------------
procedure do_check_boxesb()
#-----------------------------------------------------------------------
# do_check_boxesb / GDW / 2008-Nov-03
# 
# Purpose:
#	Automate the startup of check_boxesb by waiting for the image to be
#	written and capturing the image name automatically.
#
# Usage:
#	do_check_boxesb
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

begin
    string	image
    string	slitname
    string	boxfile
    string	dir

    print( "Waiting for image to appear...")
    wfib( "-n") 
    beep()

    # get name of last blue image...
    lastimageb() | scan( image)
    print( "Image is ", image)

    # get name of box file...
    hselect( image, "slitname", "yes") | scan( slitname)
    if ( slitname == "direct" ){
	print( "I can only determine the correct box file when you're using a slitmask!")
	print( "Please use the 'check_boxesb' command when in direct mode.")
	error( 1, "")
    }
    print( "Slitname is ", slitname)

    # check for box file...
    hselect( image, "slitname", "yes") | scan( slitname)
    breakname( image) | scan( dir)
    boxfile = dir + slitname + '.box'
    print( "Boxfile is ", boxfile)
    if ( ! access(boxfile))
	error( 1, "Unable to open box coord file "//boxfile//" -- abort!")

    # run xbox...
    print( "Running check_boxesb on image ", image)
    check_boxesb( image, boxfile)	
end

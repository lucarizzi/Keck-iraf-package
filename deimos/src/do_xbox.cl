#-----------------------------------------------------------------------
procedure do_xbox()
#-----------------------------------------------------------------------
# deimos.do_xbox / GDW / 2002-Sep-10
# 
# Purpose:
#	Automate the startup of xbox by waiting for the image to be
#	written and capturing the image name automatically.
#
# Usage:
#	do_xbox
# 
# Arguments:
#	None
# 
# Restrictions:
#	Must be run after the previous exposure has ended and before 
#	the alignment exposure has ended.
# 
# Modification history:
#	2002-Sep-10	GDW	Original version
#	2003-Jun-01	GDW	Added "-n" to wfi
#	2008-Oct-23	GDW	Use full name of last image (prepend "/s")
#-----------------------------------------------------------------------

begin
	string	image

	print( "Waiting for image to appear...")
	wfi( "-n") 
	beep()
	lastimage() | scan( image)
	image = "/s" // image
	print( "Running xbox on image ", image)
	xbox( image)	
end

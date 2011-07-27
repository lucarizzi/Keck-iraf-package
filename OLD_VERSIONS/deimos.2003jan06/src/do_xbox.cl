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
#-----------------------------------------------------------------------

begin
	string	image

	print( "Waiting for image to appear...")
	wfi() 
	beep()
	lastimage( "-tail" ) | scan( image)
	xbox( image)	
end

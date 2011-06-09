#-----------------------------------------------------------------------
procedure count_hdu( image)
#-----------------------------------------------------------------------
# Name:
#	nhdu
#
# Purpose:
#	Return the number of IMAGE extnesions in the image.  A typical
#	"flat" file will have zero.
#
# Author:
#	Gregory D. Wirth, W. M. Keck Observatory
# 
# Modification history:
#	2009-Jun-02     GDW     original version
#-----------------------------------------------------------------------
string	image	{ prompt="Name of image to scan" }
int		nhdu	{ prompt="Number of HDUs found (output)" }

begin
	string	iimage	   # internal copy of image
	int		i

	# prompt...
	iimage = image

	# count number of extensions...
	fxheader( iimage, format_file="", long_header-, count_lines-) \
		| match( "IMAGE", files="STDIN", stop-, print_file_names-, meta-) \
		| count( "STDIN") \
		| scan( i)

	# return value to calling routine...
	nhdu = i
end

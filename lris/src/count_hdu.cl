#-----------------------------------------------------------------------
procedure count_hdu( image)
#-----------------------------------------------------------------------
# Name:
#	nhdu
#
# Purpose:
#	Return the number of IMAGE extensions in the image.  A typical
#	"flat" file will have zero.
#
# Usage:
#	# count number of extensions...
#	count_hdu( image, verb+) | scan(nhdu)
#
# Author:
#	Gregory D. Wirth, W. M. Keck Observatory
# 
# Modification history:
#	2009-Jun-02     GDW     original version
#-----------------------------------------------------------------------
string	image	{ prompt="Name of image to scan" }
bool	verbose	{ no, prompt="Print result?" }
int		nhdu	{ prompt="Number of HDUs found (output)" }

begin
	string	iimage	   # internal copy of image
	int		i=-1

	# prompt...
	iimage = image

	# verify access to image...
	if ( imaccess( iimage)){

		# count number of extensions...
		fxheader( iimage, format_file="", long_header-, count_lines-) \
			| match( "IMAGE", files="STDIN", stop-, print_file_names-, meta-) \
			| count( "STDIN") \
			| scan( i)
	}

	# return value to calling routine...
	nhdu = i
	if( verbose)
		print( nhdu)
end

#-----------------------------------------------------------------------
procedure lris_add_keywords( inimage, outimage)
#-----------------------------------------------------------------------
# Name:
#	lris_add_keywords
# 
# Purpose:
#	Add selected keywords from input image header to output image
# 
# Author:
#	Gregory D. Wirth, W. M. Keck Observatory
# 
# Modification history:
#	2009-Jun-15     GDW     original version
#-----------------------------------------------------------------------
string  inimage		{ prompt="Source of keywords" }
string  outimage	{ prompt="Existing image to modify" }
string	keywords	{ prompt="File listing keywords to copy" }
bool	verbose=no	{ prompt="Give feedback?" }
struct  *ilist

begin
	string	iimage		# internal copy of iimage
	string	oimage		# internal copy of oimage
	#---
	string	keyword		# keyword name
	string	tmpprefix	# prefix for temp files
	string	keywords2	# list of keywords from existing datafile
	string	keywords3	# list of keywords from existing datafile

	# prompt...
	iimage = inimage
	oimage = outimage
	tmpprefix = "tmp$lris_add_keywords"
	keywords2 = mktemp( tmpprefix)
	keywords3 = mktemp( tmpprefix)

	# verify image access...
	if ( ! imaccess( iimage))
		error( 1, "requested input image "//iimage//" not found")
	if ( ! imaccess( oimage))
		error( 1, "requested input image "//oimage//" not found")

	# verify file access...
	if ( ! access( keywords))
		error( 1, "requested keyword list "//keywords//" not found")

	imheader( iimage//"[0]", long+, user+) \
		| match( '=', "STDIN", stop-, print-, meta-, >keywords2)

	# loop over keywords...
	ilist = keywords
	while( fscan( ilist, keyword) != EOF){

		if( verbose)
			print( "  grabbing keyword = ", keyword)
		match( "^"//keyword//"#=", keywords2, stop-, print-, meta+, 
			>>keywords3)
	}
	mkheader( oimage, keywords3, append+, verbose=verbose)

	# clean up...
	ilist = ""
	delete( keywords2, ver-)
	delete( keywords3, ver-)

end

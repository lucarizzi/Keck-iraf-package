########################################################################
procedure breakname( filename)
########################################################################
# cnoc96a.breakname 
#
# Gregory D. Wirth / UVic / 2 May 96
#
# Separates the the directory and extension names from a filename, 
# returning (1) the directory name, (2) the root name, (3) the
# extension.  The directory name ends just before the first character
# of the rootname; that is, at the last "$" or "/" character.  The
# extension begins at the first "." character after the rootname.
#
# If a component is null, the empty string "" is returned.
#
# Usage in scripts:
#
#	string	filename
#	string	dirname
#	string	rootname
#	string	extension
#	...
#	breakname( filename) | scan( dirname, rootname, extension)
#	
#
########################################################################

string	filename	{ prompt="Name of file" }

begin
	string	i_filename	# internal copy of filename
	string	buf		# string buffer
	int	i		# last character of dirname
	int	j		# last character of root
	int	l		# initial length of string
	int	m		# index position of string

	# initialize...
	i_filename = filename
	i = 0
	l = strlen( i_filename)
	j = l

	# locate last character in directory string = i
	buf = i_filename
	m = stridx( "$/", buf)
	while( m>0){
		i += m
		buf = substr( i_filename, i+1, l)
		m = stridx( "$/", buf)
	}

	# locate last character of rootname = j
	m = stridx( ".", buf)
	while( m>0){
		j = i+m-1
		buf = substr( i_filename, i+1, j)
		m = stridx( ".", buf)
	}

	# print dirname...
	if( i>0)
		printf( "%s", substr( i_filename, 1, i)
	else
		printf( "\"\"")

	# print rootname...
	if( i+1<j)
		printf( " %s", substr( i_filename, i+1, j)
	else
		printf( " \"\"")

	# print extension...
	if( j<l)
		printf( " %s\n", substr( i_filename, j+1, l)
	else
		printf( " \"\"\n")
end

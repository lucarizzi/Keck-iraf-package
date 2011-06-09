#-----------------------------------------------------------------------
procedure plotmask(input)
#-----------------------------------------------------------------------
# Name:
#	plotmask
# 
# Purpose:
#	Convert the Lick MONGO output files from DSIMULATOR into an 
#	IGI-compatible format and generate the plot.
#
# Author:
#	Gregory D. Wirth, W. M. Keck Observatory
# 
# Modification history:
#	2003-May-01	GDW	Original version
#-----------------------------------------------------------------------

string	input	{ prompt="Name of mask to plot" }
string	device="stdgraph"	{ prompt="stdplot" }
struct	line
struct	*ilist

begin
	string	i_input	# internal copy of input
	string	word1, word2
	string	tmpprefix
	string	igifile
	int		istat
	int		l

	# check required package...
	if ( !deftask("igi"))
		error( 1, "Must load STSDAS package first.")

	# prompt...
	i_input = input

	# allocate temp file...
	tmpprefix = "tmp$plotmask"
	igifile = mktemp( tmpprefix)

	# parse file...
	ilist = i_input
	while( fscan( ilist, line) != EOF){
		istat = fscan( line, word1, word2)
		if ( word1 == "psland" ){
			next 
		} else if ( word1 =="margin") {
			next
		} else if ( word1 =="submargin") {
			next
		} else if ( word1 =="square") {
			next
		} else if ( word1 =="hard") {
			next
		} else if ( word1 =="color") {
			printf( "color %s\n", word2, >>igifile)
		} else if ( word1 =="putlabel") {
			l = strlen( line)
			printf( "putlabel %s \"%s\"\n", word2, substr(line,12,l), >>igifile)
		} else {
			print( line, >>igifile)
		}
	}
	# reset...
	ilist = ""

	# run igi...
	igi( dev=device, <igifile)

	# delete...
	delete( igifile, ver-)	

end

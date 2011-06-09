#-----------------------------------------------------------------------
procedure ydistcor( input, output, ref)
#-----------------------------------------------------------------------
# nirspec$ydistcor
#
# Purpose:
#	Align emission features with image rows for low-d NIRSPEC spectrum
# 
# Author:
#   Gregory D. Wirth, W. M. Keck Observatory
# 
# Note:
#	Refer to accompanying file ydistcor.html for full documentation.
# 
# Modification history:
#   2001-Jan-04     GDW     Original version
#-----------------------------------------------------------------------

string	input		{ prompt="input image(s)" }
string	output		{ prompt="output image(s)" }
string	ref			{ prompt="Name of reference image" }
int		x1			{ min=1, max=1024, prompt="Starting column to extract" }
int		x2			{ min=1, max=1024, prompt="Ending column to extract" }
real	dy=0.		{ prompt="Feature shift in y from x1 to x2" }
bool	verbose=yes	{ prompt="Give feedback?" }
struct	*ilist
struct	*olist

begin
	string	tmpprefix	# prefix for temp files
	string	i_input		# internal copy of input
	string	i_output	# internal copy of i_output
	string	i_ref		# reference image
	string	apertimage	# image for this aperture
	string	transimage	# transformed image
	string	database="database"	# directory for identify files
	string	sxn			# image section
	string	cursor		# file for graphics cursor input to identify
	string	in			# current input image
	string	out			# current output image
	string	slit		# current list of slit coords
	string	inlist		# list of input images
	string	outlist		# list of output images
	string	solution	# name of solution file
	string	logfile		# output file for tasks
	real	shift		# pixel shift of features between adjacent samples
	int		rows		# number of rows in input image
	int		cols		# number of columns in input image
	int		maxlines=11	# maximum number of independent reident fits
	int		nsum		# number of lines to sum for ident/reident
	int		nstep		# number of lines to step for reident
	int		i			# index
	int		nslits		# number of slit images
	int		istat		# status of scan

	# queries...
	i_input = input
	i_output = output
	i_ref = ref

	# temp files...
	tmpprefix = "YDC"
	cursor = mktemp( tmpprefix)
	transimage = mktemp( tmpprefix)
	inlist = mktemp( tmpprefix)
	outlist = mktemp( tmpprefix)
	apertimage = mktemp( tmpprefix)
	solution = "fc"//i_ref

	# expand input lists...
	sections( i_input, >inlist)
	sections( i_output, >outlist)

	# define disposition of output...
	if( verbose)
		logfile = "STDOUT"
	else
		logfile = ""

	# set up cursor input file for ident...
	print( "y", >>cursor)
	print( "f", >>cursor)
	print( "q", >>cursor)
	print( "q", >>cursor)

	# convert the required global shift into a shift per aperture...	
	nsum = (x2-x1)/maxlines
	nstep = nsum
	shift = dy/(x2-x1)*nstep

	# extract the aperture...
	sxn = "[" + str(x1) + ":" + str(x2) + ",*]"
	imcopy( i_ref//sxn, apertimage, verbose+)

	# get coordinates of any emission features in middle column...
	identify( apertimage, section="middle column",
		database=database, coordlist="", nsum=nsum,
		match=5., maxfeatures=25, ftype="emission",
		fwidth=10, cradius=15, threshold=0., minsep=6,
		func="chebyshev", order=2, sample="*",
		niterate=1, low_rej=3., high_rej=3., grow=0,
		autowrite+, cursor=cursor)
	
	# re-identify them in other columns...
	reidentify( apertimage, apertimage, answer="no",
		nsum=nsum, step=nstep, interactive-, section="middle column",
		newaps+, override+, refit+, trace-, 
		shift=shift, search=0, nlost=INDEF,
		cradius=15, threshold=0., addfeatures-, coordlist="",
		match=-10., maxfeatures=50, minsep=6., logfiles="",
		database=database, plotfile="", verbose+)
	
	# derive a coordinate fit...
	fitcoords( apertimage, fitname="", interactive-,
		combine-, database=database,
		function="chebyshev", xorder=3, yorder=3, logfiles=logfile,
		plotfile="")
	imdelete( apertimage, ver-)
	
	# transform all images...
	# loop over input/output pairs...
	ilist = inlist
	olist = outlist
	while( fscan( ilist, in)!=EOF && fscan( olist, out)!=EOF){

		# check for non-existent image...
		if( ! imaccess( in)){
			print( "  WARNING: "//in//" does not exist! -- skipping")
			next
		}

		# check for clobber...
		if( imaccess( out) && out!=in){
			print( "  WARNING: "//out//" already exists! -- skipping")
			next
		}

		imcopy( in//sxn, apertimage, verbose+)
		print(' ref = ', i_ref)
		transform( apertimage, out, fitname=apertimage,
			database=database, interp="spline3",
			x1=INDEF, x2=INDEF, dx=INDEF, nx=INDEF, xlog-,
			y1=INDEF, y2=INDEF, dy=INDEF, ny=INDEF, ylog-,
			flux+, logfiles=logfile )
		imdelete( apertimage, ver-)
	}
	
	# clean up...
	delete( cursor, ver-)
	delete( inlist, ver-)
	delete( outlist, ver-)
	ilist = ""
	olist = ""

end

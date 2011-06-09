#-----------------------------------------------------------------------
procedure xdistcor( inimage, outimage)
#-----------------------------------------------------------------------
# gdwlrispkg$xdistcor / GDW / 16 Jul 95
#
# Given an input image plus a list of its apertures, this
# tasks extracts the requested apertures are removes any
# x distortion.  Note: the y-distortion must already have been removed
# from the input image.
#
# For each image to be straghtened, the program begins by reading from
# the <slitcoords> file the starting and ending rows occupied by the
# first aperture (the <slitcoords> file is created using the task
# gdwlrispkg.mkslitcoo).  The selected rows are copied to a new image,
# then longslit.identify is used to non-interactively find *any*
# bright emission features in the image.  Next, the reidentify task
# is used to find these same features in other rows of the aperture
# image, co-adding data as needed in order to have not more than
# <maxlines> independent fits for a given aperture.  The fitcoords
# task is used to find a solution which fits the coordinates with a
# polynomial which is second order in x and third order in y.  Then
# the task transform will straighten the image, after which we
# reinsert the aperture back into its original place in the output
# image.
#
# Packages required:
#	noao.twodspec	( identify )
#	images		( imcopy, etc.)
#	ctio		( imcreate )
#
#-----------------------------------------------------------------------

string	inimage		{ prompt="input image(s)" }
string	outimage	{ prompt="output image(s)" }
string	slitcoords \
		{ prompt="list(s) of slit endpoints (1 file per image)" }
struct	*ilist
struct	*olist
struct	*slist
struct	*elist

begin
	string	tmpprefix	# prefix for temp files
	string	iinimage	# internal copy of inimage
	string	ioutimage	# internal copy of ioutimage
	string	apertimage	# image for this aperture
	string	transimage	# transformed image
	string	database="database"	# directory for identify files
	string	sxn		# image section
	string	cursor		# file for graphics cursor input to identify
	string	in		# current input image
	string	out		# current output image
	string	slit		# current list of slit coords
	string	inlist		# list of input images
	string	outlist		# list of output images
	string	slitlist	# list of slit coords
	int	y1		# starting coordinates of aperture
	int	y2		# ending coordinates of aperture
	int	rows		# number of rows in input image
	int	cols		# number of columns in input image
	int	maxlines=11	# maximum number of independent reident fits
	int	nsum		# number of lines to sum for ident/reident
	int	nstep		# number of lines to step for reident
	int	i		# index
	int	nin		# number of input images
	int	nout		# number of output images
	int	nslits		# number of slit images
	int	istat		# status of scan

	# queries...
	iinimage = inimage
	ioutimage = outimage

	# temp files...
	tmpprefix = "XDC"
	cursor = mktemp( tmpprefix)
	transimage = mktemp( tmpprefix)
	inlist = mktemp( tmpprefix)
	outlist = mktemp( tmpprefix)
	slitlist = mktemp( tmpprefix)

	# caches...
	cache( "imgets")
	cache( "sections")

	# expand input lists...
	sections( iinimage, >inlist)
	nin = sections.nimages
	sections( ioutimage, >outlist)
	nout = sections.nimages
	sections( slitcoords, >slitlist)
	nslits = sections.nimages

	# check lists for errors...
	if( nin!=nout)
		error( 1, "different numbers of input and output images!")
	if( nslits!=nin)
		error( 1, "wrong number of slitcoord lists --- must have 1 per image")

	# set up cursor input file for indent...
	print( "y", >>cursor)
	print( "f", >>cursor)
	print( "q", >>cursor)
	print( "q", >>cursor)
	
	# loop over input/output pairs...
	ilist = inlist
	olist = outlist
	slist = slitlist
	while( fscan( ilist, in)!=EOF && fscan( olist, out)!=EOF && \
		fscan( slist, slit)!=EOF){

		# create dummy output image...
		hselect( in, "naxis[12]", "yes") | scan( rows, cols)
		imcreate( out, naxis=2, naxis1=rows, naxis2=cols,
			reference=in, header="copy")
	
		# read apertures from list...
		i = 0
		elist = slit
		while( fscan( elist, y1, y2)!=EOF){
	
			i+=1
			apertimage = mktemp( tmpprefix)
			nsum = (y2-y1)/maxlines
			nstep = nsum
	
			# extract the aperture...
			sxn = "[*," + str(y1) + ":" + str(y2) + "]"
			imcopy( in//sxn, apertimage, verbose+)
	
			# get coordinates of any emission features...
			longslit.identify( apertimage, section="middle line",
				database=database, coordlist="", nsum=nsum,
				match=5., maxfeatures=25, ftype="emission",
				fwidth=10, cradius=5, threshold=100, minsep=6,
				func="chebyshev", order=2, sample="*",
				niterate=1, low_rej=3., high_rej=3., grow=0,
				autowrite+, cursor=cursor)
	
			# re-identify them in other lines...
			longslit.reidentify( apertimage, apertimage, answer="no",
				nsum=nsum, step=nstep, interactive-, section="",
				newaps+, override+, refit+, trace-, nlost=INDEF,
				cradius=5, threshold=100, addfeatures-, coordlist="",
				match=5., maxfeatures=INDEF, minsep=6.,
				database=database, plotfile="", verbose+)
	
			# derive a coordinate fit...
			fitcoords( apertimage, fitname="", interactive-,
				combine-, database=database,
				function="chebyshev", xorder=3, yorder=3,
				plotfile="")
	
			# transform the subraster...
			transform( apertimage, transimage, fitname=apertimage,
				database=database, interp="spline3",
				x1=INDEF, x2=INDEF, dx=INDEF, nx=INDEF, xlog-,
				y1=INDEF, y2=INDEF, dy=INDEF, ny=INDEF, ylog-,
				flux+)
	
			# insert transformed image into rightful place...
			imcopy( transimage, out//sxn, verbose+)
	
			# clean up...
			imdelete( apertimage, ver-)
			imdelete( transimage, ver-)
		}
	}

	# clean up...
	delete( cursor, ver-)
	delete( inlist, ver-)
	delete( outlist, ver-)
	delete( slitlist, ver-)
	elist = ""
	ilist = ""
	olist = ""
	slist = ""

end

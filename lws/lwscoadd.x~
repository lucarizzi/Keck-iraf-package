include <imhdr.h>

define	SZ_LWS		128			# define LWS detector size
define	LINE		2			# dimension representing detector lines
define	CHOP_BEAM	3			# dimension representing chop beams
define	NOD_BEAM	5			# dimension representing nod beams
define	CHOP_ON		2			# index number of the chop "on" beam
define	CHOP_OFF	1			# index number of the chop "off" beam
define	NOD_ON		1			# index number of the nod "on" beam
define	NOD_OFF		2			# index number of the nod "off" beam

#-----------------------------------------------------------------------
procedure t_lwscoadd ()
#-----------------------------------------------------------------------
# Purpose:
#	Wrapper to handle user parameter I/O for lwscoadd().
#	Allows use of multiple input and output images.
#-----------------------------------------------------------------------

char	inimlist[SZ_LINE]		# list of input images
char	outimlist[SZ_LINE]		# list of output images
bool	verbose					# print operations?

char	inimage[SZ_FNAME]		# image with isophotal data
char	outimage[SZ_FNAME]		# image to contain output data
char	dirname1[SZ_FNAME]		# directory name
char	dirname2[SZ_FNAME]		# directory name
int		root_len				# length of root name
pointer ilist					# pointer to input image list
pointer olist					# pointer to output image list

pointer imtopen()				# open image list
bool	clgetb()				# get boolean input
int		imtgetim()				# get next entry from image list
int		imtlen()				# return number of elements in list
int		isdirectory()			# is this filename a directory?
int		fnldir()				# 

begin

	# query...
	call clgstr( "input", inimlist, SZ_LINE) 
	call clgstr( "output", outimlist, SZ_LINE) 
	verbose = clgetb( "verbose")

	# check whether the output string is a directory...
	if( isdirectory( outimlist, dirname2, SZ_FNAME) > 0) {

		ilist = imtopen( inimlist)
		while( imtgetim( ilist, inimage, SZ_FNAME) != EOF) {

			# Strip the image section first because fnldir recognizes it
			# as part of a directory.  Place the input image name
			# without a directory or image section in string dirname1.
			call get_root( inimage, outimage, SZ_FNAME)
			root_len = fnldir( outimage, dirname1, SZ_FNAME)
			call strcpy( outimage[root_len + 1], dirname1, SZ_FNAME)

			call strcpy( dirname2, outimage, SZ_FNAME)
			call strcat( dirname1, outimage, SZ_FNAME)
			call lwscoadd( inimage, outimage, verbose)
		}

		call imtclose( ilist)

	} else {

		# open file lists...
		ilist = imtopen( inimlist)
		olist = imtopen( outimlist)			

		# copy input images to output images...
		if( imtlen( ilist) == imtlen( olist)){	

			while( imtgetim( ilist,  inimage, SZ_FNAME) != EOF && 
				   imtgetim( olist, outimage, SZ_FNAME) != EOF)
				call lwscoadd( inimage, outimage, verbose)

		# clobber input images...
		} else if( imtlen( olist)==0){	

			while( imtgetim( ilist, inimage, SZ_FNAME) != EOF)
				call lwscoadd( inimage, inimage, verbose)

		} else {

			call imtclose( ilist)
			call imtclose( olist)
			call error( 1, "Numbers of input and output images are unequal.")
		}	

		# clean up...
		call imtclose( ilist)
		call imtclose( olist)

	}

end

#-----------------------------------------------------------------------
procedure lwscoadd( input, output, verbose)
#-----------------------------------------------------------------------
# lwscoadd / GDW / 27 Jun 1999
# 
# Purpose:
#	To combine properly the various dimensions of an LWS image.  The
#	dimensions are:
#		1) detector columns
#		2) detector rows
#		3) chpbeams
#		4) savesets
#		5) nodbeams
#		6) nodsets
# 
# Modification history:
#   1999-Jun-27  GDW  Original version
#   1999-Oct-26  GDW  Inserted code to check for clobber
#   1999-Oct-27  GDW  Included code to zero the data array
#-----------------------------------------------------------------------

char	input[ARB]				# name of input image
char	output[ARB]				# name of output image
bool	verbose					# print progress?

char	errmsg[SZ_LINE]			# text of error message
char	original[SZ_FNAME]		# name of temporary image
char	errtxt[SZ_LINE]			# error message
long	data[SZ_LWS,SZ_LWS]		# data array
long	iv[IM_MAXDIM]			# vector of input image dimensions
long	pm						# plus or minus value
int		i,j,k,l					# counters
pointer	in						# pointer to input image
pointer	out						# pointer to output image
pointer	lbuf					# line buffer

pointer	immap()
bool	envgetb()				# get boolean environment variable from CL
bool	streq()					# string comparison function
int		imgnll(), impnll()

begin

	# get input image name...
	in = immap( input, READ_ONLY, 0)

	# verify datatype...
	if( IM_PIXTYPE( in) != TY_LONG){
		call imunmap( in)
		call sprintf( errtxt, SZ_LINE, "image data type is not LONG (%s)" )
			call pargstr( input)
		call error( 1, errtxt)
	}

	# verify dimensions...
	if( IM_NDIM(in) != 6 ){
		call imunmap( in)
		call sprintf( errtxt, SZ_LINE, "image is not 6-dimensional (%s)" )
			call pargstr( input)
		call error( 1, errtxt)
	}
	
	# verify image size...
	if( IM_LEN(in,1) != SZ_LWS || IM_LEN(in,2) != SZ_LWS ){
		call imunmap( in)
		call sprintf( errtxt, SZ_LINE, "image has wrong size (%s)" )
			call pargstr( input)
		call error( 1, errtxt)
	}
	if( IM_LEN(in,CHOP_BEAM) > 2 ){
		call imunmap( in)
		call sprintf( errtxt, SZ_LINE, "image has more than 2 chopbeams (%s)" )
			call pargstr( input)
		call error( 1, errtxt)
	}
	if( IM_LEN(in,NOD_BEAM) > 2 ){
		call imunmap( in)
		call sprintf( errtxt, SZ_LINE, "image has more than 2 nodbeams (%s)" )
			call pargstr( input)
		call error( 1, errtxt)
	}
	
	# give report...
	if( verbose){
		call eprintf( "  Coadding %s -> %s\n")
			call pargstr( input)
			call pargstr( output)
	}

	# check for clobber condition...
	if( streq( input, output) && ! envgetb( "clobber")){
		call sprintf( errmsg, SZ_LINE, "requested operation would clobber existing file %s")
			call pargstr( input)
		call error( 1, errmsg)
	}

	# allow clobbering of input image...
	call xt_mkimtemp( input, output, original, SZ_FNAME)

	# erase the data buffer...
	do j=1,SZ_LWS {
		do i=1,SZ_LWS {
			data[i,j] = long(0)
		}
	}

 	# open output image and reset the number of dimensions...
 	out = immap( output, NEW_COPY, in)
	IM_NDIM( out) = 2

 	# reset vector of image dimensions...
 	call amovkl( long(1), iv, IM_MAXDIM)
	j = int(iv[LINE])
	k = iv[CHOP_BEAM]
	l = iv[NOD_BEAM]

	# loop through the input image...
	while( imgnll( in, lbuf, iv) != EOF){

		# presume line is to be added by default...
		pm = long(1)

		# determine whether the chop beam is to be added or subtracted...
		if( k == CHOP_OFF){ pm = -pm }

		# determine whether the nod beam is to be added or subtracted...
		if( l == NOD_OFF ){ pm = -pm }

		# add this line into the final array...
		do i=1,SZ_LWS { data[i,j] = data[i,j] + pm*Meml[lbuf+i-1]}

		# set line index for the next line...
		j = int(iv[LINE])
		k = iv[CHOP_BEAM]
		l = iv[NOD_BEAM]

	}

	# reset vector of image dimensions...
	call amovkl( long(1), iv, IM_MAXDIM)
	j = int(iv[LINE])

	# write out the image
	while( impnll( out, lbuf, iv) != EOF){
		do i=1,SZ_LWS { Meml[lbuf+i-1] = data[i,j] }
		j = int(iv[LINE])
	}

	# close images...
	call imunmap( in)
	call imunmap( out)
	call xt_delimtemp( output, original)
end

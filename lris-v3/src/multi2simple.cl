#-----------------------------------------------------------------------
procedure multi2simple( inimage, outimage)
#-----------------------------------------------------------------------
# Name:
#	multi2simple
# 
# Purpose:
#	Convert multi-HDU LRIS FITS image to "old-style" LRIS simple FITS image
# 
# Author:
#	Gregory D. Wirth, W. M. Keck Observatory
#
# Restrictions:
#	Must be run from a directory to which the user has write access
# 
# Examples:
#	1) Convert a single image:
#		multi2simple b091011_0000.fits flat_b091011_0000.fits
# 
#	2) Convert a set of images:
#		multi2simple @list.in @list.out
# 
# Modification history:
#	2009-May-29	GDW     original version
#	2009-Jun-19	GDW 	added trim option
#	2009-Jul-02	GDW		fixed overscan processing for images with 
#							different datasec
#	2009-Nov-25	GDW		enable multiple input/output images
#	2011-Apr-15	GDW		Fix treatment of prescan in binned case; 
#						add npskip param
#-----------------------------------------------------------------------
string  inimage		{ prompt="Name of input LRIS Multi-HDU image(s)" }
string  outimage	{ prompt="Name of output LRIS simple image(s)" }
bool	overscan=no	{ prompt="Subtract overscan?" }
bool	header=no	{ prompt="Copy image header?" }
bool	trim=yes	{ prompt="Discard pre- and overscan columns?" }
bool	verbose=yes	{ prompt="Give feedback?" }
bool	debug=no	{ prompt="Print debugging output?" }
int		npskip=3	{ prompt="Number of initial postscan pixels to discard", min=0 }
struct	*ilist
struct	*olist

begin
	string	inimages	# internal copy of inimage
	string	outimages	# internal copy of outimage
	#----
	string	iimage		# internal copy of iimage
	string	oimage		# internal copy of oimage
	string	ifile		# file listing input images
	string	ofile		# file listing output images
	string	tmpprefix	# prefix for temp files
	string	inlist		# file listing input sections
	string	offsets		# file listing offsets
	string	instrume	# name of instrument
	string	detsec		# output detector section
	string	datasec		# input (data) section
	string	ccdsec		# detector region corresponding to data
	string	ccdsum		# on-chip binning
	string	s1			# string buffer
	string	isxn		# image section
	string	biassec		# bias section
	string	tempimage1	# name of temp image
	string	tempimage2	# name of temp image
	string	outtype="ushort"	# default output type
	string	keywords	# file listing keywords
	string	im			# name of input image...
	string	keyword		# keyword name
	string	logfile		# logfile name
	int		i
	int		nhdu		# number of HDUs
	int		ix1,ix2
	int		iy1,iy2
	int		jx1,jx2
	int		xbin, ybin
	int		ibuf
	int		xoffset
	int		yoffset=0
	int		nx
	int		nxtotal
	int		x0=32768
	int		x1=-32768
	int		xend
	int		xbias1,xbias2
	int		precol
	int		postpix
	int		nin			# number of input images
	int		nout		# number of output images

	# prompt...
	inimages = inimage
	outimages = outimage

	# allocate temp file names...
	tmpprefix = "tmp$multi2simple"
	inlist = mktemp( tmpprefix)
	offsets = mktemp( tmpprefix)
	keywords = mktemp( tmpprefix)
	ifile = mktemp( tmpprefix)
	ofile = mktemp( tmpprefix)
	tempimage1 = mktemp( tmpprefix) // ".fits"
	tempimage2 = mktemp( tmpprefix) // ".fits"

	# expand input and output templates...
	sections( inimages, >ifile)
	nin = sections.nimages
	if ( nin == 0)
		error( 1, "no images in input list")

	sections( outimages, >ofile)
	nout = sections.nimages
	if ( nout == 0)
		error( 1, "no images in output list")

	# verify equal numbers of elements...
	if ( nout != nin )
		error( 1, "input and output image lists are different lengths")

	# loop over image lists...
	ilist = ifile
	olist = ofile
    while( fscan( ilist, iimage) != EOF && fscan( olist, oimage) != EOF) {
		if( debug)
			print( "start of loop")

		if( verbose)
			print( iimage, " -> ", oimage)

		# save original input image name...
		im = iimage
			
		# verify input...
		if ( ! imaccess( im))
		error( 1, "requested input image "//im//" not found")

		# verify non-existence of output image...
		if ( imaccess( oimage))
			error( 1, "operation would clobber existing output image "//oimage)
	
		# count number of extensions...
		count_hdu( im, verb+) | scan(nhdu)
		if( debug)
			print( "found ", nhdu, " HDU")
		if( nhdu < 1 )
			error( 1, "no valid HDUs found!")
	
		# check binning...
		hselect( im//"[0]", "BINNING", "yes") \
			| translit( "STDIN", ",", " ", del-, col-) \
			| scan( xbin, ybin)
		if( nscan() < 2 )
			error( 1, "Unable to read BINNING from "//im)
		if( debug)
			print( "xbin=", xbin, " ybin=", ybin)

		# create temp image...
		copy( im, tempimage1, verb-)
		im = tempimage1
	
		# fix headers...
		for ( i=1 ; i<=nhdu ; i+=1 ) {
	
			s1="["//i//"]"
	
	 		# get detsec...
			hselect( im//s1, "DETSEC", "yes") \
				| translit( "STDIN", "[:,]", " ", del-, col-) \
				| scan( ix1, ix2, iy1, iy2)
			if( nscan() < 4 )
				error( 1, "Unable to read DETSEC from "//im//s1)
	
			# correct x coords as needed...
			if( ix1 > 2048){
				ix1 -= 2048
				ix2 -= 2048
			}
	
			# construct CCDSEC...
			ccdsec = "[" + ix1 + ":" + ix2 + "," + iy1 + ":" + iy2 + "]"
			hedit( im//s1, "ccdsec", ccdsec, add+, delete-, ver-, 
				show-, update+)
			if( debug)
				print( "ccdsec=", ccdsec)
	
			# construct CCDSUM...
			ccdsum = xbin//" "//ybin
			hedit( im//s1, "ccdsum", ccdsum, add+, delete-, ver-, 
				show-, update+)
			if( debug)
				print( "ccdsum=", ccdsum)
	
		}
	
		# loop for removing overscan...
		if ( overscan ) {
	
			# read postscan values...
			s1="[0]"
			keyword = "POSTPIX"
			hselect( im//s1, keyword, "yes") | scan( postpix)
			if( nscan() != 1 )
				error( 1, "ERROR reading "//keyword//" from "//im)
	
			# build input list for DATA...
			for ( i=1 ; i<=nhdu ; i+=1 ) {
	
				s1="["//i//"]"
	
				# get image size...
				hselect( im//s1, "NAXIS1", "yes") | scan( xbias2)
	
				# build overscan region...
				xbias1 = xbias2 - postpix/xbin + 1 + npskip
				biassec = "[" + xbias1 + ":" + xbias2 + "]"
	
				# save in header...
				hedit( im//s1, "biassec", biassec, add+, delete-, ver-, 
					show-, update+)

				# debugging output...
				if( debug)
					print( "biasec[",i,"] = ", biassec)
			}
	
			# remove overscan...
			ccdproc( im, output=tempimage2, bpmasks="", ccdtype="", noproc=no,
				xtalkcor=no, fixpix=no, overscan=yes, trim=no, zerocor=no,
				darkcor=no, flatcor=no, sflatcor=no, split=no, merge=no,
				xtalkfile="", fixfile="", saturation="INDEF", sgrow=0,
				bleed="INDEF", btrail=20, bgrow=0, biassec="image",
				trimsec="", zero="", dark="", flat="", sflat="",
				minreplace=1., interactive=no, function="legendre", order=1,
				sample="*", naverage=1, niterate=1, low_reject=3.,
				high_reject=3., grow=0., fd="", fd2="")
	
			# the debiased image is now the input image for the next step...
			im = tempimage2
			outtype = "real"
		}
	
		# build input list for DATA...
		for ( i=1 ; i<=nhdu ; i+=1 ) {
	
			s1="["//i//"]"
	
			if ( debug ){
				printf( "datasec=")
				hselect( im//s1, "DATASEC", "yes")
			}
	
			# get datasec...
			hselect( im//s1, "DATASEC", "yes") \
				| translit( "STDIN", "[:,]", " ", del-, col-) \
				| scan( ix1, ix2, iy1, iy2)
	
			# get detsec...	
			hselect( im//s1, "DETSEC", "yes") \
				| translit( "STDIN", "[:,]", " ", del-, col-) \
				| scan( jx1, jx2)
	
			# measure output size...
			nx = ix2 - ix1 + 1
	
			# if detsec is flipped, then flip input section...
			if ( jx1 > jx2 ) {
				ibuf=ix1
				ix1=ix2
				ix2=ibuf
				xoffset = jx2
			} else {
				xoffset = jx1
			}
	
			# correct offset for binning factor...
			xoffset = xoffset / xbin
			yoffset = yoffset / ybin
	
			# save the lowest value of xoffset...
			if ( xoffset < x0 )
				x0 = xoffset
	
			# save the highest value of xend...
			xend = xoffset + nx - 1
			if ( xend > x1 )
				x1 = xend
	
			# define input section...
			isxn = "[" + ix1 +":" + ix2 + "," + iy1 + ":" + iy2 + "]"
	
			# save input section...
			print( im//s1//isxn, >>inlist)
	
			# save offset...
			print( xoffset, yoffset, >>offsets)
		}
	
		# add pre- and post-scan columns if requested...
		if ( ! trim ) {	
			s1 = "[0]"
	
			# read prescan values...
			keyword = "PRECOL"
			hselect( im//s1, keyword, "yes") | scan( precol)
			if( nscan() != 1 )
				error( 1, "ERROR reading "//keyword//" from "//im)
	
			# read postscan values...
			keyword = "POSTPIX"
			hselect( im//s1, keyword, "yes") | scan( postpix)
			if( nscan() != 1 )
				error( 1, "ERROR reading "//keyword//" from "//im)
	
			# make pre-scan offsets relative to the start of the first 
			# DATA section...
			xoffset = x0
	
			# add entries to inlist and offsets for PRECOL...
			for ( i=nhdu ; i>=1 ; i-=1 ) {
	
				s1="["//i//"]"
				ix1 = 1
				ix2 = precol
				isxn = "[" + ix1 +":" + ix2 + ",*]"
	
				# save input section...
				print( im//s1//isxn, >>inlist)
	
				# save offset...
				xoffset -= precol
				print( xoffset, yoffset, >>offsets)
			}
	
			# make post-scan offsets relative to the end of the last
			# DATA section...
			xoffset = x1 + 1
	
			# add entries to inlist and offsets for POSTPIX...
			for ( i=1 ; i<=nhdu ; i+=1 ) {
	
				s1="["//i//"]"
			
				# get image size...
				keyword = "naxis1"
				hselect( im//s1, keyword, "yes") | scan( nx)
				if( nscan() != 1 )
					error( 1, "ERROR reading "//keyword//" from "//im)
			
				# postscan pixels come from the end of the image...	
				ix1 = nx - postpix/xbin + 1
				ix2 = nx
				isxn = "[" + ix1 +":" + ix2 + ",*]"
	
				# save input section...
				print( im//s1//isxn, >>inlist)
	
				# save offset...
				print( xoffset, yoffset, >>offsets)
				xoffset += postpix/xbin
			}		
		}
	
		# create new image...
		if( debug) {
			printf(" inlist:\n")
			type( inlist)
			logfile = "STDOUT"
		} else {
			logfile = "dev$null"
		}
		imcombine( input="@"//inlist, output=oimage, headers="",
			bpmasks="", rejmasks="", nrejmasks="", expmasks="", sigmas="",
			logfile=logfile, combine="average", reject="none", project=no,
			outtype=outtype, outlimits="", offsets=offsets,
			masktype="none", maskvalue=0., blank=0., scale="none",
			zero="none", weight="none", statsec="", expname="",
			lthreshold=INDEF, hthreshold=INDEF, nlow=1, nhigh=1, nkeep=1,
			mclip=yes, lsigma=3., hsigma=3., rdnoise="0.", gain="1.",
			snoise="0.", sigscale=0.1, pclip=-0.5, grow=0.)  
	
		# remove DATASEC (otherwise, ds9 fails to display image properly)...
		hedit( oimage, "datasec", add-, delete+, ver-, show-, update+)
	
		# copy header info...
		if ( header ){
			hselect( im//"[0]", "INSTRUME", "yes") | scan( instrume)
			keywords = "lris$dat/lris_add_keywords." // instrume // ".dat"
			if( debug)
				print( "  Adding keywords listed in file ", keywords)
			lris_add_keywords( im, oimage, keywords=keywords, verb-)
		}
	
		# clean up...
		delete( inlist, ver-)
		delete( offsets, ver-)
		imdel( tempimage1, ver-)
		if ( im == tempimage2 )
			imdel( tempimage2, ver-)

		if( debug)
			print( "end of loop")
	}
	
	# more cleanup...
	delete( ifile, ver-)
	delete( ofile, ver-)

	# clear structs...
	ilist = ""
	olist = ""
		
end

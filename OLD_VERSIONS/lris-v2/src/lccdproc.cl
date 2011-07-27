#-----------------------------------------------------------------------
procedure lccdproc( images)
#-----------------------------------------------------------------------
# wmkolris.lccdproc
#
# Author:
#   Gregory D. Wirth [wirth@keck.hawaii.edu]
#
# Purpose:
#	Implement CCDPROC image processing for LRIS, using custom overscan 
#	subtraction routine.
# 
# Modification history:
#	16 Jun 1998		GDW		Original version
#   24 Feb 1999		GDW		Added parameters for bias subtraction
#-----------------------------------------------------------------------

string	images		{ "", prompt="List of CCD images to correct" }
string	ccdtype		{ "", prompt="CCD image type to correct" }
int		max_cache	{ 0, prompt="Maximum image caching memory (in Mbytes)", min=0}
bool	noproc		{ no, prompt="List processing steps only?\n" }
bool	fixpix		{ yes, prompt="Fix bad CCD lines and columns?" }
bool	overscan	{ yes, prompt="Apply overscan strip correction?" }
bool	trim		{ yes, prompt="Trim the image?" }
bool	zerocor		{ yes, prompt="Apply zero level correction?" }
bool	darkcor		{ yes, prompt="Apply dark count correction?" }
bool	flatcor		{ yes, prompt="Apply flat field correction?" }
bool	illumcor	{ no, prompt="Apply illumination correction?" }
bool	fringecor	{ no, prompt="Apply fringe correction?" }
bool	readcor		{ no, prompt="Convert zero level image to readout correction?"}
bool	scancor		{ no, prompt="Convert flat field image to scan correction?\n" }

string	readaxis	{ "line", prompt="Read out axis (column|line)", \
						enum="column|line"}
string	fixfile		{ "", prompt="File describing the bad lines and columns" }
string	zero		{ "", prompt="Zero level calibration image" }
string	dark		{ "", prompt="Dark count calibration image" }
string	flat		{ "", prompt="Flat field images" }
string	illum		{ "", prompt="Illumination correction images" }
string	fringe		{ "", prompt="Fringe correction images" }
real	minreplace	{ 1., prompt="Minimum flat field value" }
string	scantype	{ "shortscan", prompt="Scan type (shortscan|longscan)", \
						enum="shortscan|longscan" }
int		nscan		{ 1, prompt="Number of short scan lines\n", min=1}
bool	interactive	{ no, prompt="Fit overscan interactively?"}
string	function	{ "legendre", prompt="Fitting function" }
int		order		{ 1, prompt="Number of polynomial terms or spline pieces", min=1 }
string	sample		{ "*", prompt="Sample points to fit" }
int		naverage	{ 1, prompt="Number of sample points to combine" }
int		niterate	{ 1, prompt="Number of rejection iterations", min=0 }
real	low_reject	{ 3., prompt="Low sigma rejection factor", min=0. }
real	high_reject	{ 3., prompt="High sigma rejection factor", min=0. }
real	grow		{ 0., prompt="Rejection growing radius", min=0. }
real    rdnoise		{ INDEF, prompt="Readout Noise [electrons]"}
real    gain		{ INDEF, prompt="Gain [electrons/ADU]"}

begin
	string	i_images		# internal copy of images
	string	tmpprefix		# prefix for temp files
	string	imfile			# file listing input images

	# prompt...
	i_images = images

	# define temp files...
	tmpprefix = "lccdproc"
	imfile = mktemp( tmpprefix)

	# expand filename template...
	ccdlist( i_images, ccdtype=ccdtype, names+, long-, ccdproc="",
		>imfile)

	# remove bias via overscan...
	if( overscan)
		lrisbias( "@"//imfile, rdnoise=rdnoise, gain=gain, 
			dispaxis=-1, observat="keck", noproc=noproc,
			interactive=interactive, function=function, order=order, 
			sample=sample, naverage=naverage, niterate=niterate, 
			low_reject=low_reject, high_reject=high_reject, grow=grow)

	# perform remaining processing...
	ccdproc( i_images, ccdtype=ccdtype, max_cache=max_cache,
		noproc=noproc, fixpix=fixpix, overscan-, trim=trim,
		zerocor=zerocor, darkcor=darkcor, flatcor=flatcor, illumcor=illumcor,
		fringecor=fringecor, readcor=readcor, scancor=scancor,
		readaxis=readaxis, fixfile=fixfile, biassec="", trimsec="",
		zero=zero, dark=dark, flat=flat, illum=illum, fringe=fringe,
		minreplace=minreplace, scantype=scantype, nscan=nscan)

	# clean up...
	delete( imfile, ver-)
end

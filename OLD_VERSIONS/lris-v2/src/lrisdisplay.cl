#-----------------------------------------------------------------------
procedure lrisdisplay( image, frame)
#-----------------------------------------------------------------------
# Name:
#	lrisdisplay
# 
# Purpose:
#	Send LRIS mosaic image to IRAF imtool
# 
# Author:
#   Gregory D. Wirth, W. M. Keck Observatory
# 
# Modification history:
#   2002-Jun-12     GDW     original version
#   2002-Aug-21     GDW     Changed to using file in tmp$
#-----------------------------------------------------------------------
string  image		{ prompt="Name of input LRIS image" }
int	frame		{ prompt="frame to be written into" }
bool	check=no
bool	onepass=yes
bool	erase=yes
bool	border_erase=no
bool	select_frame=yes
bool	zscale=yes
real	contrast=0.25
bool	zrange=yes
int	order=0
real	z1=0
real	z2=1000
string	ztrans="linear"
string	lutfile=""

begin
    string  i_image     # internal copy of image
    int     i_frame	# frame number
    #----
    string  tmpprefix   # prefix for temp files
    real    iz1
    real    iz2
    int     nvidinp	# number of HDUs in LRIS image
    int     i		# counter

    # prompt...
    i_image = image
    i_frame = frame
    iz1 = z1
    iz2 = z2

    # verify input...
    if ( ! imaccess( i_image))
        error( 1, "requested input image "//i_image//" not found")

    # if full intensity range is desired, check scaling limits manually on 
    # the first HDU only and apply these to all of the other HDUs...
    if ( !zscale && zrange){
	print( "Determining min/max values off of HDU 1...")
	minmax( i_image//"[1]", force+, update-, verb-)
	iz1 = minmax.minval
	iz2 = minmax.maxval
    }

    # display image...
    mscdisplay( i_image, i_frame, xgap=0, ygap=0,
	check=check, onepass=onepass, erase=erase, border_erase=border_erase,
	select_frame=select_frame, zscale=zscale, contrast=contrast,
	zrange-, order=order, z1=iz1, z2=iz2, ztrans=ztrans, lutfile=lutfile,
	process-)

end

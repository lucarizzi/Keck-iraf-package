# QBOX: Quick box alignment script -- calls maskalign

procedure qbox (image, box_pos)

file	image {prompt = "input image"}
file	box_pos {prompt = "input box positions"}
real	fwhm {prompt = "FWHM of star (pixels)"}
real	xsz {prompt = "x-size of box (pixels)"}
real	ysz {prompt = "y-size of box (pixels)"}
int	prescan {prompt = "prescan included in input x-coordinates (pixels)"}
bool	invert {prompt = "reverse orientation? (N-top, E-right)"}

begin
{

# define local variables
file	inp, inp2
real	rotposn		# position of rotator on sky
real	angle		# position angle of image
file	qtmp1

# Construct names
inp = image
if (access (inp//".imh"))
	inp = inp//".imh"
else if (access (inp//".hhh"))	
	inp = inp//".hhh"

# Modified 5 April to allow [-*,*] to be appended
# if (! access (inp)) {
#	beep
#	print ""
#	print ("Error: "//inp//" does not exist ! -- exiting")
#	bye
# }

inp2 = box_pos
if (! access (inp2)) {
	beep
        print ""
        print ("Error: "//inp2//" does not exist ! -- exiting")
        bye
}

qtmp1 = mktemp ("tmpq")

mbox (inp, inp2, qtmp1, xsz=xsz, ysz=ysz, xfwhm=fwhm, yfwhm=fwhm, prescan=presca
n)

# Get the solution (can also call ALIGN)
#
print ("\nMake sure MASKALIGN parameters are correct:^G")
print (maskalign.zpa, maskalign.pixscale, maskalign.xrot, maskalign.yrot, "(LRIS
 defaults: +90.  0.215  1066  1024)\n")

# Get the rotator angle --> PA
hselect (inp, "ROTPOSN", yes) | scan (rotposn)
if (rotposn <= 90.) {
	angle = rotposn + 90.
} else {
	angle = rotposn -270.
}
print ("\nPosition angle found:  ", angle)

maskalign (qtmp1, pa=angle, invert=invert, rotposn=rotposn)
delete (qtmp1, veri-)

}
end

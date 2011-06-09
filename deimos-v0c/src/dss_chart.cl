# DSS_CHART: Quick script to convert DSS image to DEIMOS guider image

procedure dss_chart (in_image, out_image)

file	in_image {prompt = "input image"}
file	out_image {prompt = "output image"}
real	pa_mask {prompt = "pa of mask (deg)"}

begin
{

# define local variables
file	inp, out
real	ang, mag
file	qtmp1

# Construct names
inp = in_image
if (access (inp//".imh"))
	inp = inp//".imh"
else if (access (inp//".hhh"))	
	inp = inp//".hhh"
else if (access (inp//".fit"))	
	inp = inp//".fit"
else if (access (inp//".fits"))	
	inp = inp//".fits"

if (! access (inp)) {
	beep
	print ""
	print ("Error: "//inp//" does not exist ! -- exiting")
	bye
}

out = out_image

qtmp1 = mktemp ("tmpq")

ang = 91.4 - pa_mask
mag = 0.203

print ("transforming "//inp//" --> "//out)
imlintran (inp, qtmp1, ang, ang, mag, mag, ncol=1024, nlin=1024, interp="linear", verb-)
imcopy (qtmp1//"[*,-*]", out, verb-)

imdel (qtmp1)

}
end

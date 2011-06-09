#------------------------------------------------------------------------
procedure add_fits_chunks( in, chunks, out)
#------------------------------------------------------------------------

string	in	{ prompt="Input image(s)" }
string	chunks	{ "default", prompt="FITS chunks to insert" }
string	out	{ prompt="Output image(s)" }
struct  *ilist

begin
	string	a,b,c
	string	buf
	string	barcode
        string  inlist          # file listing input images
	string	tmpprefix	# temp file prefix

	# prompt...
	a = in
	b = chunks
	c = out

	# allocate temp files...
        tmpprefix = "tmp$add_fits_chunks"
        inlist = mktemp( tmpprefix)

	# check for default...
	if ( b == "default" ) {
		hselect( a//"[0]", "slmskbar", "yes") | scan( barcode)
		b = "Barcode." + barcode + ".fits"
	}

	# build list of files to copy...
	buf = a
	buf = buf + "," + b + "[1]"
	buf = buf + "," + b + "[2]"
	buf = buf + "," + b + "[3]"
	buf = buf + "," + b + "[4]"
	buf = buf + "," + b + "[5]"
	buf = buf + "," + b + "[6]"
	buf = buf + "," + b + "[7]"

	# copy image...
	fxcopy( buf, c, new+)

end

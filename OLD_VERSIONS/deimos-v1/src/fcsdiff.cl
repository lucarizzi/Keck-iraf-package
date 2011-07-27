#-----------------------------------------------------------------------
procedure fcsdiff( image1, image2 )
#-----------------------------------------------------------------------
# deimos.fcsdiff
# 
# Purpose:
#       Take two fcs images, difference the images, and display them.
#       If images are not specified, it will take the current fcs image
#		and the current fcs reference image
#
# Usage:
#       fcsdiff
# 
# Arguments:
#       image1	-	default is reference image
#		image2	-   default is the last fcs image
# 
# Restrictions:
# 
# Modification history:
#       2006-Sep-27     MK GW     Original version
#-----------------------------------------------------------------------

string	image1="reference"	{ prompt="Name of DEIMOS FCS reference image in use ('reference' to use the current reference image)" }
string	image2="last"	{ prompt="Name of a DEIMOS FCS image ('last' for most recent)" }
bool   display=yes    { prompt="Option to display the image"}
string outfile=""     { prompt="Output file name"}

begin
  string i_img1
  string i_img2
  string cmd
  string host="polo"
  string outdir
  string outdumb
  bool   outdel=no

  #promt the user for input if it was not provided.
i_img1 = image1
i_img2 = image2


# get the current fcs outdir...
  cmd =  "show -s deifcs -terse outdir"
  rsh(host, cmd) | scan (outdir)  
  if ( nscan() != 1 ) 
		error (1,"Could not identify an FCS outdir")

# translate the name of the first image...
# if it is the default use the fcs reference image
	if( i_img2 == "last" ) {
      # get disk name of latest file...
      cmd =  "\ls -1t "//outdir//" | awk '/\.fits$/{print;exit 1}'"
      rsh(host, cmd) | scan (i_img2)  
      if ( nscan() != 1 ) 
		error (1,"Could not find an FCS image in the output directory")
    }

# translate the name of the second image...
# if it is the default use the 'last' available fcs image.
	if( i_img1 == "reference" ) {
      cmd =  "show -s deifcs -terse fcsimgfi"
      rsh(host, cmd) | scan (i_img1)  
      if ( nscan() != 1 ) 
		error (1,"Unable to identify an FCS reference image.")
    }

# Set the output file
# detect whether an output file is specified.
	if(outfile == "") {
		outdumb = mktemp("fcsdiff")//".fits"
		outdel = yes
	} else {  
		outdumb = outdir
	}

if( !imaccess(i_img1)) {
  print("Could not find "//i_img1)
  i_img1 = "/s"//outdir//"/"//i_img1
  if( !imaccess(i_img1)) {
	error(1,"First image does not exist: "//i_img1)
  } else {
	print("Found image "//i_img1)
  }
}

if( !imaccess(i_img2)) {
  print("Could not find "//i_img2)
  i_img2 = "/s"//outdir//"/"//i_img2
  if( !imaccess(i_img2)) {
	error(1,"Second image does not exist: "//i_img2)
  } else {
	print("Found image "//i_img2)
  }
}  

if( imaccess(outdumb)) {
  print("The file "//outdumb//" already exists")
  imdel(outdumb, verify+)
  if( imaccess(outdumb)) 
    error(1, "Outfile exists. Would not over write the output file.")
}

imarith(i_img1,"-",i_img2,outdumb)

if( display ) {
	print("Displaying image "//outdumb)
	display(outdumb, frame=1)
}

if(outdel) {
	imdel(outdumb, verify- )
    print("Deleted temporary file "//outdumb)
}

end


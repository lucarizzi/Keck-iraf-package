# DISTEVAL:  evaluate distortion models -- currently Sutin's DEIMOS distortion
# map
#

include	<math.h>
include	<error.h>


procedure	t_disteval()

char	input[SZ_FNAME]			# input file name
pointer	fda				# in file descriptors

char	tchar
double	thetamin, radinch
double	foclen
pointer	bufx, bufy

int	i, ndx, npts

int	fscan(), nscan()
real	clgetr()
pointer	open()

begin
	call clgstr ("input", input, SZ_FNAME)
	
	foclen = clgetr ("dfoc_len") + 15.0D0

        fda = open (input, READ_ONLY, TEXT_FILE)

# Count the entries
	ndx = 0
	while (fscan (fda) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#' || nscan() == 0) {
			next
		}
		call reset_scan()
		ndx = ndx + 1
	}
	npts = ndx * 2
	call seek (fda, BOF)

# Allocate arrays
	call malloc (bufx, npts, TY_DOUBLE)
	call malloc (bufy, npts, TY_DOUBLE)

	ndx = 0
	while (fscan(fda) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#' || nscan() == 0) {
			next
		}
		call reset_scan()
		call gargd (thetamin)
		call gargd (radinch)
		if (nscan() < 2)
		    call fatal (0, "poorly formatted")

# Convert to radians
		Memd[bufx+ndx] = DEGTORAD(thetamin / 60.)
		Memd[bufy+ndx] = atan (radinch / foclen)
		ndx = ndx + 1
	}

	call adivd (Memd[bufy], Memd[bufx], Memd[bufy], ndx)
	call amulkd (Memd[bufx], -1.D0, Memd[bufx+ndx], ndx)
	call amovd (Memd[bufy], Memd[bufy+ndx], ndx)

	do i = 0, ndx*2-1 {
		call printf (" %10.7f %10.7f\n")
			call pargd (Memd[bufx+i])
			call pargd (Memd[bufy+i])
	}

	call mfree (bufy, TY_DOUBLE)
	call mfree (bufx, TY_DOUBLE)
end



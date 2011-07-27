include	<math.h>

DEFINE	MMPI		25.4		# mm per inch
DEFINE	MAX_Y		330.		# Max y (mm) for slits
DEFINE	DOT_SIZE	0.35		# EXPAND size for dots
DEFINE	SZ_ID		20		# Max char size of ID

define	CODE_AS		-2		# Code for alignment stars

#
# GEN_IGI: produce an input file for igi, to produce truescale drawing
# of the slit-mask. Uses autoslit output as input.
# NB use stsdas.stplot.igi:
#	> igi dev=stdplot ; gflush

procedure t_gen_igi()

char	input[SZ_FNAME]			# Autoslit "setup" file
char	output[SZ_FNAME]		# file of IGI plot commands
char	title[SZ_LINE]
real	pgwidth, pgheight		# Height, width of paper (landscape)
real	slit_wid			# slit width (mm)
bool	true_scale			# True-scale drawing?
bool	label
pointer	fda, fdb

char	tchar
char	id[SZ_ID]
real	xs, ys, y1, y2, delpa
real	x1, x2, scale
real	x1a, x1b, x2a, x2b, width
real	arg6, arg7, arg9
int	pcode

bool	clgetb()
int	fscan(), nscan()
real	clgetr()
pointer	open()

begin
	call clgstr ("input", input, SZ_FNAME)
	fda = open (input, READ_ONLY, TEXT_FILE)
	call clgstr ("output", output, SZ_FNAME)
	fdb = open (output, NEW_FILE, TEXT_FILE)

	pgwidth = clgetr ("width_page")
	pgheight = clgetr ("height_page")
	slit_wid = clgetr ("slit_width")

	call clgstr ("title", title, SZ_LINE)

	true_scale = clgetb ("true_scale")
	label = clgetb ("label")

	if (true_scale)
		scale = 1.
	else
		scale = MAX_Y / (pgwidth * MMPI)

# Setup the initial IGI commands:
	call fprintf (fdb, "physical %f %f %f %f\n")
		call pargr (0.)
		call pargr (pgwidth)
		call pargr (0.)
		call pargr (pgheight)

	call fprintf (fdb, "vpage 0 1 0 1\n")
	call fprintf (fdb, "location 0 1 0 1\n")

	call fprintf (fdb, "fontset igi\n")		# Otherwise "gio"

	call fprintf (fdb, "limits 0 %f 0 %f\n")	# Full page for labels
		call pargr (pgwidth * MMPI)
		call pargr (pgheight * MMPI)
	call fprintf (fdb, "relocate %f %f\n")
		call pargr (pgwidth * MMPI - 5.)
		call pargr (pgheight * MMPI - 5.)
	call fprintf (fdb, "angle -180.\n")
	call fprintf (fdb, "putlabel 9 %s: %s\n")
		call pargstr (input, SZ_FNAME)
		call pargstr (title, SZ_LINE)
	call fprintf (fdb, "angle 0.\n")

	call fprintf (fdb, "limits 0 %f 0 %f\n")	# Scaled page
		call pargr (pgwidth * MMPI * scale)
		call pargr (pgheight * MMPI * scale)

	call fprintf (fdb, "expand %f\n")
		call pargr (DOT_SIZE)
	call fprintf (fdb, "angle 90.\n")
	call fprintf (fdb, "ptype 10 3\n")

# Set the labels in "hardware font mode"
	call fprintf (fdb, "fontset gio\n")		# Otherwise "igi"

# Get the input entries
	while (fscan (fda) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#' || nscan() == 0) {
			next
		}
		call reset_scan()
		call gargr (xs)
		call gargr (ys)
		call gargr (y1)
		call gargr (y2)
		call gargr (delpa)
		call gargr (arg6)
		call gargr (arg7)
		call gargi (pcode)
		if (nscan() < 8) {
			call eprintf ("WARNING: input line skipped\n")
			next
		}
		if (delpa == INDEF)
			delpa = 0.
		if (delpa < -90.) {
			delpa = delpa + 180.
		} else if (delpa > 90.) {
			delpa = delpa - 180.
		}
		if (label) {
			call gargr (arg9)
			call gargwrd (id, SZ_ID)
			if (nscan() < 10)
				id[1] = EOS
		}

#	For each slit, write the IGI commands
		if (pcode == CODE_AS) {
			width = y1 - y2
			x1 = xs
			x2 = xs
		} else {
			width = slit_wid
			x1 = xs + (y1 - ys) * tan (DEGTORAD(delpa))
			x2 = xs + (y2 - ys) * tan (DEGTORAD(delpa))
		}
		x1a = x1 - width / 2.
		x1b = x1 + width / 2.
		x2a = x2 - width / 2.
		x2b = x2 + width / 2.
		call fprintf (fdb, "rel  %f %f\n")
			call pargr (y1)
			call pargr (x1a)
		call fprintf (fdb, "draw %f %f\n")
			call pargr (y2)
			call pargr (x2a)
		call fprintf (fdb, "draw %f %f\n")
			call pargr (y2)
			call pargr (x2b)
		call fprintf (fdb, "draw %f %f\n")
			call pargr (y1)
			call pargr (x1b)
		call fprintf (fdb, "draw %f %f\n")
			call pargr (y1)
			call pargr (x1a)
		call fprintf (fdb, "rel  %f %f\n")
			call pargr (ys)
			call pargr (xs)
		call fprintf (fdb, "dot\n")
		if (label) {
			call fprintf (fdb, "rel  %f %f\n")
				call pargr (ys)
				call pargr (xs + 1.5)
			call fprintf (fdb, "expand %f\n")
				call pargr (1.6*DOT_SIZE)
			call fprintf (fdb, "putlabel 6 %s\n")
				call pargstr (id, SZ_ID)
			call fprintf (fdb, "expand %f\n")
				call pargr (DOT_SIZE)
		}
	}

# If true-scale, print metric lines
	if (true_scale) {
		call fprintf (fdb, "fontset igi\n")	# Otherwise "gio"
#	X-scale metric
		call fprintf (fdb, "expand 0.8 \n")
		call fprintf (fdb, "rel   3.5   5.\n")
		call fprintf (fdb, "draw  6.5   5.\n")
		call fprintf (fdb, "rel   3.5 185.\n")
		call fprintf (fdb, "draw  6.5 185.\n")
		call fprintf (fdb, "rel   5.    5.\n")
		call fprintf (fdb, "draw  5.  185.\n")
		call fprintf (fdb, "rel   3.  95.\n")
		call fprintf (fdb, "angle 90.\n")
		call fprintf (fdb, "putlabel 5 180 mm\n")

#	Y-scale metric
		call fprintf (fdb, "rel   10. 3.5\n")
		call fprintf (fdb, "draw  10. 6.5\n")
		call fprintf (fdb, "rel  250. 3.5\n")
		call fprintf (fdb, "draw 250. 6.5\n")
		call fprintf (fdb, "rel   10. 5.\n")
		call fprintf (fdb, "draw 250. 5.\n")
		call fprintf (fdb, "rel  130. 3.\n")
		call fprintf (fdb, "angle 180.\n")
		call fprintf (fdb, "putlabel 5 240 mm\n")
		call fprintf (fdb, "expand 1  \n")

#	Outer border
		call fprintf (fdb, "rel  -1000. -1000.\n")
		call fprintf (fdb, "draw  1000. -1000.\n")
		call fprintf (fdb, "draw  1000.  1000.\n")
		call fprintf (fdb, "draw -1000.  1000.\n")
		call fprintf (fdb, "draw -1000. -1000.\n")
	}


# End the IGI commands
	call fprintf (fdb, "end\n")

	call close (fdb)
	call close (fda)
end


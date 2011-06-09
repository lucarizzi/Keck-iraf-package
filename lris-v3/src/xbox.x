include <imhdr.h>
include <gset.h>
include <math.h>
include "futil.h"

define          REL_HT  0.4             # relative height for crossings
define          ID_CHSZ 9               # Character size of ID string
define          LAB_SZ 80               # Char size of title/label strings
define          X_HWID  1.5             # x half-width of triangle function
define          Y_HWID  1.5             # y half-width of triangle function
define          X_CRAD  13              # x centering radius
define          Y_CRAD  13              # y centering radius
define          FRAC    0.333           # fractional lev. in sorted list for sky

define          PRE_COL 21              # Number of prescan columns per amp

# Neither used:
# define                KEYSFILE1       "lris$lib/scr/boxfind.key"
# define                KEYSFILE2       "lris$lib/scr/maskalign.key"

define          ASECPIXR  0.2109        # Arcsec/pixel [LRIS red side]
define          ASECPIXB  0.2151        # Arcsec/pixel [LRIS blue side]
define          PA_ACCEPT 0.027         # this is 0.1" at 1000 px
define          SZ_FITS   40            # number of characters in keyword

#-----------------------------------------------------------------------
# Name:
#       xbox
# 
# Purpose:
#       locate boxes on LRIS slitmasks and measure centroids of stars 
#       within them, and solve for the instrument rotation and telescope
#       translation required to move the stars to the box centers.
#       This program is a combination of two previous tasks, mboxfind
#       and maskalign.
#
#       Since the boxes are relatively small, read in section and deal 
#       with the 2-D image rather than cuts.
#
#       We currently have a kludge for robustness -- runs box-center twice,
#       the second time after recentering image section on initial box center.
# 
# Author:
#       A. C. Phillips (UCO/Lick)
#
# Modification history:
#       Date unknown    ACP     Original version (0.0)
#       2000-Oct-24     GDW     Modified for blue side (v1.0):
#                               - added logic to determine pixel scale
#                                 from INSTRUME keyword
#                               - changed 'Send DCS commands' prompt to
#                                 ensure that it happens even when the 
#                                 task in invoked via ":go"
#       2001-Sep-09     GDW     Modified for use with broken DCS:
#                               - added prompt to have OA disable
#                               guiding before rotating.
#                               - added prompt to have OA enable
#                               guiding after rotating.
#-----------------------------------------------------------------------

procedure t_xbox()

char    image[SZ_FNAME]                 # Input Image 
char    input[SZ_FNAME]                 # Input Box coordinates
char    coordfile[SZ_FNAME]             # Box-Star coordinates (pre-measured)
pointer im
int     nxbox, nybox                    # Size of subraster
real    xsz, ysz                        # Full size of boxes (pix)
real    xfwhm, yfwhm                    # FWHM of star (pix)
real    xoff, yoff                      # x,y shifts to apply to input box coord
pointer fda, fdb
bool    find_coord                      # find box/star coordinate pairs?

char    tchar, idstr[ID_CHSZ]
int     ncols, nlines
int     npt, ndx, i, j
int     istat
bool    stat

int     namp, prepix
int     nx, ny
real    xs, ys
real    xb, yb, xstar, ystar
int     x1, y1, x2, y2
pointer bufx, bufy, bufzx, bufzy, buftx, bufty
pointer buf

real    rot, el, xflex, yflex           # parameters for flexure correction
real    asecpix                         # pixel scale arcsec/px

# parameters for aligning
real    xrot, yrot
bool    invert
# int   nx, ny
int     niter
real    bxs, def1err, def2err
bool    dcs, dcspa

char    poname[SZ_FITS]                         # Pointing Origin name
char    instrument[SZ_FITS]                     # instrument [LRIS|LRISBLUE]
char    cmdline[SZ_LINE]                        # Command line for rsh
real    rotposn                                 # PA = rotposn+90
real    coeff[2,3], err[3]
real    sina, cosa
real    earcsec, narcsec, eerr, nerr, pdeg, perr
pointer xbuf1, ybuf1, xbuf2, ybuf2
pointer ebufx1, ebufy1, ebufx2, ebufy2
pointer rbufx, rbufy, wbuf

bool    box_graph()
bool    imaccf(), streq(), strne(), yesno()
int     clgeti(), imgeti()
int     fscan(), nscan(), scan()
int     oscmd()
real    clgetr(), imgetr()
pointer immap(), open(), imgs2r()

begin
        call clgstr ("image", image, SZ_FNAME)
        nxbox = clgeti ("nxbox")
        nybox = clgeti ("nybox")
        xsz = clgetr ("xsz")
        ysz = clgetr ("ysz")
        xfwhm = clgetr ("fwhm")
        yfwhm = clgetr ("fwhm")
        xoff  = clgetr ("xoff")
        yoff  = clgetr ("yoff")

        im = immap (image, READ_ONLY, 0)
        ncols = IM_LEN(im,1)
        nlines = IM_LEN(im,2)

        call clgstr ("pairs", coordfile, SZ_FNAME)
        if (streq (coordfile, "")) {
                find_coord = true
        } else {
                fdb = open (coordfile, READ_ONLY, TEXT_FILE)
                find_coord = false
        }

# Get relevant keywords ...
# ... for the number of prescan columns     # ?? FIX - need string
        if (imaccf (im, "NUMAMPS")) {
                namp = imgeti (im, "NUMAMPS")
        } else {
                call eprintf ("NUMAMPS missing; ONE AMP assumed!\n")
                namp = 1
        }

        if (imaccf (im, "PREPIX")) {
                prepix = namp * imgetr (im, "PREPIX")
        } else {
                call eprintf ("PREPIX missing; %2d assumed!\n")
                        call pargi (PRE_COL)
                prepix = namp * PRE_COL
        }

# ... for the rotator angle and elevation for flexure mapping
        if (imaccf (im, "ROTPPOSN")) {
                rot = imgetr (im, "ROTPPOSN")
        } else {
                call eprintf ("ROTPPOSN missing; 90 (LRIS stow) assumed!\n")
                rot = 90.
        }

        if (imaccf (im, "EL")) {
                el = imgetr (im, "EL")
        } else {
                call eprintf ("EL missing; 0 (LRIS stow) assumed!\n")
                el = 0.
        }
# get the flexure estimates
        call flex_corr (rot, el, xflex, yflex)

# ... and for the PA and pointing origin
        if (imaccf (im, "ROTPOSN")) {
                rotposn = imgetr (im, "ROTPOSN")
        } else {
                call eprintf ("ROTPOSN missing; -180 (LRIS stow) assumed\n")
                rotposn = -180.
        }

        if (imaccf (im, "PONAME")) {
                call imgstr (im, "PONAME", poname, SZ_FITS)
        } else {
                call strcpy ("UNKNOWN", poname, SZ_FITS)
        }

# Assign the proper rotation center for the PO
        if (streq (poname, "LRIS")) {
                xrot = 1024.5
                yrot = 1024.5
        } else if (streq (poname, "Pickoff")) {
                xrot = 1327
                yrot = 890
        } else {
                if (strne (poname, ""))
                        call eprintf ("Unknown Pointing Origin!\n")
                xrot = 1024.5
                yrot = 1024.5
        }

# determine which side of the spectrograph is in use and set the
# pixel scale accordingly...
        asecpix = 0.
        if (imaccf (im, "INSTRUME")) {
                call imgstr (im, "INSTRUME", instrument, SZ_FITS)
                if (streq ( instrument, "LRIS")){
                        asecpix = ASECPIXR
                        call eprintf ("Image taken with RED side of LRIS.\n")
                } else if (streq ( instrument, "LRISBLUE")){
                        asecpix = ASECPIXB
                        call eprintf ("Image taken with BLUE side of LRIS.\n")
                }
        }
        if ( asecpix == 0. ) {
                asecpix = ASECPIXR
                call eprintf ("Invalid INSTRUME keyword value in image header!\n")
        }
        call eprintf ("Using pixel scale of %f\n")
                        call pargr( asecpix)

# Print out summary of image:
        call printf ("  Prescan col.   flexure     offset     Total\n")
        call printf ("X:     %3d         %5.1f      %5.1f     %5.1f\n")
                call pargi (prepix)
                call pargr (xflex)
                call pargr (xoff)
                call pargr (prepix+xflex+xoff)
        call printf ("Y:     %3d         %5.1f      %5.1f     %5.1f\n\n")
                call pargi (0)
                call pargr (yflex)
                call pargr (yoff)
                call pargr (yflex+yoff)

        xoff = prepix+xflex+xoff
        yoff = yflex+yoff

# READY: DO WE WANT TO FIND BOX/STAR COORDS?
        if (find_coord) {

# Get and open the input file
            call clgstr ("input", input, SZ_FNAME)
            fda = open (input, READ_ONLY, TEXT_FILE)

# Allocate arrays for marginal plots
            call malloc (bufx, nxbox, TY_REAL)
            call malloc (bufy, nybox, TY_REAL)
            call malloc (bufzx, nxbox, TY_REAL)
            call malloc (bufzy, nybox, TY_REAL)
            call malloc (buftx, nxbox, TY_REAL)
            call malloc (bufty, nybox, TY_REAL)

# Count entries in input file
            npt = 0
            while (fscan(fda) != EOF)
                npt = npt + 1
            call seek (fda, BOF)

# Allocate the coordinate arrays (others deferred until needed)
            call malloc (xbuf1, npt, TY_REAL)           # ref coord
            call malloc (ybuf1, npt, TY_REAL)
            call malloc (xbuf2, npt, TY_REAL)           # input coord
            call malloc (ybuf2, npt, TY_REAL)

# Get the input entries
            ndx = 0
            while (fscan (fda) != EOF) {
                call gargwrd (tchar, 1)
                if (tchar == '#' || nscan() == 0) {
                        next
                }
                call reset_scan()
                call gargr (xs)
                call gargr (ys)
                if (nscan() < 2) {
                        call eprintf ("WARNING: input line skipped\n")
                        next
                }

# adjust values for offsets:
                xs = xs + xoff
                ys = ys + yoff

                call gargwrd (idstr, ID_CHSZ)
                if (nscan() < 3)
                        call strcpy ("(no ID)", idstr, ID_CHSZ)

# Find the Box; we do a kludge -- run this twice, recentered on box the 2nd time
                xb = xs
                yb = ys
                do j = 1, 2 {
                    x1 = xb - nxbox/2
                    x2 = x1 + nxbox - 1
                    y1 = yb - nybox/2
                    y2 = y1 + nybox - 1

# checks on out-of-bounds
                    x1 = max (1, x1)
                    x2 = min (ncols, x2)
                    y1 = max (1, y1)
                    y2 = min (nlines, y2)

# Actual box length 
                    nx = x2 - x1 + 1
                    ny = y2 - y1 + 1

# Get the image section
                    buf = imgs2r (im, x1, x2, y1, y2)

# Fill position vectors
                    do i = 0, nx-1 {
                        Memr[bufx+i] = i + x1
                    }
                    do i = 0, ny-1 {
                        Memr[bufy+i] = i + y1
                    }

# Get the box position
                    call box_center (Memr[bufx], Memr[bufy], Memr[buf], nx, ny,
                        Memr[bufzx], Memr[bufzy], xsz, ysz, xb, yb)
# On second time through loop recenter on box
                }

                call printf ("Box center:  %6.2f %6.2f  (%4.1fx,%4.1fy removed) (del:%4.1f,%4.1f)\n")
                        call pargr (xb-xoff)
                        call pargr (yb-yoff)
                        call pargr (xoff)
                        call pargr (yoff)
                        call pargr (xb-xs)
                        call pargr (yb-ys)

# Now get star location
                stat = box_graph (Memr[bufx], Memr[bufy], Memr[bufzx],
                        Memr[bufzy], Memr[buftx], Memr[bufty], nx, ny, idstr,
                                xb, yb, xsz, ysz, xfwhm, yfwhm, xstar, ystar)

# Store box/star coordinates
# Note -- unsure if better to subtract xoff,yoff (negligible if small)
                if (stat == true) {
                        Memr[xbuf1+ndx] = xstar - prepix
                        Memr[ybuf1+ndx] = ystar
                        Memr[xbuf2+ndx] = xb - prepix
                        Memr[ybuf2+ndx] = yb
                        ndx = ndx + 1
                }
            }
            npt = ndx           # npt reduced by bad pairs, comment lines

# Done with the box-finding: clean up
            call mfree (bufty, TY_REAL)
            call mfree (buftx, TY_REAL)
            call mfree (bufzy, TY_REAL)
            call mfree (bufzx, TY_REAL)
            call mfree (bufy, TY_REAL)
            call mfree (bufx, TY_REAL)
            call close (fda)
            call imunmap (im)

        } else {
# ... OR READ IN COORDINATE PAIRS?

# Count entries in input file
            npt = 0
            while (fscan(fdb) != EOF)
                npt = npt + 1
            call seek (fdb, BOF)

# Allocate the coordinate arrays (others deferred until needed)
            call malloc (xbuf1, npt, TY_REAL)           # ref coord
            call malloc (ybuf1, npt, TY_REAL)
            call malloc (xbuf2, npt, TY_REAL)           # input coord
            call malloc (ybuf2, npt, TY_REAL)

# Get the input entries
            ndx = 0
            while (fscan (fdb) != EOF) {
                call gargwrd (tchar, 1)
                if (tchar == '#' || nscan() == 0) {
                        next
                }
                call reset_scan()
                call gargr (xb)
                call gargr (yb)
                call gargr (xstar)
                call gargr (ystar)
                if (nscan() < 4) {
                        call eprintf ("WARNING: input line skipped\n")
                        next
                }
                Memr[xbuf1+ndx] = xstar - prepix
                Memr[ybuf1+ndx] = ystar
                Memr[xbuf2+ndx] = xb - prepix
                Memr[ybuf2+ndx] = yb
                ndx = ndx + 1
            }
            npt = ndx           # npt reduced by bad pairs, comment lines
            call close (fdb)
            call imunmap (im)
        }

# ON TO ALIGNING: setup
        nx = ncols                      # TMP: should fix; also symmetrize
        ny = nlines                     # TMP   ... maskalign
        bxs = -2. * clgetr ("box_size")
        niter = clgeti ("niter")
        def1err = clgetr ("def_ref_err")
        def2err = clgetr ("def_err")

# Allocate the remaining arrays needed
        call malloc (ebufx1, npt, TY_REAL)              # error arrays
        call malloc (ebufy1, npt, TY_REAL)
        call malloc (ebufx2, npt, TY_REAL)
        call malloc (ebufy2, npt, TY_REAL)
        call malloc (rbufx, npt, TY_REAL)               # residuals
        call malloc (rbufy, npt, TY_REAL)
        call malloc (wbuf, npt, TY_REAL)                # weights

        call amovkr (1., Memr[wbuf], npt)

# Store default error estimates (someday could be individual...)
        call amovkr (def1err, Memr[ebufx1], npt)
        call amovkr (def1err, Memr[ebufy1], npt)
        call amovkr (def2err, Memr[ebufx2], npt)
        call amovkr (def2err, Memr[ebufy2], npt)

# square the errors, as required for get_lsqf:
        call amulr (Memr[ebufx1], Memr[ebufx1], Memr[ebufx1], npt)
        call amulr (Memr[ebufy1], Memr[ebufy1], Memr[ebufy1], npt)
        call amulr (Memr[ebufx2], Memr[ebufx2], Memr[ebufx2], npt)
        call amulr (Memr[ebufy2], Memr[ebufy2], Memr[ebufy2], npt)

# Reference coordinates to center of rotation axis
        call aaddkr (Memr[xbuf1], -xrot, Memr[xbuf1], npt)
        call aaddkr (Memr[ybuf1], -yrot, Memr[ybuf1], npt)
        call aaddkr (Memr[xbuf2], -xrot, Memr[xbuf2], npt)
        call aaddkr (Memr[ybuf2], -yrot, Memr[ybuf2], npt)

        call maskalign (Memr[xbuf1], Memr[ybuf1], Memr[xbuf2], Memr[ybuf2],
                Memr[ebufx1], Memr[ebufy1], Memr[ebufx2], Memr[ebufy2],
                Memr[rbufx], Memr[rbufy], Memr[wbuf], npt, nx, ny, niter,
                                                xrot, yrot, bxs, coeff, err)
# Print out results
        call printf ("#  obj-x  obj-y  slit-x   xres   slit-y   yres     w \n")
        do i = 0, npt-1 {
                call printf ("#%7.2f%7.2f  %7.2f (%5.3f) %7.2f (%5.3f) %6.2f\n")
                        call pargr(Memr[xbuf1+i])
                        call pargr(Memr[ybuf1+i])
                        call pargr(Memr[xbuf2+i])
                        call pargr(Memr[rbufx+i])
                        call pargr(Memr[ybuf2+i])
                        call pargr(Memr[rbufy+i])
                        call pargr(Memr[wbuf+i])
        }

        call printf ("\n\n#  x-xform: %7.5fx + %7.5fy  + %7.3f (%5.3f)\n")
                call pargr (coeff[1,1])
                call pargr (coeff[1,2])
                call pargr (coeff[1,3])
                call pargr (err[2])
        call printf ("#  y-xform: %7.5fx + %7.5fy  + %7.3f (%5.3f)\n")
                call pargr (coeff[2,1])
                call pargr (coeff[2,2])
                call pargr (coeff[2,3])
                call pargr (err[3])

# Now adjust for existing PA (convert x,y into E,N):
# Note that "angle" = (angle of "PA" wrt x-axis) - (PA = rotposn+90)
#        ...        = 90 - (rotposn+90)
        invert = false                     # Left for historical reasons
        if (invert) {
                call eprintf ("Sorry -- not yet tested\n")
                pdeg = RADTODEG( atan2(coeff[1,2], coeff[1,1]))
                perr = RADTODEG(asin(err[1]))
                cosa = cos (DEGTORAD(-rotposn))
                sina = sin (DEGTORAD(-rotposn))
                narcsec = asecpix * (-cosa * coeff[1,3] + sina * coeff[2,3])
                earcsec = asecpix * ( sina * coeff[1,3] + cosa * coeff[2,3])
                nerr = asecpix * sqrt ((cosa * err[2])**2 + (sina * err[3])**2)
                eerr = asecpix * sqrt ((sina * err[2])**2 + (cosa * err[3])**2)
        } else {
                pdeg = RADTODEG(-atan2(coeff[1,2], coeff[1,1]))
                perr = RADTODEG(asin(err[1]))
                cosa = cos (DEGTORAD(-rotposn))
                sina = sin (DEGTORAD(-rotposn))
                narcsec = asecpix * ( cosa * coeff[1,3] + sina * coeff[2,3])
                earcsec = asecpix * (-sina * coeff[1,3] + cosa * coeff[2,3])
                nerr = asecpix * sqrt ((cosa * err[2])**2 + (sina * err[3])**2)
                eerr = asecpix * sqrt ((sina * err[2])**2 + (cosa * err[3])**2)
        }

# PA offset is backwards from rotation:
        call printf ("\n====================================================\n")
        Call printf ("\n*** MOVE TELESCOPE/ROTATOR by the following offsets:\n")
        call printf ("\n Offset PA by %5.2f (%4.2f) degree\n")
                call pargr (pdeg)
                call pargr (perr)
        call printf (" Offsets: %6.2f\" EAST (%4.2f)\t %6.2f\" NORTH (%4.2f)\n\n")
                call pargr (earcsec)
                call pargr (eerr)
                call pargr (narcsec)
                call pargr (nerr)
        call printf ("====================================================\n\n")

        dcs = yesno( "[CONFIRM]  Send DCS commands?", SZ_LINE, true)
        if (dcs) {
# Check rotator error for magnitude
            if (abs (pdeg) < PA_ACCEPT) {
                call printf ("PA error is small...  ")
                dcspa = yesno( "[CONFIRM]  Apply rotator update?", SZ_LINE, false)
            } else {
                dcspa = true
            }
# Apply rotation
            if (dcspa) {

                # verify that guiding is OFF...
                call eprintf( "\n\007")
                call eprintf( "------------------------------------------------------------------------------\n")
                call eprintf( "Have the OA STOP guiding while we rotate LRIS, then press ENTER to continue...\n")
                call eprintf( "------------------------------------------------------------------------------\n")
                istat = scan()

                # rotate...
                call sprintf (cmdline, SZ_LINE,
                  "%s modify -s dcs ROTDEST=%.2f ROTMODE=1")    # (1=pos.angle)
                        call pargstr ("rsh punaluu")
                        call pargr (rotposn+pdeg)
                call eprintf ("\n(sending)  %s \n")
                        call pargstr (cmdline)
                istat = oscmd (cmdline)
                if (istat != OK) {
                        call eprintf ("command failed!  (%d)\n")
                                call pargi (istat)
                }

                call sprintf (cmdline, SZ_LINE, "%s waitfor -s dcs ROTSTAT=8")
                        call pargstr ("rsh punaluu")    # (8=tracking)
                call eprintf ("\n(sending)  %s \n")
                        call pargstr (cmdline)
                istat = oscmd (cmdline)
                if (istat != OK) {
                        call eprintf ("command failed!  (%d)\n")
                                call pargi (istat)
                }

                # verify that guiding is back ON...
                call eprintf( "\n\007")
                call eprintf( "------------------------------------------------------------------------------\n")
                call eprintf( "Have the OA RESUME guiding, then press ENTER to continue...\n")
                call eprintf( "------------------------------------------------------------------------------\n")
                istat = scan()
            }

# ... and apply translation
                call sprintf (cmdline, SZ_LINE,
                  "%s modify -s dcs RAOFF=%.2f DECOFF=%.2f REL2CURR=1")
                        call pargstr ("rsh punaluu")
                        call pargr (earcsec)
                        call pargr (narcsec)
                call eprintf ("\n(sending)  %s \n")
                        call pargstr (cmdline)
                istat = oscmd (cmdline)
                if (istat != OK) {
                        call eprintf ("command failed!  (%d)\n")
                                call pargi (istat)
                }

                call sprintf (cmdline, SZ_LINE, "%s waitfor -s dcs AXESTAT=64")
                        call pargstr ("rsh punaluu")    # (64=tracking)
                call eprintf ("\n(sending)  %s \n")
                        call pargstr (cmdline)
                istat = oscmd (cmdline)
                if (istat != OK) {
                        call eprintf ("command failed!  (%d)\n")
                                call pargi (istat)
                }

                call eprintf ("... done! \007 \n")
        }

end


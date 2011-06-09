include <math.h>
include "lris.h"
#-----------------------------------------------------------------------
# Name:
#       mshiftb
# 
# Purpose:
#       Execute a telescope offset in pixel coordinates, given the
#       current coordinates, desired (i.e., "reference") coordinates,
#       current LRIS position angle, and LRIS camera (red or blue)
# 
# Author:
#       A. C. Phillips (UCO/Lick)
#
# Modification history:
#       Date unknown    ACP     Original version (0.0)
#       2000-Oct-24     GDW     Modified for blue side (v1.0):
#                               - added 'beam' parameter for red or blue
#                               - added logic to determine pixel scale
#                                 from beam name
#                               - changed 'Send DCS commands' prompt to
#                                 ensure that it happens even when the 
#                                 task in invoked via ":go"
#       2002-June-01    AES     Modified for new blue camera (v2.0)
#                               - added 'xgap' parameter which is x-value
#                                 at start of RH blue chip
#                               - added new blue pixel scale    
#                               - rotated by 90 degrees for blue side
#                               - added in coord transformations to go
#                                 from detector to pixel coords for RH blue chip
#-----------------------------------------------------------------------

procedure t_mshiftb ()

real    x, y                                    # measured x and y (pixels)
real    xref, yref                              # desired x and y (pixels)
real    pa                                      # pa in degrees
real    narcsec, earcsec                        # returned offsets, N and E
real    pscale                                  # pixel scale [arcsec/px]
bool    invert                                  # invert orient. (x-flipped)?
bool    dcs                                     # print DCS commands
char    beam[SZ_LINE]                           # name of current LRIS beam

char    cmdline[SZ_LINE]                        # command string for DCS
int     stat                                    # command stat
bool    answer                                  # user reply

bool    clgetb(), streq(), yesno()
int     oscmd()
real    clgetr()

begin
        x = clgetr ("xobs")
        y = clgetr ("yobs")
        xref = clgetr ("xref")
        yref = clgetr ("yref")
        pa = clgetr ("pa")
        invert = clgetb ("invert")
        call clgstr ("beam", beam, SZ_LINE)
        dcs = clgetb ("dcs")

        # assign the appropriate pixel scale...
        if (streq (beam, "red")) {
                pscale = ASECPPR
        } else if (streq (beam, "blue")) {
                pscale = ASECPPNB
        } else {
                call eprintf ("Invalid BEAM; RED side assumed!\n")
                pscale = ASECPPR
        }
        call printf( "\nUsing pixel scale of %f arcsec/pixel.\n\n")
                call pargr( pscale)

        call mshiftb (x, y, xref, yref, pa, pscale, narcsec, earcsec, invert)

        call printf ("====================================================\n\n")
        call printf (
                " To shift (%5.1f,%5.1f) to (%5.1f,%5.1f) at mask PA=%5.1f:\n")
                call pargr (x)
                call pargr (y)
                call pargr (xref)
                call pargr (yref)
                call pargr (pa)
        if (invert)
                call printf (" (one-amp, x-inverted mode)\n")
        call printf ("\n   MOVE TELECSOPE   %5.2f''E  and %5.2f''N \n\n")
                call pargr (earcsec)
                call pargr (narcsec)
        call printf ("====================================================\n\n")

        if (dcs) {

                # verify...
                answer = yesno( "[CONFIRM]  Send DCS commands?", SZ_LINE, true)
                if ( ! answer ) { call error( 1, "Aborted.") }

                call sprintf (cmdline, SZ_LINE,
                  "%s modify -s dcs RAOFF=%.2f DECOFF=%.2f REL2CURR=1")
                        call pargstr ("rsh punaluu")
                        call pargr (earcsec)
                        call pargr (narcsec)
                call eprintf ("\n %s \n")
                        call pargstr (cmdline)
                stat = oscmd (cmdline)
                if (stat != OK) {
                        call eprintf ("command failed!  (%d)\n")
                                call pargi (stat)
                }

                call sprintf (cmdline, SZ_LINE, "%s waitfor -s dcs AXESTAT=64")
                        call pargstr ("rsh punaluu")    # (64=tracking)
                call eprintf ("\n %s \n")
                        call pargstr (cmdline)
                stat = oscmd (cmdline)
                if (stat != OK) {
                        call eprintf ("command failed!  (%d)\n")
                                call pargi (stat)
                }

                call eprintf ("... done! \n")
        }
end

#
# MSHIFT: work out pixel shifts for a given MASK PA
#

procedure mshiftb (x, y, xref, yref, pa, pixscale, narcsec, earcsec, invert)

real    x, y                                    # measured x and y (pixels)
real    xref, yref                              # desired x and y (pixels)
real    pa                                      # pa in degrees
real    pixscale                                # arcsec/pixel
real    narcsec, earcsec                        # returned offsets, N and E
bool    invert                                  # invert orient. (x-flipped)?

real    cosa, sina
real    xgap
real    clgetr()

begin
        cosa = cos (DEGTORAD(90.-(pa+90.)))
        sina = sin (DEGTORAD(90.-(pa+90.)))
	xgap=clgetr("xgap")

	if (x>xgap) {
	    x=x+96.627+0.00162*y
	    y=y+2.982
	}
	if (xref>xgap) {
	    xref=xref+96.627+0.00162*yref
	    yref=yref+2.982
	}
        if (!invert) {
                narcsec = pixscale * ( cosa * (x-xref) + sina * (y-yref))
                earcsec = pixscale * (-sina * (x-xref) + cosa * (y-yref))
        } else {
                narcsec = pixscale * (-cosa * (x-xref) + sina * (y-yref))
                earcsec = pixscale * ( sina * (x-xref) + cosa * (y-yref))
        }
end


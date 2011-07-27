#-----------------------------------------------------------------------
# Name:
#       lris.h
# 
# Purpose:
#       Define constants used for mshift, lris_util, and simulator tasks
# 
# Author:
#       A. C. Phillips (UCO/Lick)
#
# Modification history:
#       Date unknown    ACP     Original version (0.0)
#       2000-Oct-24     GDW     Modified for blue side (v1.0):
#                               - updated ASECPPX parameter
#                               - added ASECPPR and ASECPPB parameters
#-----------------------------------------------------------------------
define          CCD_NX  2048            # X-Size of CCD (pix)
define          CCD_NY  2048            # X-Size of CCD (pix)
define          CCD_X0  1024            # X-Center of CCD (pix)
define          CCD_Y0  1024            # Y-Center of CCD (pix)
define          GUID_X1  20             # Lower X-limit to guider box (mm)
define          GUID_X2 196             # Upper X-limit to guider box (mm)
define          GUID_Y1 386             # Lower Y-limit to guider box (mm)
define          GUID_Y2 416             # Upper Y-limit to guider box (mm)
define          ASECPMM 1.378           # Arcsec/mm
define          ASECPPX 0.2109          # Arcsec/pix
define          ASECPPR 0.2109          # Arcsec/pix
define          ASECPPB 0.2151          # Arcsec/pix
define          ASECPPNB 0.1353          # Arcsec/pix
define          OBS_LAT 19.8            # Keck latitude
define          CAM_ANG 44.08           # Angle of camera wrt collimator
define          GEN_ZMX 75.0            # General max zenith angle
define          SRV_AZ1 185.0           # Azimuth of Service tower exclusion K2
define          SRV_AZ2 332.0           # Azimuth of Service tower exclusion K2
define          SRV_ZMX 53.2            # Z-angle of Service tower exclusion K2
# define        SRV_AZ1 1.6             # Azimuth of Service tower exclusion
# define        SRV_AZ2 151.0           # Azimuth of Service tower exclusion
# define        SRV_ZMX 54.             # Z-angle of Service tower exclusion

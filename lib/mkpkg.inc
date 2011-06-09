# Global MKPKG definitions for the KECK package.

$set XFLAGS     = "$(XFLAGS) -p keck"
$set XVFLAGS    = "$(XVFLAGS) -p keck"
$set LFLAGS     = "$(LFLAGS) -p keck"

# Uncomment and modify the following to add special file list entries for
# various machine architectures and Fortran compilers.

#$ifeq (MACH, ssun) then
#    $include "keck$lib/mkpkg.sf.SSUN"
#$else $ifeq (MACH, redhat) then
#    $include "keck$lib/mkpkg.sf.REDHAT"
#$else $ifeq (MACH, sparc) then
#    $include "keck$lib/mkpkg.sf.SUN4"
#$end

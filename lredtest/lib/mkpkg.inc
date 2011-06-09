# I don't know what this is all about....
# Global MKPKG definitions for the LREDTEST package.

$set XFLAGS	= "$(XFLAGS) -p lredtest"
$set XVFLAGS	= "$(XVFLAGS) -p lredtest"
$set LFLAGS	= "$(LFLAGS) -p lredtest"

# Uncomment and modify the following to add special file list entries for
# various machine architectures and Fortran compilers.

# $ifeq (MACH, sparc) then
# $include "lredtest$lib/mkpkg.sf.sun4"
# $else $ifeq (MACH, vms) then
# $include "lredtest$lib/mkpkg.sf.vms"
# $end

#{ KECK - The W. M. Keck Observatory suite of packages for IRAF
#
# Author:
#   Gregory D. Wirth [wirth@keck.hawaii.edu]
#
#-----------------------------------------------------------------------

# load global defs...
cl < "keck$lib/zzsetenv.def"

package keck, bin = keckbin$

	# define subpackages...
	task	$deimos.pkg	= "deimos$deimos.cl"
##	task	$dproto.pkg	= "dproto$dproto.cl"
	task	$lris.pkg	= "lris$lris.cl"
	task	$lris2.pkg	= "lris2$lris2.cl"
	task	$lred.pkg	= "lred$lred.cl"
	task	$lredtest.pkg	= "lredtest$lredtest.cl"
	task	$lws.pkg	= "lws$lws.cl"
	task	$nirspec.pkg	= "nirspec$nirspec.cl"

	# common tasks...
	task breakname		= "keck$src/breakname.cl"

	hidetask breakname

	type "keck$motd"

clbye()

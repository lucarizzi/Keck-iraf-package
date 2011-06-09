#{ Package LRED is Drew Phillips' lris software, based on DEIMOS, for LRED


print ("\n LREDTEST (ver.0a for IRAF 2.14.1) -- Software in development -- User assumes risk\n")

cl < "lredtest$lib/zzsetenv.def"
package	lredtest, bin = lredtestbin$

task	lbox = "lredtest$src/x_lredtest.e"


clbye()

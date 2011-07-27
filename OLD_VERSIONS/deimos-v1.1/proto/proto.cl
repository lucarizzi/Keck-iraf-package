#{ Package DEIMOS is Drew Phillips' deimos software


cl < "proto$lib/zzsetenv.def"
package	proto, bin = protobin$

task	dsimulator,
	simtest,
	refl,
	trace,
	msolve,
	qextr,
	qtrmap,
	qmodel,
	qrevmod,
	disteval,
	mkmosaic,
	mkmontage,
	ccpkfind,
	vern,
	rtest,
	$xlate,
	$tvalign,
	$tmos,
	xbox,
	$get_boxes,
	decode,
	dastrom,
	cohu	= "proto$src/x_proto.e"

task dss_chart="proto$src/dss_chart.cl"


clbye()

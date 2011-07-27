#{ Package DEIMOS is Drew Phillips' deimos software

# load requried packages...
mscred	# needed for dmosdisplay
;

cl < "deimos$lib/zzsetenv.def"
package	deimos, bin = keckbin$

    # script tasks...
	task	$cdata		= "deimos$src/cdata.cl"
	task	dss_chart	= "deimos$src/dss_chart.cl"
    task	$do_xbox	= "deimos$src/do_xbox.cl"
    task	$lastimage	= "$rsh polo lastimage"
    task	$wfi		= "$rsh polo wfi"
    task	get_seeing	= "deimos$src/get_seeing.cl"
    task	tune_qmodel	= "deimos$src/tune_qmodel.cl"
    task 	check_boxes	= "deimos$src/check_boxes.cl"
    task 	dmosdisplay	= "deimos$src/dmosdisplay.cl"
    task 	plotmask	= "deimos$src/plotmask.cl"
    task 	fcsdiff		= "deimos$src/fcsdiff.cl"

task	dsimulator,
	refl,
	trace,
	qtrmap,
	qmodel,
	qrevmod,
	xbox,
	$get_boxes	= "deimos$src/x_deimos.e"

	# hide these tasks
    hidetask refl,trace,qtrmap,qrevmod

    type "deimos$motd"

clbye()

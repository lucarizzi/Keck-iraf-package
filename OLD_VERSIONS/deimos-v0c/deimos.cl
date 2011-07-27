#{ Package DEIMOS is Drew Phillips' deimos software

# load requried packages...
mscred	# needed for dmosdisplay
;

print ("\n DEIMOS (ver.0c for IRAF 2.12) -- Software in development -- User assumes risk\n")

cl < "deimos$lib/zzsetenv.def"
package	deimos, bin = deimosbin$

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

task	dsimulator,
	refl,
	trace,
	qtrmap,
	qmodel,
	qrevmod,
	xbox,
	$get_boxes	= "deimos$src/x_deimos.e"

    type "deimos$motd"

clbye()

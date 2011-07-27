#{ Package DEIMOS is Drew Phillips' deimos software

# load requried packages...
mscred	# needed for dmosdisplay
;

cl < "dproto$lib/zzsetenv.def"
package	dproto, bin = dprotobin$

    # script tasks...
	task	$cdata		= "dproto$src/cdata.cl"
	task	dss_chart	= "dproto$src/dss_chart.cl"
    task	$do_xbox	= "dproto$src/do_xbox.cl"
    task	$lastimage	= "$rsh polo lastimage"
    task	$wfi		= "$rsh polo wfi"
    task	get_seeing	= "dproto$src/get_seeing.cl"
    task	tune_qmodel	= "dproto$src/tune_qmodel.cl"
    task 	check_boxes	= "dproto$src/check_boxes.cl"
    task 	dmosdisplay	= "dproto$src/dmosdisplay.cl"
    task 	plotmask	= "dproto$src/plotmask.cl"

task	dsimulator,
	refl,
	trace,
	qtrmap,
	qmodel,
	qrevmod,
	xbox,
	$get_boxes	= "dproto$src/x_deimos.e"

    type "dproto$motd"

clbye()

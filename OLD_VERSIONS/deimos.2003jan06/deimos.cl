#{ Package DEIMOS is Drew Phillips' deimos software

# load requried packages...
mscred	# needed for dmosdisplay
;

package	deimos

    # script tasks...
    task 	dmosdisplay	= "deimos$src/dmosdisplay.cl"
    task 	check_boxes	= "deimos$src/check_boxes.cl"
    task	$wfi		= "$rsh polo wfi"
    task	$lastimage	= "$rsh polo lastimage"
    task	$do_xbox	= "deimos$src/do_xbox.cl"
    task	get_seeing	= "deimos$src/get_seeing.cl"
    task	tune_qmodel	= "deimos$src/tune_qmodel.cl"

    # compiled tasks...
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
		$xlate,
		$tvalign,
		$tmos,
		xbox,
		$get_boxes,
		cohu	= "deimos$src/x_deimos.e"

    type "deimos$motd"

clbye()

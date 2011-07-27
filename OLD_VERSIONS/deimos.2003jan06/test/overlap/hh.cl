trace htest.dat "" "" "" coll_ang=0. t2=0. o3=0. roll=0. x_opt=0 y_opt=0 mos_r=0. nord=1 gmm=900 refw=6736.
trace htest.dat "" "" "" coll_ang=0. t2=0. o3=0. roll=0. x_opt=0 y_opt=0 mos_r=0. nord=1 gmm=900 refw=7000.
trace htest.dat "" "" "" coll_ang=0. t2=0. o3=0. roll=0. x_opt=0 y_opt=0 mos_r=0. nord=1 gmm=900 refw=7258.
trace htest.dat "" "" "" coll_ang=0. t2=0. o3=0. roll=0. x_opt=0 y_opt=0 mos_r=0. nord=1 gmm=900 refw=7510.

del tmpa,tmpb,tmpc,tmpd
trace xym90 "" "" "" coll_ang=0. t2=0. o3=0. roll=0. x_opt=0 y_opt=0 mos_r=0. nord=1 gmm=900 refw=6736. > tmpa
trace xy0 "" "" "" coll_ang=0. t2=0. o3=0. roll=0. x_opt=0 y_opt=0 mos_r=0. nord=1 gmm=900 refw=7000. > tmpb
trace xyp90 "" "" "" coll_ang=0. t2=0. o3=0. roll=0. x_opt=0 y_opt=0 mos_r=0. nord=1 gmm=900 refw=7258. > tmpc
trace xyp180 "" "" "" coll_ang=0. t2=0. o3=0. roll=0. x_opt=0 y_opt=0 mos_r=0. nord=1 gmm=900 refw=7510. > tmpd

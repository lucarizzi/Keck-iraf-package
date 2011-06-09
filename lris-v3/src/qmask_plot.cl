# QMASK_PLOT: Quick mask-plot generating script -- calls gen_igi

procedure qmask_plot (input)

file	input {prompt = "set-up file"}
string	title {prompt = "name of mask"}
real	ht_page  {prompt = "height of plot area (in.)(portrait)"}
real	wid_page {prompt = "width of plot area (in.)(portrait)"}
real	slit_wid {prompt = "slit width in mm"}
bool	true_scale {prompt = "true-scale plot?"}

begin
{

# define local variables
file	inp, tmp_plot

# Construct names
inp = input
tmp_plot = mktemp ("tmpq")

if (! access (inp)) {
	beep
        print ""
        print ("Error: "//inp//" does not exist ! -- exiting")
        bye
}

# Generate the plot
#
gen_igi (inp, tmp_plot, tit=title, height=wid_page, wid=ht_page,\
slit_wid=slit_wid, true_scale=true_scale)
stsdas
igi (dev="stdplot", < tmp_plot)
gflush
delete (tmp_plot, veri-)
}
end

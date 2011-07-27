# COHU: swap lines in COHU images

include	<imhdr.h>

procedure	t_cohu()

char	input[SZ_FNAME]
char	output[SZ_FNAME]

int	nx, ny
int	i
pointer	im1, im2

pointer	imgl2r(), impl2r()
pointer	immap()

begin
	call clgstr ("input", input, SZ_FNAME)
	call clgstr ("output", output, SZ_FNAME)
	im1 = immap (input, READ_ONLY, 0)
	im2 = immap (output, NEW_COPY, im1)

	nx = IM_LEN (im1,1)
	ny = IM_LEN (im1,2)

	do i = 1, ny-1, 2 {
		call amovr (Memr[imgl2r (im1, i+1)], Memr[impl2r (im2, i)], nx)
		call amovr (Memr[imgl2r (im1, i)], Memr[impl2r (im2, i+1)], nx)
	}
	call imunmap (im2)
	call imunmap (im1)
end

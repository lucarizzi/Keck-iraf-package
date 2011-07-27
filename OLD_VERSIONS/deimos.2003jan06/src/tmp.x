
# CCDGEOM: specify the geometry of the CCDs in the mosaic
# Probably want to add FCS devices also (note that FDS CCDs are rotated 90deg)
# Order is x(pix), y(pix), theta(deg)

procedure ccdgeom (a)

real	a[NCCD,3]

begin
#	coeff    pix-off   nom.gap    adjustment

	a[1,1] = -4096.  - 100.
	a[1,2] = -4096.  -   3.3333
	a[1,3] = 0.

	a[2,1] = -2048.  -  33.3333
	a[2,2] = -4096.  -   3.3333
	a[2,3] = 0.

	a[3,1] = 0.      +  33.3333
	a[3,2] = -4096.  -   3.3333
	a[3,3] = 0.

	a[4,1] = 2048.   + 100.
	a[4,2] = -4096.  -   3.3333
	a[4,3] = 0.

	a[5,1] = -4096.  - 100.
	a[5,2] = 0.      +   3.3333
	a[5,3] = 0.

	a[6,1] = -2048.  -  33.3333
	a[6,2] = 0.      +   3.3333
	a[6,3] = 0.

	a[7,1] = 0.      +  33.3333
	a[7,2] = 0.      +   3.3333
	a[7,3] = 0.

	a[8,1] = 2048.   + 100.
	a[8,2] = 0.      +   3.3333
	a[8,3] = 0.

	a[9,1] = -4096.
	a[9,2] = -4096.
	a[9,3] = 0.
end

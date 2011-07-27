#
# ETC.X: Various small procedures that don't fit anywhere else:
#	bksp_strcat:	back-space, then string concat
#	hdr_mulr:	multiply header field by some factor (real)
#

# BKSP_STRCAT: A version of strcat that backspaces over the last character --
# usually a newline; especially useful for image history comments.
# eg: call bksp_strcat (", CR cleaned (crr)\n", IM_HISTORY(im), SZ_IMHIST)
#
	procedure	bksp_strcat (str, outstr, maxch)

char	str[ARB]			# string to append
char	outstr[maxch]			# output string
int	maxch				# max length of outstr

int	ndx

int	gstrcat(), strlen()

begin
	ndx = max (strlen (outstr), 1)
	outstr[ndx] = EOS
	ndx = gstrcat (str, outstr, maxch)
end

#
# HDR_MULR: multiply a header field (or alternate) by a real factor)
#
	procedure hdr_mulr (im, factor, field, afield, add)

pointer	im				# image pointer
real	factor				# multiplicative factor
char	field[ARB]			# field to be modified
char	afield[ARB]			# alternate field to be modified
bool	add				# add the alternate if not there?

bool	imaccf()
real	imgetr()
begin
	if (imaccf (im, field))
		call imputr (im, field, factor * imgetr (im, field))
	else if (imaccf (im, afield))
		call imputr (im, afield, factor * imgetr (im, afield))
	else
		call imaddr (im, afield, factor)
end

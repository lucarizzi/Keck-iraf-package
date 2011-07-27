# Miscellaneous useful routines go here.

# int procedure	line_count (fd)	-- returns count
# int procedure req_valr (str, curr, new, vmin, vmax) -- returns stat
# int procedure req_vali (str, curr, new, vmin, vmax) -- returns stat
# int procedure req_valb (str, curr, new) -- returns stat
# procedure get_valr (str, curr, new, vmin, vmax)
# procedure get_vali (str, curr, new, vmin, vmax)
# procedure get_valb (str, curr, new)


#
# LINE_COUNT: (I know I have another version of this around somewhere ...)
#

int procedure	line_count (fd)

pointer	fd			# file descriptor of open file

char	tchar
int	ndx

int	fscan(), nscan()
begin

# Count the entries
	ndx = 0
	while (fscan (fd) != EOF) {
		call gargwrd (tchar, 1)
		if (tchar == '#' || nscan() == 0) {
			next
		}
		ndx = ndx + 1
	}
	call seek (fd, BOF)

	return (ndx)
end


define	SZ_VALSTR 80
define	SZ_FMTSTR 80

#
# REQ_VALR: request a real value, allowing a default
#

int	procedure req_valr (str, curr, new, vmin, vmax)

char	str[ARB]				# value description
real	curr					# current value	## SPEC
real	new					# new value	## SPEC
real	vmin, vmax				# acceptible range	## SPEC

char	valstr[SZ_VALSTR]
real	val						## SPECIFIC

int	fscan(), nscan()

begin
	call printf ("\n%s (%g): ")			## SPEC
		call pargstr (str)
		call pargr (curr)			## SPEC
	call flush (STDOUT)
	if (fscan(STDIN) != EOF) {
		call gargwrd (valstr, SZ_VALSTR)
		if (nscan() < 1) {
			new = curr
			return (OK)
		} else {
			call sscan (valstr)
			call gargr (val)			## SPEC
			if (nscan() == 1) {
			    if (val <= vmax && val >= vmin) {
				new = val
				return (OK)
			    } else {
				call printf (
					"Value '%g' out of range [%g,%g]!\007  ")
					call pargr (val)	## SPEC
					call pargr (vmin)
					call pargr (vmax)
				call flush (STDOUT)
			    }
			} else {
			    call printf ("Invalid data-type!\007  ")
			    call flush (STDOUT)
			}
		}
	}
	return (ERR)
end


#
# GET_VALI: get an integer value, allowing a default
#


int	procedure req_vali (str, curr, new, vmin, vmax)

char	str[ARB]				# value description
int	curr					# current value	## SPEC
int	new					# new value	## SPEC
int	vmin, vmax				# acceptible range	## SPEC

char	valstr[SZ_VALSTR]
int	val						## SPECIFIC

int	fscan(), nscan()

begin
	call printf ("\n%s (%d): ")			## SPEC
		call pargstr (str)
		call pargi (curr)			## SPEC
	call flush (STDOUT)
	if (fscan(STDIN) != EOF) {
		call gargwrd (valstr, SZ_VALSTR)
		if (nscan() < 1) {
			new = curr
			return (OK)
		} else {
			call sscan (valstr)
			call gargi (val)			## SPEC
			if (nscan() == 1) {
			    if (val <= vmax && val >= vmin) {
				new = val
				return (OK)
			    } else {
				call printf (
					"Value '%d' out of range [%d,%d]!\007  ")
					call pargi (val)	## SPEC
					call pargi (vmin)
					call pargi (vmax)
				call flush (STDOUT)
			    }
			} else {
			    call printf ("Invalid data-type!\007  ")
			    call flush (STDOUT)
			}
		}
	}
	return (ERR)
end


#
# REQ_VALB: request a boolean value, allowing a default
#

int	procedure req_valb (str, curr, new)

char	str[ARB]				# value description
bool	curr					# current value	## SPEC
bool	new					# new value	## SPEC

char	valstr[SZ_VALSTR]

bool	itob(), streq()
int	fscan(), nscan()

begin
	call printf ("\n%s (%b): ")
		call pargstr (str)
		call pargb (curr)
	call flush (STDOUT)
	if (fscan(STDIN) != EOF) {
		call gargwrd (valstr, SZ_VALSTR)
		if (nscan() < 1) {
			new = curr
			return (OK)
		} else if (streq (valstr, "yes") || streq (valstr, "y")) {
			new = itob (YES)
			return (OK)
		} else if (streq (valstr, "no") || streq (valstr, "n")) {
			new = itob (NO)
			return (OK)
		} else {
			call printf ("Invalid; must be (y|n)\007  ")
			call flush (STDOUT)
		}
	}

	return (ERR)
end


#
# GET_VALR: get a real value
#

procedure	get_valr (str, curr, new, vmin, vmax)

char	str[ARB]				# value description
real	curr					# current value	## SPEC
real	new					# new value	## SPEC
real	vmin, vmax				# acceptible range	## SPEC

int	req_valr()
begin
	while (req_valr (str, curr, new, vmin, vmax) != OK)
		;
end


#
# GET_VALI: get an integer value
#

procedure	get_vali (str, curr, new, vmin, vmax)

char	str[ARB]				# value description
int	curr					# current value	## SPEC
int	new					# new value	## SPEC
int	vmin, vmax				# acceptible range	## SPEC

int	req_vali()
begin
	while (req_vali (str, curr, new, vmin, vmax) != OK)
		;
end


#
# GET_VALB: get a bool value
#

procedure	get_valb (str, curr, new)

char	str[ARB]
bool	curr
bool	new

int	req_valb()
begin
	while (req_valb (str, curr, new) != OK)
		;
end

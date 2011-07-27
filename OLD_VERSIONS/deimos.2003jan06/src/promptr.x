# promptr.x / GDW / 24 Oct 2000
# 
# Purpose:
#	Prompt the user for a real value, with default value and 
#	optional range checking.
#-----------------------------------------------------------------------

real procedure promptr( prompt, maxch, default, minval, maxval)

char	prompt[ARB]		# prompt to print
int	maxch			# max length of prompt
real	default			# default value
real	minval			# minimum allowed value
real	maxval			# maximum allowed value
#---
real	answer			# user input
int	stat			# status of scan

int	scan()

begin

    # start loop...
    while( true ){

        # clean input buffer...
        call reset_scan()

        # print prompt...
        call eprintf( "\n%s [%g]: ")
            call pargstr( prompt, maxch)
            call pargr( default)
        call flush( STDERR)

        # get input...
        stat = scan() 

        # quit if no input...
        if( stat == EOF){ return default }

        # parse input and convert to lowercase...
        call gargr( answer)

        # check for valid response..
        if ( (minval != INDEF && answer < minval) 
            || (maxval != INDEF && answer > maxval)){
            if ( minval == INDEF ) {
                call eprintf("ERROR: value must be <= %g!  ")
                call pargr( maxval)
            } else if ( maxval == INDEF ) {
                call eprintf("ERROR: value must be >= %g!  ")
                call pargr( minval)
            } else {
                call eprintf("ERROR: value must be in range [%g,%g]!  ")
                call pargr( minval)
                call pargr( maxval)
            }
        } else {
            return answer
        }

    }

end

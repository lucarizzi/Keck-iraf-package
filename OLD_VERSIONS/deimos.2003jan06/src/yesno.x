# yesno.x / GDW / 24 Oct 2000
# 
# Purpose:
#	Prompt the user for a yes or no answer.  Return a boolean.
#-----------------------------------------------------------------------

bool procedure yesno( prompt, maxch, default)

char	prompt[ARB]		# prompt to print
int	maxch			# max length of prompt
bool	default			# default value
#---
char	answer			# user input
int	stat			# status of scan
char	def			# default answer (char)

bool	streq()
int	scan()
char	chrlwr()

begin

	# define default answer...
	if ( default ){
		call strcpy( "y", def, 1)
	} else {
		call strcpy( "n", def, 1)
	}

	# start loop...
	while( true ){

		# clean input buffer...
		call reset_scan()

		# print prompt...
		call eprintf( "\n%s (y/n) [%s]: ")
			call pargstr( prompt, maxch)
			call pargstr( def, 1)
		call flush( STDERR)

		# get input...
		stat = scan() 

		# quit if no input...
		if( stat == EOF){ return default }

		# parse input and convert to lowercase...
		call gargstr( answer, 1)
		answer = chrlwr( answer)

		# check for valid response..
		if ( streq( answer, "y")){
			return true
		} else if ( streq( answer, "n")){
			return false
		} else if ( streq( answer, "")){
			return default
		}

	}

end

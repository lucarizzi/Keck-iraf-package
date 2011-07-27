#-----------------------------------------------------------------------
procedure cdata()
#-----------------------------------------------------------------------
# Purpose:
#	Change current directory to the DEIMOS data directory on polo
#-----------------------------------------------------------------------

begin
	string	host="deimosserver"
	string	dir

	rsh( host, " show -s deiccd -terse outdir") | scan( dir)
	dir = "/s" // dir
	chdir( dir)
end

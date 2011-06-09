#-----------------------------------------------------------------------
procedure cdata()
#-----------------------------------------------------------------------
# Purpose:
#	Change current directory to the lris data directory on polo
#-----------------------------------------------------------------------

begin
	string	host="lrisserver"
	string	dir

	rsh( host, " show -s lris -terse outdir") | scan( dir)
	dir = "/s" // dir
	chdir( dir)
	path
end

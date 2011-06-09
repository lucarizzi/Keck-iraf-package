# Name:
#	Makefile
# 
# Purpose:
# 	Generate files for NIRSPEC package
# 
# Type:
#	GNU Makefile format
# 
# Usage:
#	gmake
# 
# Author:
#	Gregory D. Wirth
# 
# Note:
#	This file is only required for development of the package, 
#	not for installation at other sites.  Please refer to the
#	file README.install for instructions on installing the package
#	for use in IRAF.
# 
# Revision history:
#	2001-Apr-17	GDW	Original version
#-----------------------------------------------------------------------

# define variables...
doc_dir = doc/
html_files = $(wildcard *.html)
lroff_files = $(subst html,hlp,$(html_files))
tar_file = ../NIRSPEC.tar.gz

# generate IRAF-compatible documentation from HTML files and 
# store in documenation directory...
documentation: $(html_files)
	html2lroff $^
	-mkdir $(doc_dir)
	mv $(subst html,hlp,$^) $(doc_dir)

# generate file listing OH lines...
lowd_ir_ohlines.dat:	ir_ohlines.dat lowd_ir_ohlines.awk
	awk -f lowd_ir_ohlines.awk ir_ohlines.dat > lowd_ir_ohlines.dat

distrib: clean
	cd .. ; mktarfile NIRSPEC nirspec

# remove clutter...
clean:
	-rm \#*\# .*~ *~

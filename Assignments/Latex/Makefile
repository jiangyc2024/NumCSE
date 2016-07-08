all: setup NCSEProblems.pdf clean

svn_url="https://svn.id.ethz.ch/sam/Numcourses/ncse/ncse_new"

#
# Download packages, numcse_macro and aux file from svn
#
setup: NCSENEW_macros.tex samlecture_macros.tex samlecturenotes.sty NCSENEW.aux

NCSENEW_macros.tex:
	svn export $(svn_url)/$@
	chmod 555 $@

samlecture_macros.tex:
	svn export $(svn_url)/$@
	chmod 555 $@

samlecturenotes.sty:
	svn export $(svn_url)/$@
	chmod 555 $@

NCSENEW.aux:
	svn export $(svn_url)/$@
	chmod 555 $@

# todo: download images

#
# Compile latex
#
NCSEProblems.dvi: NCSEProblems.tex
	latex -interaction=scrollmode -shell-escape $<
	latex -interaction=scrollmode -shell-escape $<

%.pdf: %.dvi
	dvipdf $<

#
# Cleanup
#
clean:
	rm -f *.bbl *.aux *.log *.blg *~ *.bak *.ps *.dvi *.out *.tmp *.mtc *.mtc0 *.mw *.maf *.edx *.idx *.ddx *.sdx

cleanexternal:
	rm -f NCSENEW_macros.tex samlecture_macros.tex samlecturenotes.sty

cleanall: clean cleanexternal
	rm -f NCSEProblems.pdf
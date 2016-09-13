all: setup NCSEProblems.pdf clean

svn_url="https://svn.id.ethz.ch/sam/Numcourses/ncse/ncse_new"
slidemode="4"

#
# Download packages, numcse_macro and aux file from svn
#
setup: NCSENEW_macros.tex samlecture_macros.tex samlecturenotes.sty NCSENEW.aux NCSEProblems.tex rev.tex

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

rev.tex:
	svn export $(svn_url)/$@
	chmod 555 $@
NCSEProblems.tex:
	svn export $(svn_url)/$@
	chmod 555 $@


# todo: download images

#
# Compile latex
#
define LATEXCOMPILE
latex -interaction=scrollmode -shell-escape "\newcommand{\slidemode}{$(slidemode)} \input{$1}"
endef

NCSEProblems.dvi: NCSEProblems.tex
	$(call LATEXCOMPILE,$<)
	$(call LATEXCOMPILE,$<)

%.pdf: %.dvi
	dvipdf $<

#
# Cleanup
#
clean:
	find *.aux ! -name 'NCSENEW.aux' -type f -exec rm -f {} +
	rm -f *.bbl *.log *.blg *~ *.bak *.ps *.dvi *.out *.tmp *.mtc *.mtc0 *.mw *.maf *.edx *.idx *.ddx *.sdx

cleanexternal:
	rm -f NCSENEW_macros.tex samlecture_macros.tex samlecturenotes.sty NCSENEW.aux rev.tex

cleanall: clean cleanexternal
	rm -f NCSEProblems.pdf

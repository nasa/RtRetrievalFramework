all: l2.pdf

l2.pdf: l2.dvi
	dvipdf $^

l2.dvi: l2.tex
	latex $^
	latex $^

clean:
	-/bin/rm l2.aux l2.dvi l2.log l2.toc *~

distclean: clean
	-/bin/rm l2.pdf


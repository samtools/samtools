.SUFFIXES: .gp .tex .eps .pdf .eps.gz

.eps.pdf:
		epstopdf --outfile $@ $<

.eps.gz.pdf:
		gzip -dc $< | epstopdf --filter > $@

all:samtools.pdf

samtools.pdf:samtools.tex ahmm.pdf
		pdflatex samtools; pdflatex samtools

clean:
		rm -fr *.toc *.aux *.bbl *.blg *.idx *.log *.out *~ notes.pdf samtools.pdf

# latex Makefile
texpath=/usr/texbin/
PDFLATEX=${texpath}pdflatex -halt-on-error -synctex=1 --interaction=nonstopmode
LATEX=${PDFLATEX}
BIBTEX=bibtex
DVIPS=dvips
PS2PDF=ps2pdf

all: brick_letter

brick_letter: 
	@rm -f brick_letter*.aux brick_letter*.bbl brick_letter*.blg brick_letter*.dvi brick_letter*.log brick_letter*.lot brick_letter*.lof
	${LATEX} brick_letter.tex
	${BIBTEX} brick_letter
	${LATEX} brick_letter.tex
	${BIBTEX} brick_letter
	${LATEX} brick_letter.tex

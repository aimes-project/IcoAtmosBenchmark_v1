###########################################################
#
# Makefile for compiling tex document
#
###########################################################

# target file name
TEXNAME  = IAB_DYNAMICO_kernels

SRCS = $(wildcard *.tex)

all: pdf

pdf: $(TEXNAME).pdf
	install $< ../..

clean:
	rm -f *.pdf *.out *.aux *.bbl *.blg *.dvi *.fdb_latexmk *.fls *.log *.synctex.gz *.toc

$(TEXNAME).pdf:
	latexmk $(TEXNAME)

.SUFFIXES:
.SUFFIXES: .tex .bbl .dvi .pdf .bib .html

.PHONY : pdf clean

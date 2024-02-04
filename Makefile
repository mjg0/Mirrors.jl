.PHONY: all thesis clean distclean

all: thesis

thesis: thesis.pdf

thesis.pdf: thesis/thesis.tex
	@mkdir -p thesis/out
	@cd thesis && pdflatex -output-directory=out thesis.tex
	@cp thesis/out/thesis.pdf .

clean:
	@rm -f thesis.pdf

distclean: clean
	@rm -rf thesis/out
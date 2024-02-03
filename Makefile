.PHONY: thesis all clean

all: thesis

thesis: thesis.pdf

thesis.pdf:
	@mkdir -p thesis/out
	@cd thesis && pdflatex -output-directory=out thesis.tex
	@cp thesis/out/thesis.pdf .

clean:
	@rm -f thesis.pdf
	@rm -rf thesis/out

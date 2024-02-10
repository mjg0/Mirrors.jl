.PHONY: all thesis clean distclean

DIR_STRUCTURE = thesis/code-directory-structure.txt

all: thesis

thesis: thesis.pdf

$(DIR_STRUCTURE):
	@(cd .. && find Mirrors.jl -type f -not -path '*/examples/*' -iname '*.jl' -or -iname '*.toml' | sort) > $@

thesis.pdf: thesis/thesis.tex $(DIR_STRUCTURE)
	@mkdir -p thesis/out
	@cd thesis && pdflatex -shell-escape -output-directory=out thesis.tex
	@cp thesis/out/thesis.pdf .

clean:
	@rm -f thesis.pdf

distclean: clean
	@rm -rf thesis/out $(DIR_STRUCTURE)
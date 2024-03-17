.PHONY: all thesis clean distclean

DIR_STRUCTURE = thesis/code-directory-structure.txt
CODE_TEX = thesis/code.tex

all: thesis

thesis: thesis.pdf

$(DIR_STRUCTURE):
	@(cd .. && find Mirrors.jl -type f -not -path '*/examples/*' -iname '*.jl' -or -iname Project.toml | sort) > $@

$(CODE_TEX): $(DIR_STRUCTURE) src/*.jl ext/*.jl test/*.jl Project.toml
	@bash thesis/generate-code-tex.sh $< $@

thesis.pdf: thesis/thesis.tex $(DIR_STRUCTURE) $(CODE_TEX) thesis/BYUPhys.cls thesis/references.bib
	@mkdir -p thesis/out
	@cd thesis && lualatex -shell-escape -output-directory=out thesis.tex \
	           && bibtex out/thesis.aux \
	           && lualatex -shell-escape -output-directory=out thesis.tex \
	           && lualatex -shell-escape -output-directory=out thesis.tex
	@cp thesis/out/thesis.pdf .

clean:
	@rm -f thesis.pdf

distclean: clean
	@rm -rf thesis/out $(DIR_STRUCTURE) $(CODE_TEX)
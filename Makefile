RMD_FILES := $(patsubst %.Rmd, %.md ,$(wildcard paper_sections/*.Rmd))
TEX_FILES := $(patsubst %.md, %.tex ,$(wildcard paper_sections/*.md))

tex: md
main: tex
all: main

md: $(RMD_FILES)
%.md: %.Rmd
				R --vanilla --slave -e "knitr::knit('$<', output='$@')"

tex: $(TEX_FILES)
%.tex: %.md
				pandoc --natbib $< -o $@



main: paper.Rnw
				R --vanilla --slave -e "knitr::knit('paper.Rnw')"
				pdflatex paper.tex
				bibtex paper
				pdflatex paper.tex
				pdflatex paper.tex
				rm paper.tex


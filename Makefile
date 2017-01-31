RMD_FILES := $(patsubst %.Rmd, %.md ,$(wildcard paper_sections/*.Rmd))
TEX_FILES := $(patsubst %.md, %.tex ,$(wildcard paper_sections/*.md))

all: main tex md
main: tex md
tex: md

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
				rm paper_sections/*.md paper_sections/*.tex


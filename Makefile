all: main example_data md fig tex compare
main: example_data md fig tex 
compare: main
makezip: main

run_comparison = false

example_data: data/bird_move.csv data/zooplankton_example.csv
				R --vanilla --slave -e "source('code/bird_example_data.R')"
				R --vanilla --slave -e "source('code/cleaning_zooplankton_data.R')"

md: paper_sections/01-intro.Rmd paper_sections/02-gams.Rmd paper_sections/03-hierarchical_gams.Rmd paper_sections/04-examples.Rmd paper_sections/05-computational_and_statistical_issues.Rmd 

fig: figures/temp_growth_example.pdf code/alternate_models.R 
				R --vanilla --slave -e "source('code/alternate_models.R')"

tex: paper_sections/bibliography.bib paper_sections/peerj.csl paper_sections/preamble.sty

main: paper_sections/full_document.Rmd
				R --vanilla --slave -e "library(knitr); purl('paper_sections/full_document.Rmd',documentation =0, output = 'compiled_paper/supplemental_code.R')"
				R --vanilla --slave -e "library(rmarkdown); render('paper_sections/full_document.Rmd',output_file = 'full_document.tex',knit_root_dir='../')"
				sed -i.bak '/xcolor/d' paper_sections/full_document.tex
				R --vanilla --slave -e "setwd('paper_sections'); library(tools); texi2pdf('full_document.tex',clean = TRUE)"
				mv paper_sections/full_document.pdf compiled_paper/full_document.pdf 
				mv paper_sections/full_document.tex compiled_paper/full_document.tex 
				
compare: compiled_paper/full_document.tex compiled_paper/prior_submission.tex
        ifeq ($(run_comparison),true)
					latexdiff --config="PICTUREENV=(?:picture|DIFnomarkup|table)[\w\d*@]*" compiled_paper/prior_submission.tex compiled_paper/full_document.tex  > compiled_paper/diff.tex
					-pdflatex -interaction=nonstopmode compiled_paper/diff.tex -output-directory compiled_paper
					-pdflatex -interaction=nonstopmode compiled_paper/diff.tex -output-directory compiled_paper
					rm compiled_paper/diff.log compiled_paper/diff.aux compiled_paper/diff.out
        endif

makezip: compiled_paper/supplemental_code.R data/bird_move.csv data/zooplankton_example.csv data/bird_move_global.csv mixed-effect-gams.Rproj
				zip compiled_paper/supplemental_info.zip data/bird_move.csv data/zooplankton_example.csv data/bird_move_global.csv mixed-effect-gams.Rproj
				cp compiled_paper/supplemental_code.R code/supplemental_code.R 
				zip -ur compiled_paper/supplemental_info.zip code/supplemental_code.R
				rm -f compiled_paper/supplemental_code.R

				
				
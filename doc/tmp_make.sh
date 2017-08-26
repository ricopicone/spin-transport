#!/bin/bash
set -x  # display all commands in output

# Safe execution of a Unix command: exit if failure
function system {
  "$@"
  if [ $? -ne 0 ]; then
    echo "make.sh: unsuccessful command $@"
    echo "abort!"
    exit 1
  fi
}
\nsystem doconce format html doc  --html_style=bootswatch_readable --html_code_style=inherit --html_pre_style=inherit --toc_depth=2 --pygments_html_style=default --html_style=bootswatch_journal --html_output=doc \n\nsystem doconce format ipynb doc  \n\nsystem doconce format pdflatex doc  \n\nsystem rm -f doc.aux\n\nsystem doconce format pdflatex doc --latex_style=std --latex_title_layout=std --latex_font=palatino --latex_admon=yellowicon --latex_admon_color=yellow!5 --latex_fancy_header --latex_code_style=pyg-gray --minted_latex_style=friendly --latex_section_headings=black --latex_colored_table_rows=blue --latex_preamble=latex_preamble.tex \nsystem doconce replace \title{ \title\{\\sffamily\\bfseries\{\}  doc.tex\nsystem doconce replace \author{ \author{\\sffamily{}  doc.tex\nsystem doconce replace \date{ \date{\\sffamily{}  doc.tex\nsystem pdflatex -file-line-error -interaction nonstopmode doc\nsystem pdflatex -file-line-error -interaction nonstopmode doc\nsystem pdflatex -file-line-error -interaction nonstopmode doc\n
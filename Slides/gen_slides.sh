#!/bin/sh

pdflatex slides_graphene.tex
pdflatex slides_graphene.tex
rm *.aux *.log *.out *.toc *.nav *.snm

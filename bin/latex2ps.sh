#!/bin/bash

input=$1 
base=${input%.tex}

latex $input
dvips -o $base.ps $base.dvi
rm $base.out $base.snm $base.aux $base.nav $base.toc $base.dvi

#!/bin/bash

for i in $1 ; do { f=${i%.eps} ; convert -density 300 -rotate 90 $i $f.jpg ; } ; done
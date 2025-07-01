#!/bin/bash
#this script batch converts pdf files into pngs.

#to run this you must install poppler: https://poppler.freedesktop.org/

#inputs
PDF_DIR="/Users/gnegreira/Dropbox/work/ITM/Scripts/Mine/Scripts for My papers/2025 - Atrandi paper/atrandi_paper"
RESOLUTION=600

#script
cd "$PDF_DIR"
for file in *.pdf; 
    do 
    pdftoppm -png -r $RESOLUTION "$file" "${file%.pdf}"; 
done
 
#set png resolution
RESOLUTION=600

#run any R script containing 'Figure' in its name
find . -maxdepth 1 -type f -name '*Figure*.R' -exec Rscript {} \;

#convert any PDF file to png.
for file in ./*.pdf; do
  if [ -e "$file" ]; then
    pdftoppm -png -r "$RESOLUTION" "$file" "${file%.pdf}"
  fi
done

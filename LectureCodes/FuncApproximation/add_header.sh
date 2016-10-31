#! /bin/bash

shopt -s globstar extglob

echo Header:
cat header.txt
echo "Press enter to continue..."
read

printf "\n" >> header.txt

for f in **/*.@(cpp|gpp); do
  echo Processing $f
  cat header.txt $f > $f.new
  mv $f.new $f
done

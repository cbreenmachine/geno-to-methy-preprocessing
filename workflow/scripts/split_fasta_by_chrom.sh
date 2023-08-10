#!/bin/bash
# c/o https://crashcourse.housegordon.org/split-fasta-files.html

# Parse inputs
ifile="${1}"


odir="$(dirname "$ifile")/by-chrom"
mkdir ${odir}

csplit -s ${ifile} '/>/' '{100}'

for i in xx*
do
  n=$(sed 's/>// ; s/ .*// ; 1q' "$i")
  mv "$i" "$n.fa"
done

mv *fa ${odir}

#END
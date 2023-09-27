for i in *.fa; do
  j="/home/marinas/rbpmap/UCSC/hg38/${i%.*}.nib";
  "../faToNib" "/home/marinas/rbpmap/UCSC/hg38/${i}" "${j}";
  done

#local change

for i in *.fa; do
  j="/home/marina/rbpmap/UCSC/hg38/${i%.*}.nib";
  "../faToNib" "/home/marina/rbpmap/UCSC/hg38/${i}" "${j}";
  done

#local change

for i in *.fa; do
  j="/home/hanah/RBPmap_1.2/UCSC/hg38/${i%.*}.nib";
  "../faToNib" "/home/hanah/RBPmap_1.2/UCSC/hg38/${i}" "${j}";
  done

#local change

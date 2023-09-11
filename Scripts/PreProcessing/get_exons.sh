#!/bin/bash


json_file="paths.json"

field1="\"Reference genome\""
field2="\"Reference fasta\""

# Get the paths to reference files
gff3=$(jq -r ".[\"$1\"][$field1]" "$json_file")
fasta=$(jq -r ".[\"$1\"][$field2]" "$json_file")

#prepare genome file
samtools faidx $fasta

# get all transcripts of protein coding genes with TSL 1-3
grep -E "transcript_support_level=[123]" $gff3 | 
awk -F'\t' -v OFS='\t' '$3 == "transcript" && $9 ~ /protein_coding/ {
    match($9, /gene_id=([^;]+)/, gene_id);
    match($9, /transcript_id=([^;]+)/, transcript_id);
    print $1, $4, $5, gene_id[1], transcript_id[1], $7;
}' > ${1}_0_Files/gencode_TSS.tsv


# get TSS
awk 'BEGIN {OFS="\t"} {
    if ($6 == "+") {ed
        $7 = $2
    } else if ($6 == "-") {
        $7 = $3
    }
    print
}' ${1}_0_Files/gencode_TSS.tsv > ${1}_0_Files/TSS.tsv

#get TSS coordinates as bed file
awk 'BEGIN {OFS="\t"} {print $1,$7, $7+1, "TSS", ".", $6}' ${1}_0_Files/TSS.tsv | sort | uniq > ${1}_0_Files/TSS.bed

# Get all exons of transcripts with TSL 1-3
grep -E "transcript_support_level=[123]" $gff3 | grep "protein_coding" | awk -F'\t' -v OFS='\t' '$3 == "exon" { print $1, $4, $5, "Exon", ".", $7}' | sort | uniq > ${1}_0_Files/all_exons.bed

#extend exon body by 200bp
bedtools slop -i ${1}_0_Files/all_exons.bed -g $fasta.fai  -b 200 > ${1}_0_Files/exons_flanked.bed

#Drop flanked exons overlapping with TSS
bedtools intersect -wa -a ${1}_0_Files/exons_flanked.bed -b ${1}_0_Files/TSS.bed -s -v > ${1}_0_Files/filtered_flanked_exons.bed

# get exon coordinates (remove flanking regions)
bedtools slop -i ${1}_0_Files/filtered_flanked_exons.bed -g $fasta.fai -l -200 -r -200 -s > ${1}_0_Files/exon_coords.bed

#cleanup
rm ${1}_0_Files/*flanked*.bed
rm ${1}_0_Files/*.tsv

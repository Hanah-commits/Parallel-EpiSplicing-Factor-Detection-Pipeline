library(GenomicFeatures)
library(jsonlite)

# get the gff3 file path
json_data <- fromJSON("paths.json")
json_data <- lapply(json_data, function(x) {
  x[sapply(x, is.null)] <- NA
  unlist(x)
})

do.call("rbind", json_data)
gff3.file <- json_data[["Reference genome"]]


## make a transcript database from the input gff3 file
gff3 <- GenomicFeatures::makeTxDbFromGFF(organism = "Homo sapiens", 
                                    format = "gff3",
                                    file =  gff3.file)

## extract exons from the transcript database
exons <- exons(gtf)
exons
## convert exonsanges object into a dataframe
df.exon <- data.frame(seqnames=seqnames(exons),
                 starts=start(exons)-1,
                 ends=end(exons),
                 names=c(rep(".", length(exons))),
                 scores=c(rep(".", length(exons))),
                 strands=strand(exons))

# write exon coordinates into bed file
write.table(df.exon[, c(1,2,3)], file="0_Files/exon_coords.bed", quote=F, sep="\t", row.names=F, col.names=F)

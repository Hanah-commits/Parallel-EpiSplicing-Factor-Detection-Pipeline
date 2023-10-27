library(edgeR)
library(limma)
library(jsonlite)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("Process name must be supplied", call.=FALSE)
} 

proc = args[1]

json_data <- fromJSON("paths.json") 
proc_json_data <- json_data[[proc]]

proc_json_data <- lapply(proc_json_data, function(x) {
  x[sapply(x, is.null)] <- NA
  unlist(x)
})

do.call("rbind", proc_json_data)


tissue1 <- proc_json_data$tissue1
tissue2 <- proc_json_data$tissue2
data_dir <- proc_json_data$"RNASeq files"
tissue1_count <- length(list.files(data_dir, pattern=paste(tissue1,".+bam$", sep=""), full.names = TRUE, recursive = FALSE))
tissue2_count <- length(list.files(data_dir, pattern=paste(tissue2,".+bam$", sep=""), full.names = TRUE, recursive = FALSE))

file = paste(proc, "_0_Files/" , tissue1 , '_' , tissue2 , '_counts', sep="")

counts <- read.delim(file, skip=1, header=TRUE, sep='\t', row.names = 1)
counts = subset(counts, select=-c( Chr, Start, End, Strand, Length))


# Get normalization factors
dge <- DGEList((counts))
dge <- calcNormFactors(dge)

# specify model to be fitted
sample <- factor(rep(c(tissue1, tissue2), times = c(tissue1_count, tissue2_count)))
design.mat <- model.matrix(~0+sample)
colnames(design.mat) <- levels(sample)

# log fold-changes are obtained as contrasts
contrast.mat <- makeContrasts(
  Diff = get(tissue1) - get(tissue2), 
  levels = design.mat
)

# Counts are transformed to logCPM
y <- voom(dge, design.mat)


fit = lmFit(y, design.mat)
fit = contrasts.fit(fit, contrast.mat)
fit = eBayes(fit)

write.table(coef(fit),paste(proc, "_0_Files/logFC.csv"),sep="\t", col.names = FALSE)

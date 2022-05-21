library(edgeR)
library(limma)
library(jsonlite)

json_data <- fromJSON("paths.json") 
json_data <- lapply(json_data, function(x) {
  x[sapply(x, is.null)] <- NA
  unlist(x)
})

do.call("rbind", json_data)


tissue1 <- json_data$tissue1
tissue2 <- json_data$tissue2
tissue1_count <- as.numeric(json_data$tissue1_count)
tissue2_count <- as.numeric(json_data$tissue2_count)


file = paste("0_Files/" , tissue1 , '_' , tissue2 , '_counts', sep="")

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

write.table(coef(fit),"0_Files/logFC.csv",sep="\t", col.names = FALSE)

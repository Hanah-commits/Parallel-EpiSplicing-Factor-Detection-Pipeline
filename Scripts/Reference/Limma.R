library(edgeR)
library(limma)
library(rstudioapi)

# Getting the path of current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))

counts <- read.delim("ectodermalcell_H1_counts", skip=1, header=TRUE, sep='\t', row.names = 1)
counts = subset(counts, select=-c( Chr, Start, End, Strand, Length))
print(head(counts))
#write.table(counts, file='count_matrix.csv',  sep="\t")

# Get normalization factors
dge <- DGEList((counts))
dge <- calcNormFactors(dge)

# specify model to be fitted
sample <- factor(rep(c('ectodermalcell', 'H1'), times = c(6, 2)))
design.mat <- model.matrix(~0+sample)
colnames(design.mat) <- levels(sample)

# log fold-changes are obtained as contrasts
contrast.mat <- makeContrasts(
  Diff = ectodermalcell - H1,
  levels = design.mat
)

# Counts are transformed to logCPM
y <- voom(dge, design.mat)

fit = lmFit(y, design.mat)
fit = contrasts.fit(fit, contrast.mat)
fit = eBayes(fit)

print(head(coef(fit)))
write.table(coef(fit),"logFC.csv",sep="\t", col.names = FALSE)

top.table <- topTable(fit, n = Inf)
write.table(top.table,"limma.csv",sep="\t", col.names = TRUE)




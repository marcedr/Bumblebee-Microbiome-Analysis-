### ======================================================================= ###
### Title: Bumblebee Microbiome Analysis - Differential Abundance Analysis

### Description: This script outlies the DESeq2 used in 
### "Field-relevant stressors alter the bumblebee gut microbial community".
### Date of the last modification: Oct 6th, 2026
### ======================================================================= ###

### Housekeeping

rm(list=ls())

### Installing packages and loading libraries

# biocLite("DESeq2")
# source("https://bioconductor.org/biocLite.R")
library(phyloseq)
library(DESeq2)



### ======================================================================= ###
### Read in table
AllMergedFull <- read.csv(file.choose())
str(AllMergedFull)

### Give sample names
sample_names <- as.character(AllMergedFull[1, 2:ncol(AllMergedFull)])
length(sample_names)
sample_names

### Give treatment names
treatment <- as.character(AllMergedFull[2, 2:ncol(AllMergedFull)])
length(treatment)
treatment

### Prepare factors for DESeq2
Factors4Deseq <- data.frame(Treatment = treatment, row.names = sample_names)
Factors4Deseq$Treatment <- factor(Factors4Deseq$Treatment)
dim (Factors4Deseq)

### Prepare table for DESeq2 analysis
Counts4Deseq <- as.matrix(AllMergedFull[4:nrow(AllMergedFull), 2:ncol(AllMergedFull)])
dim(Counts4Deseq)

rownames(Counts4Deseq) <- AllMergeFull[4:nrow(AllMergedFull), 1]
colnames(Counts4Deseq) <- as.character(AllMergedFull[1, 2:ncol(AllMergedFull)])
head(Counts4Deseq)

### Check for negative values
Counts4Deseq <- matrix(
  as.numeric(trimws(Counts4Deseq)),
  nrow = nrow(Counts4Deseq),
  ncol = ncol(Counts4Deseq),
  dimnames = dimnames(Counts4Deseq)
)
class(Counts4Deseq)
typeof(Counts4Deseq)
min(Counts4Deseq)

### Sanity check
all(colnames(Counts4Deseq) == rownames(Factors4Deseq))
table(Factors4Deseq$Treatment)

### ======================================================================= ###
### Run DESeq2
LasiDESeq2Obj <- DESeqDataSetFromMatrix(
  countData = Counts4Deseq,
  colData = Factors4Deseq,
  design = ~Treatment
)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(LasiDESeq2Obj), 1, gm_mean)
diagdds = estimateSizeFactors(LasiDESeq2Obj, geoMeans = geoMeans)
DESeq2diagdds = DESeq(diagdds, fitType="local")


treatments=c("Ctrl", "Par", "PFP", "NFP")
combinations=t(combn(treatments, 2))

Results=as.list(rep(NA, nrow(combinations)))
for (i in seq_len(length(Results))) {
	names(Results)[i]=paste(combinations[i,1],combinations[i,2], sep="vs")
}

for (i in seq_len(length(Results))) {
	Results[[i]]=lfcShrink(DESeq2diagdds,
						contrast=c("AllMerged.Treatment", combinations[i,1], combinations[i,2]),
						type="normal"
					)
}


Results2=lapply(Results, function(x) cbind(as(x,"data.frame")))

for (i in seq_len(length(Results2))) {
write.csv(Results2[[i]], file=paste(names(Results2)[i], "_ALL.txt", sep="") )
write.csv(Results2[[i]][which(Results2[[i]]$padj<0.05),], file=paste(names(Results2)[i], "_Significant.txt", sep="") )
}


sink("FDR_to_001.txt")
lapply(Results2, function(x) x[which(x$padj<0.01),])
sink()

#################################
### Fix NAs in DESeq2 results ###
#################################

lapply(Results2, function(X) anyNA(X$pvalue))
lapply(Results2, function(X) anyNA(X$padj))
anyNA(Results2$CtrlvsPar$padj)
# Only one case containing NA

cbind(p.adjust(Results2$CtrlvsPar$pvalue),Results2$CtrlvsPar$padj)

cbind(p.adjust(Results2$CtrlvsPFP$pvalue, method="BH"),
      Results2$CtrlvsPFP$padj)

p.adjust(Results2$CtrlvsPFP$pvalue, method="BH")==
  Results2$CtrlvsPFP$padj
# Using BH you get correct p values (double-check just in case)

Results2$CtrlvsPar$padj=p.adjust(Results2$CtrlvsPar$pvalue, method="BH")
# applying correction

write.csv(Results2$CtrlvsPar,
          file=paste("CtrlvsPar_noNA", "_ALL.txt", sep="") )
write.csv(Results2$CtrlvsPar[which(Results2$CtrlvsPar$padj<0.05),],
          file=paste("CtrlvsPar_noNA", "_Significant.txt", sep="") )
# Export significant

save.image(compress=TRUE)

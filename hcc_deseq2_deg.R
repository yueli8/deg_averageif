library(DESeq2)
library(data.table)
library(limma)
library(calibrate)

#XENA-TCGAGTEx-LIHC
setwd("~/nanjing/tcgagtex")
DF<-read.table("lihc01.txt", head=TRUE)
b<-aggregate(DF[, -c(1:2)], by=list(DF$EntrezID, DF$Name), mean)
cts<-b[,-1]
write.table(cts,"cts")
cts<-read.table("cts",head=TRUE,row.names = "Group.2")
cts<-cts[,-1]
cts_round<-round(cts)
countData<-cts_round
condition <- factor(c(rep("p",110),rep("m",421)), levels = c("p","m"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
countData[countData==-1]<-0 #Error in DESeqDataSet(se, design = design, ignoreRank) : some values in assay are negative
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("m","p"))
dds$condition <- relevel(dds$condition, ref = "p")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"tcgagtex_norm.csv")
dds2<-DESeq(dds)
resultsNames(dds2)
res <- results(dds2)
res <- results(dds2, contrast=c("condition","m","p"))
resultsNames(dds2)
write.csv(res,"tcgagtex_deg.csv")

res <- results(dds2, contrast=c("condition","p","m"))
resultsNames(dds2)
write.csv(res,"tcgagtex_degp.csv")


setwd("~/nanjing/gse248562")
#b_vs_a
cts<-read.table("bvsa",head=TRUE) 
#setup DESeqDataSetFromMatrix
cts_round<-round(cts)
countData<-cts_round
condition <- factor(c(rep("b",65),rep("a",65)), levels = c("b","a"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("b","a"))
dds$condition <- relevel(dds$condition, ref = "a")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"b_vs_a_norm.csv")
dds2<-DESeq(dds)
resultsNames(dds2)
res <- results(dds2)
res <- results(dds2, contrast=c("condition","b","a"))
resultsNames(dds2)
write.csv(res,"b_vs_a_deg.csv")

#metastasis_vs_primary
cts<-read.table("metastasis_vs_primary",head=TRUE) 
#setup DESeqDataSetFromMatrix
cts_round<-round(cts)
countData<-cts_round
condition <- factor(c(rep("m",3),rep("p",22)), levels = c("m","p"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
countData[countData==-1]<-0 #Error in DESeqDataSet(se, design = design, ignoreRank) : some values in assay are negative
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("m","p"))
dds$condition <- relevel(dds$condition, ref = "p")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"m_vs_p_norm.csv")
dds2<-DESeq(dds)
resultsNames(dds2)
res <- results(dds2)
res <- results(dds2, contrast=c("condition","m","p"))
resultsNames(dds2)
write.csv(res,"m_vs_p_deg.csv")

#metastasis_vs_primary
cts<-read.table("m_vs_p_14_11",head=TRUE) 
#setup DESeqDataSetFromMatrix
cts_round<-round(cts)
countData<-cts_round
condition <- factor(c(rep("m",14),rep("p",11)), levels = c("m","p"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
countData[countData==-1]<-0 #Error in DESeqDataSet(se, design = design, ignoreRank) : some values in assay are negative
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("m","p"))
dds$condition <- relevel(dds$condition, ref = "p")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"m_vs_p_norm1411.csv")
dds2<-DESeq(dds)
resultsNames(dds2)
res <- results(dds2)
res <- results(dds2, contrast=c("condition","m","p"))
resultsNames(dds2)
write.csv(res,"m_vs_p_deg1411.csv")

#3vs3
cts<-read.table("3vs3",head=TRUE)
cts_round<-round(cts)
countData<-cts_round
condition <- factor(c(rep("m",3),rep("p",3)), levels = c("m","p"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
countData[countData==-1]<-0 #Error in DESeqDataSet(se, design = design, ignoreRank) : some values in assay are negative
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("m","p"))
dds$condition <- relevel(dds$condition, ref = "p")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"3vs3.csv")
dds2<-DESeq(dds)
resultsNames(dds2)
res <- results(dds2)
res <- results(dds2, contrast=c("condition","m","p"))
resultsNames(dds2)
write.csv(res,"3vs3deg.csv")

#8vs8
cts<-read.table("8vs8",head=TRUE)
cts_round<-round(cts)
countData<-cts_round
condition <- factor(c(rep("m",8),rep("p",8)), levels = c("m","p"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
countData[countData==-1]<-0 #Error in DESeqDataSet(se, design = design, ignoreRank) : some values in assay are negative
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("m","p"))
dds$condition <- relevel(dds$condition, ref = "p")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"8vs8.csv")
dds2<-DESeq(dds)
resultsNames(dds2)
res <- results(dds2)
res <- results(dds2, contrast=c("condition","m","p"))
resultsNames(dds2)
write.csv(res,"8vs8deg.csv")

#gse157905
setwd("~/nanjing/gse157905")
cts<-read.table("tmp01",head=TRUE)
cts_round<-round(cts)
countData<-cts_round
condition <- factor(c(rep("m",12),rep("p",4)), levels = c("m","p"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
countData[countData==-1]<-0 #Error in DESeqDataSet(se, design = design, ignoreRank) : some values in assay are negative
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("m","p"))
dds$condition <- relevel(dds$condition, ref = "p")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"gse157905.csv")
dds2<-DESeq(dds)
resultsNames(dds2)
res <- results(dds2)
res <- results(dds2, contrast=c("condition","m","p"))
resultsNames(dds2)
write.csv(res,"157905deg.csv")

#gse117623
setwd("~/nanjing/gse117623")
cts<-read.table("tmp02",head=TRUE)
#cts<-aggregate(a[, -c(1:2)], by=list(DF$EntrezID, DF$Name), mean)
cts_round<-round(cts)
countData<-cts_round
condition <- factor(c(rep("m",52),rep("p",65)), levels = c("m","p"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
countData[countData==-1]<-0 #Error in DESeqDataSet(se, design = design, ignoreRank) : some values in assay are negative
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("m","p"))
dds$condition <- relevel(dds$condition, ref = "p")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"gse117623.csv")
dds2<-DESeq(dds)
resultsNames(dds2)
res <- results(dds2)
res <- results(dds2, contrast=c("condition","m","p"))
resultsNames(dds2)
write.csv(res,"gse117623deg.csv")
















#c_vs_a
cts<-read.table("cvsa",head=TRUE) 
#setup DESeqDataSetFromMatrix
cts_round<-round(cts)
countData<-cts_round
condition <- factor(c(rep("c",6),rep("a",5)), levels = c("c","a"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("c","a"))
dds$condition <- relevel(dds$condition, ref = "a")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"c_vs_a_norm.csv")
dds2<-DESeq(dds)
resultsNames(dds2)
res <- results(dds2)
res <- results(dds2, contrast=c("condition","c","a"))
resultsNames(dds2)
write.csv(res,"c_vs_a_deg.csv")

#d_vs_a
cts<-read.table("dvsa",head=TRUE) 
#setup DESeqDataSetFromMatrix
cts_round<-round(cts)
countData<-cts_round
condition <- factor(c(rep("d",5),rep("a",5)), levels = c("d","a"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("d","a"))
dds$condition <- relevel(dds$condition, ref = "a")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"d_vs_a_norm.csv")
dds2<-DESeq(dds)
resultsNames(dds2)
res <- results(dds2)
res <- results(dds2, contrast=c("condition","d","a"))
resultsNames(dds2)
write.csv(res,"d_vs_a_deg.csv")


#b_vs_c
cts<-read.table("bvsc",head=TRUE) 
#setup DESeqDataSetFromMatrix
cts_round<-round(cts)
countData<-cts_round
condition <- factor(c(rep("b",6),rep("c",6)), levels = c("b","c"))
coldata<-data.frame(row.names=colnames(countData),condition)
condition
coldata
dds<-DESeqDataSetFromMatrix(countData=countData, colData=coldata, design=~condition)
head(dds)
featureData <- data.frame(gene=rownames(cts))
mcols(dds) <- DataFrame(mcols(dds), featureData)
mcols(dds)
keep <- rowSums(counts(dds)) >= 10#keep sum counts >=10
dds <- dds[keep,]
dim(dds)
dds$condition <- factor(dds$condition, levels = c("b","c"))
dds$condition <- relevel(dds$condition, ref = "c")
dds$condition <- droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
dds_norm<-counts(dds,normalized=TRUE)
head(dds_norm)
write.csv(dds_norm,"b_vs_c_norm.csv")
dds2<-DESeq(dds)
resultsNames(dds2)
res <- results(dds2)
res <- results(dds2, contrast=c("condition","b","c"))
resultsNames(dds2)
write.csv(res,"b_vs_c_deg.csv")

#tpm
a<-"counts"
counts01<-data.frame(fread(a),check.names=FALSE,row.names=1)
counts<-as.matrix(counts01)
head(counts)[,1:9]

featureLength01<-read.table("feature_length01",head=TRUE)
featureLength02<-t(featureLength01)
featureLength<-as.vector(featureLength02)


Counts_to_tpm <- function(counts, featureLength) {
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- featureLength
  
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen)
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}

tmp_counts<-Counts_to_tpm(counts,featureLength)
write.table(tpm_counts,"tpm_counts")

#volcano plot delta11_vs_mock_DEG
res <- read.csv("delta11_vs_mock_deg.csv", header=TRUE)
head(res)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-5,7)))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
with(subset(res, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, pvalue<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)
#dev.off()

#volcano plot delta11_vs_wt_DEG
res <- read.csv("delta11_vs_wt_deg.csv", header=TRUE)
head(res)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-5,6)))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
with(subset(res, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, pvalue<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)
#dev.off()

#volcano plot mock_vs_wt_DEG
res <- read.csv("mock_vs_wt_deg.csv", header=TRUE)
head(res)
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-6,7)))
# Add colored points: red if pvalue<0.05, orange of log2FC>1, green if both)
with(subset(res, pvalue<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
with(subset(res, abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="orange"))
with(subset(res, pvalue<.05 & abs(log2FoldChange)>1), points(log2FoldChange, -log10(pvalue), pch=20, col="green"))
#with(subset(res, P.Value<.05 & abs(log2FC)>1), textxy(logFC, -log10(P.Value), labs=Gene, cex=.6))
abline(h=1.3,v=1,lty=3)
abline(v=-1,lty=3)
#dev.off()

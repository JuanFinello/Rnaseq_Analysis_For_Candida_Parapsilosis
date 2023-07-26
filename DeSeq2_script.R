### Difential expresion analisis 

setwd("/media/juan/01D9AF894E0C2760/WORKSHOP/rnaseq_fungi/")
install.packages("DESeq2")
library("DESeq2")

# Import csv file with counts table

Count_data = read.table(file = "merged_counts.csv", header = T, sep = ",",row.names=1,check.names = FALSE)
dim(Count_data)

# Import Table with experiment conditions by column
Col_data = read.table(file = "condition.csv", header = T, sep = ",",row.names = 1)
dim(Col_data)
rownames(Count_data)
colnames(Count_data)
rownames(Col_data)
colnames(Col_data)

# The number of columns in te counnts table and conditions table must be the same
all(rownames(Col_data)==colnames(Count_data))

#boxplot(Count_data)
hist(Count_data[,1]) # Plotting only the first sample (column 1)



##Expresion analisis

#install.packages("DESeq2")
library(DESeq2) # load the DESeq2 package
#count no of NA values in matrix
class(Count_data)
(is.na(Count_data))
which(is.na(Count_data),arr.ind=TRUE)
sum(is.na(Count_data))

#replace missing values in matrix with rowsums
install.packages("zoo")
library(zoo)
Count_data[]<-t(na.aggregate(t(Count_data)))
Count_data

#replace NA values in matrix with zero
Count_data[is.na(Count_data)] <- 0
#removing genes with all zero values
df1[rowSums(df1[])>0,]
?DESeq


dds = DESeqDataSetFromMatrix(countData = round(Count_data),
                             colData = Col_data,
                             design = ~ condition) # we're testing for the different condidtions
dds$condition <- relevel(dds$condition, ref = " planktonic ")
dds
dds <- DESeq(dds)
res1 <- results(dds)
summary(res1)

###keep only sig results, padj<0.05 and log2FoldChange >1
resSigUp <- subset(res1, padj < 0.05 & log2FoldChange >1)
write.csv(resSigUp, "Upregulatedd.csv")

###keep only sig results, padj<0.05 and log2FoldChange < -1
resSigDown <- subset(res1, padj < 0.05 & log2FoldChange < -1)
write.csv(resSigDown, "Downregulated.csv")

###keep UP and Down in one file with padj<0.05
resSig <- subset(res1, log2FoldChange >1 & padj <0.05 | log2FoldChange < -1 & padj < 0.05)
write.csv(resSig, "DE.csv")

##Write for volcano plot
##write.table(res1, file="Volcano", row.names=F, sep="\t")

#Enhanced Volcano
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)

rpl <- read.table("Volcano", header = TRUE, sep = "\t")
x<- res1$log2FoldChange
y<- res1$pvalue

EnhancedVolcano(res1,
                lab = "",
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = NULL,
                title = NULL,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                xlim = c(-8,8),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylim = c(0,12),
                ylab = bquote(~-Log[10]~italic(P)),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                #transcriptPointSize = 0.5,
                #transcriptLabSize = 4.0,
                colAlpha = 1,
                shape = 19,
                subtitle = NULL,
                legendPosition = 'top',
                legendLabSize = 12,
                legendIconSize = 4.0,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                drawConnectors = FALSE,
                widthConnectors = 0.2,
                colConnectors = 'grey50',
                border = 'full' )


# Density plot of count per gene

tpm <- t(t(Count_data)/colSums(Count_data))*1e6
inlog <- log(tpm)
colLabel <- c(rep("#E41A1C", 3), rep("#377EB8", 3))
colTy <- c(rep(1:3, 3), rep(1:3, 3))
plot(density(inlog[,1]), ylim=c(0,0.4), main="Density plot of
counts per gene", lty=colTy[1], xlab="Log of TPM per gene",
     ylab="Density", col=colLabel[1])
for(i in 2:ncol(tpm)){
  lines(density(inlog[,i]), lty=colTy[i], col=colLabel[i])
}
legend("topright", legend=colnames(tpm), lty=colTy, col=
         colLabel)
colLabel

#PCoA plot:

  d <- dist(t(tpm))
fit=cmdscale(d, eig=TRUE, k=2)
x=fit$points[,1]
y=fit$points[,2]
plot(x, y, type="p", pch=20)
text(x, y, labels=row.names(t(tpm)), cex=1, adj=c(-0.25,-0.25))


# PCoA plot:

d <- dist(t(tpm))
fit <- cmdscale(d, eig = TRUE, k = 2)
x <- fit$points[, 1]
y <- fit$points[, 2]

# Create a vector of colors for each sample

sample_colors <- c("#E41A1C", "#E41A1C", "#E41A1C", "#377EB8", "#377EB8", "#377EB8")
#http://127.0.0.1:33937/graphics/plot_zoom_png?width=442&height=482
plot(x, y, type = "p", pch = 20, col = sample_colors)

# Add text labels for the first three samples next to each point
text(x[1:3], y[1:3], labels = row.names(t(tpm))[1:3], cex = 0.6, pos = 4, offset = 0.5)

# Add text labels for the last three samples to the left of each point
text(x[4:6], y[4:6], labels = row.names(t(tpm))[4:6], cex = 0.6, pos = 2, offset = 0.5)

# Add the title
title("Principal Coordinate Analysis (PCoA)")





# Overall Log2 Fold Change

hist(res1[!is.na(res1$padj) & res1$padj <= 0.05,"log2FoldChange"],
     breaks = seq(-15, 15, 0.25),
     xlab = "Log2 Fold Change",
     main = "Overall Log2 Fold Change",
     col=c(rep("tomato", 56), rep("white", 8), rep("tomato", 56)))




# Generating HTML Reports

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install

BiocManager::install("ReportingTools")

library("ReportingTools")

#Use the functions HTMLReport() and publish() to generate a website from the results data.frame created by DESeq2.

htmlRep <- HTMLReport(shortName="P-vs-B_results", reportDirectory = "./reports")

publish(cbind(GeneID=rownames(res1),as.data.frame(res1)),htmlRep)

finish(htmlRep)


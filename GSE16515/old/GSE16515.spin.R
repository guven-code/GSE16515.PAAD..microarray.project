#############################################################################
# Gene expression analysis of expression data from human bone
# marrow hematopoietic stem cells
# Organism	Homo sapiens
# Dataset: GSE32719 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32719)
# Pang WW, Price EA, Sahoo D, Beerman I et al. Human bone marrow hematopoietic stem cells 
# are increased in frequency and myeloid-biased with age.
# Proc Natl Acad Sci U S A 2011 Dec 13;108(50):20012-7.
# PMID: 22123971 (http://www.ncbi.nlm.nih.gov/pubmed/22123971)
# R code: Emine Guven
##############################################################################

library(gplots)
library(biomaRt)
library(matrixStats)
library(Biobase)
library(GEOquery)

# # load series and platform data from GEO
# 
# gset <- getGEO("GSE16515", GSEMatrix =TRUE, getGPL=FALSE)
# if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
# gset <- gset[[idx]]
# 
# # set parameters and draw the plot
# 
# dev.new(width=4+dim(gset)[[2]]/5, height=6)
# par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
# title <- paste ("GSE16515", '/', annotation(gset), " selected samples", sep ='')
# boxplot(exprs(gset), boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
# 
# gse16515exprs<-exprs(gset)

#saveRDS(gse16515exprs, "gset16515.rds")
GSEdata<-readRDS("gset16515.rds")

#geneIDs=GSEdata$X
#conditions=colnames(GSEdata[,-1])

# Check the loaded dataset
dim(GSEdata) # Dimension of the dataset
head(GSEdata) # First few rows
tail(GSEdata) # Last few rows

rownames(GSEdata)
colnames(GSEdata)

#raw_data=data.matrix(GSEdata,rownames.force = NA)
#rownames(raw_data)=geneIDs

pdf("plots/histogramGSE.pdf")
hist(GSEdata, col = "gray", main="GSE16515 - Histogram")
dev.off()


# Check the behavior of logdata
pdf("plots/histogramLogGSE.pdf")
hist(log2(GSEdata), col = "gray", main="GSE16515 - Histogram")
dev.off()

data2=log2(GSEdata)

# Log2 transformation (why?)

###################
# Exploratory plots
###################




pdf("plots/boxPlotRawData.pdf")
boxplot(data2, col=c("red", "red" , "red", "red","blue","red", "blue",
                     "red","blue","red","blue","red","red","red","blue",
                     "red","blue","red","blue","red","red","red","red",
                     "blue","red","red","red","red","red","blue","red",
                     "blue","red","blue","red","red","red","red","red","blue", 
                     "red","red","blue","red","red","red","blue","red",
                     "blue","red","red","blue") , main="GSE16515 - boxplots", las=2,cex.axis=0.6)

legend("center", inset=.02, legend=c("tumor", "normal"),
       fill=c("red", "blue"),cex=0.8,text.font=4)
dev.off()

# Hierarchical clustering of the "samples" based on
# the correlation coefficients of the expression values
hc = hclust(as.dist(1-cor(data2)))
pdf("plots/HclustRawData.pdf")
plot(hc, main="GSE15616 - Hierarchical Clustering",cex=0.6)
dev.off()
#######################################
# Differential expression (DE) analysis
#######################################

# Separate each of the conditions into three smaller data frames
tumor = data2[,c(1:4,6,8,10,12:14,16,18,20:23,25:29,31,33,35:38,40:42,44:46,48,50:51)]
normal = data2[,-c(1:4,6,8,10,12:14,16,18,20:23,25:29,31,33,35:38,40:42,44:46,48,50:51)]


# Compute the means of the samples of each condition
tumor.mean = apply(tumor, 1, mean)
normal.mean = apply(normal, 1, mean)

head(tumor.mean)

head(normal.mean)

# Just get the maximum of all the means
limit = max(tumor.mean, normal.mean)

# Scatter plot
pdf("plots/Scatter.pdf")
plot(normal.mean ~ tumor.mean, xlab = "tumor", ylab = "normal",
     main = "GSE16515 - Scatter", xlim = c(0, limit), ylim = c(0, limit))
# Diagonal line
abline(0, 1, col = "red")
dev.off()

# Compute fold-change (biological significance)
# Difference between the means of the conditions
fold = tumor.mean - normal.mean

# Histogram of the fold differences
pdf("plots/hist.pdf")
hist(fold, col = "gray")
dev.off()
# Compute statistical significance (using t-test)
pvalue = NULL # Empty list for the p-values
tstat = NULL # Empty list of the t test statistics

for(i in 1 : nrow(GSEdata)) { # For each gene : 
  x = tumor[i,] # tumor of gene number i
  y = normal[i,] # normal of gene number i
  
  # Compute t-test between the two conditions
  t = t.test(x, y)
  
  # Put the current p-value in the pvalues list
  pvalue[i] = t$p.value
  # Put the current t-statistic in the tstats list
  tstat[i] = t$statistic
}

head(pvalue)

# Histogram of p-values (-log10)
pdf("plots/pvalueHist.pdf")
hist(-log10(pvalue), col = "gray")
dev.off()
# Volcano: put the biological significance (fold-change)
# and statistical significance (p-value) in one plot
pdf("plots/volcano.pdf")
plot(fold, -log10(pvalue), main = "GSE16515 - Volcano")

fold_cutoff = log2(1)
pvalue_cutoff = 0.01
abline(v = fold_cutoff, col = "blue", lwd = 3)
abline(v = -fold_cutoff, col = "red", lwd = 3)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 3)
dev.off()
# Screen for the genes that satisfy the filtering criteria

# Fold-change filter for "biological" significance
filter_by_fold = abs(fold) >= fold_cutoff
dim(data2[filter_by_fold, ])

# P-value filter for "statistical" significance
filter_by_pvalue = pvalue <= pvalue_cutoff
dim(data2[filter_by_pvalue, ])

# Combined filter (both biological and statistical)
filter_combined = filter_by_fold & filter_by_pvalue

filtered = data2[filter_combined,]
dim(filtered)

head(filtered)


# Let's generate the volcano plot again,
# highlighting the significantly differential expressed genes
pdf("plots/DEGs.pdf")
plot(fold, -log10(pvalue), main = "GSE16515 - Volcano #2")
points (fold[filter_combined], -log10(pvalue[filter_combined]),
        pch = 16, col = "red")
dev.off()

# Highlighting up-regulated in red and down-regulated in blue
pdf("plots/up&down.pdf")
plot(fold, -log10(pvalue), main = "GSE16515 - Volcano #3")
points (fold[filter_combined & fold < 0],
        -log10(pvalue[filter_combined & fold < 0]),
        pch = 16, col = "red")
points (fold[filter_combined & fold > 0],
        -log10(pvalue[filter_combined & fold > 0]),
        pch = 16, col = "blue")
dev.off()

# Cluster the rows (genes) & columns (samples) by correlation
rowv = as.dendrogram(hclust(as.dist(1-cor(t(filtered)))))
colv = as.dendrogram(hclust(as.dist(1-cor(filtered))))

# Generate a heatmap
pdf("plots/heatmap.pdf")
heatmap(filtered, Rowv=rowv, Colv=colv, cexCol=0.7)
dev.off()
# install.packages("gplots")		# Uncomment if not already installed
# install.packages("RColorBrewer")	# Uncomment if not already installed

library(gplots)

# Enhanced heatmap
pdf("plots/heatmap2.pdf")
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
          col = rev(redblue(256)), scale = "row")
dev.off()
# Save the heatmap to a PDF file
pdf ("plots/GSE16515_DE_Heatmap.pdf")
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
          col = rev(redblue(256)), scale = "row")
dev.off()

# Save the DE genes to a text file
write.table (filtered, "GSE16515_DE.txt", sep = "\t",
             quote = FALSE)

#############################################################################

n = nrow(filtered)

cor.table = NULL
x = NULL
y = NULL
cor.val = NULL
cor.sig = NULL

for (i in 1 : (n-1)) {
  x_name = rownames(filtered)[i]
  x_exps = filtered[i, ]	
  
  for (j in (i+1) : n) {
    y_name = rownames(filtered)[j]
    y_exps = filtered[j, ]
    
    output = cor.test(x_exps,y_exps)
    
    x = c(x, x_name)
    y = c(y, y_name)
    cor.val = c(cor.val, output$estimate)
    cor.sig = c(cor.sig, output$p.value)
  }
}

cor.table = data.frame (x, y, cor.val, cor.sig)

dim(cor.table)
head(cor.table)

sig_cutoff = 0.001

cor.filtered = subset (cor.table, cor.sig < sig_cutoff)

dim(cor.filtered)
head(cor.filtered)
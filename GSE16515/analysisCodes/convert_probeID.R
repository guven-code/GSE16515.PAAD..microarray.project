source("http://bioconductor.org/biocLite.R")
biocLite("GEOquery")
biocLite("limma")
biocLite("Biobase")
biocLite("affy")

library(GEOquery)
library(biomaRt)


gset <- getGEO("GSE16515", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

class(gset)

length(gset)

slotNames(gset)

dim(gset)

str(gset)

#require("biomaRt")
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("hsapiens_gene_ensembl", mart)
annotLookup <- getBM(
  mart=mart,
  attributes=c(
    "affy_hg_u133_plus_2",
    "ensembl_gene_id",
    "gene_biotype",
    "external_gene_name"),
  filter = "affy_hg_u133_plus_2",
  values = rownames(exprs(gset)), uniqueRows=TRUE)

indicesLookup <- match(rownames(gset), annotLookup$affy_hg_u133_plus_2)
rownames(gset) <- paste(annotLookup[indicesLookup, "external_gene_name"], c(1:length(indicesLookup)), sep="_")

head(rownames(gset),20)
extGeneNames<-gsub("_[0-9]*$", "", rownames(gset))
gene.exp<-exprs(gset)

rownames(gene.exp)<-extGeneNames

write.csv(gene.exp,"GSE16515_externalGeneNames_mat.csv")




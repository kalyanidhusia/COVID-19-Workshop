##############################################################################
#                                                                            #   
# DIFFERENTIAL GENE EXPRESSION OF SARS-COV-2 INFECTION OF PRIMARY HUMAN LUNG #
# EPITHELIUM FOR COVID 19.                                                   #
#                                                                            #
#  RNA-SEQ DATA FROM  GES160435                                              #
#  https://www.ncbi.nlm.nih.gov//geo/query/acc.cgi?acc=GSE160435             #
#                                                                            #
#  NUMBER OF MOCK SAMPLES = 5                                                #                           
#   NUMBER OF INFECTED SAMPLES 5                                             #
#   TOTAL NUMBER OF SAMPLES = 10                                             #
#                                                                            #
#                                                                            #
#                                                                            #
# AUTHOR = DANILE K MARRI                                                    #                   
#                                                                            #
#                                                                            #
#                                                                            #  
#                                                                            #
##############################################################################

#SET WORKING DIRECTORY (CHANGE THIS TO THE DIRECTORY YOU HAVE THE RNASEQ DATASET )

setwd("C:/Users/marri/OneDrive/Desktop/COVID 19 WORKSHOP FOLDER")


#IMPORTING THE DATASET

RNASEQ_DATA <- read.table("GSE160435_TPM.txt", header=TRUE, row.names=1)

# CONVERT DATAFRAME TO MATRIX
RNASEQ_DATA <- as.matrix(RNASEQ_DATA)
head(RNASEQ_DATA)

#PLOTTING THE DATASET
#barplot(colSums(RNASEQ_DATA)/1e5, las=3)

plot(log(RNASEQ_DATA[,5]),log(RNASEQ_DATA[,10]))


# ASSIGNING CONDIRIONS TO THE DATASET (FIRST FIVE COV INFECTED, REMAINING FIVE MOCK DATASET)
(condition <- factor(c(rep("exp", 5), rep("ctl", 5))))


# ANALYSIS WITH DESeq2

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)
library(RColorBrewer)
library(ggplot2)
library(dplyr)


# Create a coldata frame and instantiate the DESeqDataSet.
(coldata <- data.frame(row.names=colnames(RNASEQ_DATA), condition))


#CReating the DESeq Dataset. 
dds <- DESeqDataSetFromMatrix(countData=round(RNASEQ_DATA), colData=coldata, design=~condition)
dds

# Run the DESeq pipeline (Creating the DESeq object)
dds <- DESeq(dds)

#checking our dataset
nrow(dds)


#Remove genes with count transcript <3

dds = dds[rowSums(counts(dds))>3,]


nrow(dds)

sizeFactors(dds)


# Plot dispersions
png("qc-dispersions.png")
plotDispEsts(dds, main="Dispersion plot")
dev.off()

plotDispEsts(dds, main="Dispersion plot")


#PRINCIPAL COMPONENTS PLOT
rld = rlog(dds)
png("PCA_PLOT.png")
plotPCA(rld, intgroup = "condition")
dev.off()

plotPCA(rld, intgroup = "condition")

# Function to detect the sample groups:

detectGroups <- function (x){  # x are col names
  tem <- gsub("[0-9]*$","",x) # Remove all numbers from end
  #tem = gsub("_Rep|_rep|_REP","",tem)
  tem <- gsub("_$","",tem); # remove "_" from end
  tem <- gsub("_Rep$","",tem); # remove "_Rep" from end
  tem <- gsub("_rep$","",tem); # remove "_rep" from end
  tem <- gsub("_REP$","",tem)  # remove "_REP" from end
  return( tem )
}

#detectGroups(colnames(RNASEQ_DATA))


#Function for calculating distance between samples for the heatmap:
# distance function = 1-PCC (Pearson's correlation coefficient)
  
dist2 <- function(x, ...)as.dist(1-cor(t(x), method="pearson"))


hclust2 <- function(x, method="average", ...)  # average linkage in hierarchical clustering
  hclust(x, method=method, ...)

n=100 # number of top genes by standard deviation

x = assay(rld)
if(n>dim(x)[1]) n = dim(x)[1] # max	as data

x = x[order(apply(x,1,sd),decreasing=TRUE),]  # sort genes by standard deviation

x = x[1:n,]   # only keep the n genes

# this will cutoff very large values, which could skew the color 
x=as.matrix(x[1:n,])-apply(x[1:n,],1,mean)
cutoff = median(unlist(x)) + 4*sd (unlist(x)) 
x[x>cutoff] <- cutoff
cutoff = median(unlist(x)) - 4*sd (unlist(x)) 
x[x< cutoff] <- cutoff

groups = detectGroups(colnames(x) )
groups.colors = rainbow(length(unique(groups) ) )


lmat = rbind(c(5,4),c(0,1),c(3,2))
lwid = c(1.5,4)
lhei = c(1,.2,4)

png("RNASEQ_HEATMAP.png")
heatmap.2(x, distfun = dist2,hclustfun=hclust2,
          col=greenred(75), density.info="none", trace="none", scale="none", keysize=.5
          ,key=T, symkey=F
          ,ColSideColors=groups.colors[ as.factor(groups)]
          ,margins=c(8,12)
          ,cexRow=1
          ,srtCol=45
          ,cexCol=1.  # size of font for sample names
          ,lmat = lmat, lwid = lwid, lhei = lhei
)

dev.off()

heatmap.2(x, distfun = dist2,hclustfun=hclust2,
          col=greenred(75), density.info="none", trace="none", scale="none", keysize=.5
          ,key=T, symkey=F
          ,ColSideColors=groups.colors[ as.factor(groups)]
          ,margins=c(8,12)
          ,cexRow=1
          ,srtCol=45
          ,cexCol=1.  # size of font for sample names
          ,lmat = lmat, lwid = lwid, lhei = lhei
)

#DIFFERENTIAL GENE EXPRESSION OF THE DATASET.


#GETTING THE DIFFERENTIAL EXPRESSION RESULTS..
# WE CAN USE resultNames function to get the name of the DGE results

DGE_res = results(dds)
head(DGE_res)
summary(DGE_res)

#CUSTOMIZING THE LFC(LOG FOLD CHANGE) with a threshold
DGE_ress = results(dds, lfcThreshold = 0.01)
head(DGE_ress)
summary(DGE_ress)


#MA plot to visualize to visualize the significant of the various genes.

png("RNASEQ_MAPLOT.png")
DESeq2::plotMA(DGE_ress, ylim=c(-5,5))
dev.off()

#PLOT THE VOLCANO PLOT

res = as.data.frame(DGE_ress)

#Mutate to add signicance column to the dataframe base on the adjusted p-value

res_mut = mutate(res, significant=ifelse(res$padj<0.1, "Sig","Not_Sig"))
res_mut[which(abs(res$log2FoldChange) < 1.0), "significant"]= "Not_Sig"

#plot the results


png("RNASEQ_SIG_NON_SIG_GENES.png")
ggplot(res_mut, aes(log2FoldChange, -log10(padj))) + geom_point(aes(col = significant))+scale_color_manual(values= c("red","black"))
dev.off()

# SORTING OUT THE HIGHEST LOG FOLD GENES

res_mut = res_mut[order(abs(res_mut$log2FoldChange), decreasing = TRUE), ]
topgene = rownames(res_mut)[1]

#plotcounts for the top gene

plotCounts(dds, gene = topgene, intgroup = "condition")

# Write the dataframe as a csv file

write.csv(res_mut, file = "DGE_results.csv")
























































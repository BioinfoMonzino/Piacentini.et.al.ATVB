##########################################################################################################################################
###                                                                                                                                    ###
###      Whole-blood transcriptome unveils altered immune response in acute myocardial infarction patients with aortic valve sclerosis ###
###                                                                                                                                    ###
##########################################################################################################################################

## Clean Global Environment:
rm(list=ls())

## Load libraries:
library(EDASeq)
library(RUVSeq)
library(edgeR)
library(RColorBrewer)
library(DESeq2)
library(DaMiRseq)

## set working directory:
# setwd("<path to dir>")

## Load matrix of raw counts:
Raw_Counts <- read.delim("<path to>/Raw_Counts.txt") # change path
# coerce to matrix:
Raw_Counts<-as.matrix(Raw_Counts)

## Load metadata:
metadata <- read.delim("<path to>/metadata.txt", stringsAsFactors=TRUE) # change path

################################################################################

###################
## RNA-Seq power ##
###################

library(RNASeqPower)

## calculate power for actual sample size:
rnapower(depth=10, 
         n=66, n2=44,
         cv=0.2, # estimated from previous work
         effect=1.4,
         alpha= .001)

# depth=c(2000, 245, 37, 3.7, 0.3) # c(1000, 100, 27, 10, 1)
# power=c(.7, .8, .9)
# effect=c(1,1.25, 1.5, 2, 3)

################################################################################

################################################################################


######################################################
##     Differential Expression Analysis (DEA)       ##
######################################################

# filter count matrix:
retained_genes<- rowSums(Raw_Counts>10) >= 21
Raw_Counts_filt<-as.matrix(Raw_Counts[retained_genes,])
# check matrix:
dim(Raw_Counts_filt)
head(Raw_Counts_filt, n=10)

# create 'SeqExpressionSet' object with raw counts and metadata:
set <- newSeqExpressionSet(Raw_Counts_filt, phenoData = metadata, row.names=colnames(Raw_Counts_filt))
set
# normalize between samples:
set <- betweenLaneNormalization(set, which="upper", offset=TRUE)

# First pass DEA by DESeq2 (with AMI+AVSc):
dds <- DESeqDataSetFromMatrix(countData = counts(set),
                              colData = pData(set),
                              design = ~AMI+AVSc)

# test the full model against a null model:
dds <- DESeq(dds, test="LRT", reduced=as.formula("~1"))
res <- results(dds)
# order by pvalues:
res<-res[order(res$pvalue),]
res

# set empirical genes:
empirical <- rownames(res)[which(!(rownames(res) %in% rownames(res)[1:15000]))]
empirical

# performs factor analysis of the read counts by RUVg method; test different number of k to be estimated and see differences:
setk1 <- RUVg(set, empirical, k=1)
setk5 <- RUVg(set, empirical, k=5)
setk13 <- RUVg(set, empirical, k=13)
setk15 <- RUVg(set, empirical, k=15)
setk20 <- RUVg(set, empirical, k=20) # max number of k tested, include W_1 to W_20 "latent" variables

# # Look at RLE and PCA plots after data adjustement for a different number of k ("latent" variables); example with k=13
# # set colors for RLE and PCA plots:
# colors <- brewer.pal(4, "Set2")
# # plot RLE (change first argument to test different "k" adjustment):
# plotRLE(setk13@assayData$normalizedCounts, outline=FALSE, col=colors[metadata$AVSc], main="k13", las=2, cex.axis=0.25)
# # plot PCA (change first argument to test different "k" adjustment):
# plotPCA(setk13@assayData$normalizedCounts, outline=FALSE, col=colors[metadata$AVSc], main="k13")

## Second pass DEA with adjustment for latent variables:
# running analysis for Mod2 (adjustment includes 13 "latent" variables (W_1 to W_13), AMI, age and Previous AMI/PCI/CABG); AVSc is the variable of interest; 
# change design according to different models to be tested.
dds <- DESeqDataSetFromMatrix(countData = counts(setk20),
                              colData = pData(setk20),
                              design = ~W_1+W_2+W_3+W_4+W_5+W_6+W_7+W_8+W_9+W_10+W_11+W_12+W_13+age+previous+AMI+AVSc) 

dds <- DESeq(dds)
res <- results(dds, contrast=c("AVSc","#1","#0"))
res<-res[order(res$pvalue),]
res

# plot histogram of pvalues:
hist(res[,5], main="AVSc vs. no-AVSC", breaks=20, col="green", xlab="pval") 

# MA-plot from base means and log fold changes:
plotMA(res, ylim=c(-2.5,2.5))

################################################################################

################################################################################

#################################
##     Feature selection       ##
#################################

## Extract normalized counts:
Norm_Counts<-setk13@assayData$normalizedCounts

# library (DaMiRseq) 

# add class (=AVSc) column label to metadata as per package requirement:
metadata$class<-metadata$AVSc

# proceed by creating a 'SummarizedExperiment' object and log transform normalized counts:
SE<-DaMiR.makeSE(Norm_Counts, metadata)
data_norm <- DaMiR.normalization(SE, minCounts=0, fSample=0.1, type="vst",
                                 hyper = "no")


## Load library:
library(GARS)

# set.seed to reproduce result:
set.seed(1)
# run GARS to select 100 informative features (or a different number of features changing chr.len):
# Note: performing this step with the following setting is quite time consuming (~1 hour) for a common workstation, e.g. with an esa core CPU (2.40 GHz, 16 GB RAM) 
res_GA <- GARS_GA(data=data_norm,
                  classes = colData(data_norm),
                  chr.num = 1000,
                  chr.len = 100,
                  generat = 500,
                  co.rate = 0.8,
                  mut.rate = 0.01,
                  n.elit = 10,
                  type.sel = "RW",
                  type.co = "one.p",
                  type.one.p.co = "I.quart",
                  n.gen.conv = 150,
                  plots="yes",
                  verbose="yes")

# plot fitness scores vs. generations:
fitness_scores <- FitScore(res_GA)
GARS_PlotFitnessEvolution(fitness_scores)

# bubble chart to assess the usage of n features across generations:
Allpopulations <- AllPop(res_GA)
GARS_PlotFeaturesUsage(Allpopulations,
                       rownames(data_norm),
                       nFeat = 10)

# Extract the expression matrix of features selected:
fs_GARS <- MatrixFeatures(res_GA)

# Visualize feature selected data by MDS and hierarchical clustering:
DaMiR.MDSplot(fs_GARS, as.data.frame(colData(data_norm)))
DaMiR.Clustplot(fs_GARS, as.data.frame(colData(data_norm)[, c(2,4,5,6)])) # class=AVSc

# # Run the following chunck if needed:
# # replace features characters that may affect computation in the following R codes:
# colnames(fs_GARS)<-gsub("[-]", "_", colnames(fs_GARS))
# colnames(fs_GARS)<-gsub("[.]", "_", colnames(fs_GARS))
# # inspect for other possible character replacement:
# colnames(fs_GARS)


## Testing prediction performance of the selected features:

# library(DaMiRseq)

# set.seed to reproduce result:
set.seed(12345)
# perform prediction (generate prediction metrics and plots for the chosen classifiers and the Ensemble):
Classification.res <- DaMiR.EnsembleLearning(as.data.frame(fs_GARS),
                                             as.factor(metadata$class),
                                             fSample.tr = 0.7,
                                             fSample.tr.w = 0.7, 
                                             iter = 500, 
                                             cl_type = c("RF", "SVM", "PLS")
)


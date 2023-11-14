##########################################################################################################################################
###                                                                                                                                    ###
###      Whole-blood transcriptome unveils altered immune response in acute myocardial infarction patients with aortic valve sclerosis ###
###                                                                                                                                    ###
##########################################################################################################################################
## Clean  Environment:
rm(list=ls())
graphics.off()

## Load libraries:
library(EDASeq)
library(RUVSeq)
library(edgeR)
library(RColorBrewer)
library(DESeq2)
library(DaMiRseq)
library(GARS)
library(RNASeqPower)
## set working directory:
# setwd("<path to dir>")

################################################################################
###################
## RNA-Seq power ##
###################
## calculate power for actual sample size:
rnapower(depth=10, # number of reads to call a gene expressed
         n=66, # number of samples in the first group
         n2=44,# number of samples in the second group
         cv=0.2, # estimated from previous work
         effect=1.4, # effect size, i.e., the fold change (linear scale)
         alpha= .001)

################################################################################
###################
##  Data Import  ##
###################
## Load matrix of raw counts and convert into matrix:
Raw_Counts <- read.delim("<path to>/Raw_Counts.txt") # change path
Raw_Counts <- as.matrix(Raw_Counts)

## Load metadata:
metadata <- read.delim("<path to>/metadata.txt", stringsAsFactors = TRUE, na.strings = "") # change path

################################################################################
######################################################
##     Differential Expression Analysis (DEA)       ##
######################################################
# filter count matrix (keep genes whose count is greater than 10 reads in at least 21 samples):
retained_genes <- rowSums(Raw_Counts>10) >= 21
Raw_Counts_filt <-as.matrix(Raw_Counts[retained_genes,])

# create 'SeqExpressionSet' object with raw counts and metadata:
set <- newSeqExpressionSet(counts = Raw_Counts_filt, 
                           phenoData = metadata,
                           row.names = colnames(Raw_Counts_filt))

# normalize between samples:
set <- betweenLaneNormalization(x = set, which="upper", offset=TRUE)

######################## First pass DEA 
# DEA by DESeq2 (with AMI+AVSc):
dds <- DESeqDataSetFromMatrix(countData = counts(set),
                              colData = pData(set),
                              design = ~AMI+AVSc)

# test the full model against a null model:
dds <- DESeq(object = dds, test="LRT", reduced = as.formula("~1"))
res <- results(dds)
res <- res[order(res$pvalue),]

# set empirical genes, to estimate the latent variables:
empirical <- rownames(res)[which(!(rownames(res) %in% rownames(res)[1:15000]))]

# identify latent features by RUVg method, trying different k values
setk1 <- RUVg(set, empirical, k=1)
setk5 <- RUVg(set, empirical, k=5)
setk13 <- RUVg(set, empirical, k=13)
setk15 <- RUVg(set, empirical, k=15)
setk20 <- RUVg(set, empirical, k=20) 

# Diagnostic plots (example with k=13):
colors <- brewer.pal(4, "Set2")
# plot RLE
plotRLE(setk13@assayData$normalizedCounts, outline=FALSE, col=colors[metadata$AVSc], main="k13", las=2, cex.axis=0.25)
# plot PCA 
plotPCA(setk13@assayData$normalizedCounts, outline=FALSE, col=colors[metadata$AVSc], main="k13")

###################### Second pass DEA with adjustment for latent variables
# running analysis for Mod2 (adjustment includes 13 "latent" variables (W_1 to W_13), AMI, age and Previous AMI/PCI/CABG); AVSc is the variable of interest:
dds <- DESeqDataSetFromMatrix(countData = counts(setk20),
                              colData = pData(setk20),
                              design = ~W_1+W_2+W_3+W_4+W_5+W_6+W_7+W_8+W_9+W_10+W_11+W_12+W_13+age+previous+AMI+AVSc) 

dds <- DESeq(dds)
res <- results(dds, contrast=c("AVSc","#1","#0"))
res <- res[order(res$pvalue),]

# plot histogram of pvalues:
hist(res[,5], main="AVSc vs. no-AVSC", breaks=20, col="green", xlab="pval") 

# MA-plot:
plotMA(res, ylim=c(-2.5,2.5))

################################################################################
#################################
##     Feature selection       ##
#################################
## Extract normalized counts:
Norm_Counts<-setk13@assayData$normalizedCounts

# add 'class' column label to metadata as per package requirement:
metadata$class<-metadata$AVSc

# creating a 'SummarizedExperiment' object and log transform normalized counts:
SE<-DaMiR.makeSE(Norm_Counts, metadata)
data_norm <- DaMiR.normalization(SE, 
                                 minCounts=0, 
                                 fSample=0.1, 
                                 type="vst",
                                 hyper = "no")

# set.seed to reproduce result:
set.seed(1)

# Implement a Genetic Algorithm to select 100 informative features:
res_GA <- GARS_GA(data = data_norm,
                  classes = colData(data_norm),
                  chr.num = 1000, 
                  chr.len = 100, # number of expected informative genes
                  generat = 500, 
                  co.rate = 0.8, 
                  mut.rate = 0.01,
                  n.elit = 10,
                  type.sel = "RW",
                  type.co = "one.p",
                  type.one.p.co = "I.quart",
                  n.gen.conv = 150,
                  plots = "yes",
                  verbose = "no")

# Diagnostic plot (fitness evolution over generations):
fitness_scores <- FitScore(res_GA)
GARS_PlotFitnessEvolution(fitness_scores)

# bubble chart to assess the usage of the first 10 features across generations:
Allpopulations <- AllPop(res_GA)
GARS_PlotFeaturesUsage(Allpopulations,
                       rownames(data_norm),
                       nFeat = 10)

# Extract the expression matrix of features selected:
fs_GARS <- MatrixFeatures(res_GA)

# Visualize feature selected data by MDS and hierarchical clustering:
DaMiR.MDSplot(fs_GARS, as.data.frame(colData(data_norm)))
DaMiR.Clustplot(fs_GARS, as.data.frame(colData(data_norm)[, c(2,4,5,6)])) 

# # Run the following chunck if needed:
# # replace features characters that may affect computation in the following R codes:
# colnames(fs_GARS)<-gsub("[-]", "_", colnames(fs_GARS))
# colnames(fs_GARS)<-gsub("[.]", "_", colnames(fs_GARS))
# # inspect for other possible character replacement:
# colnames(fs_GARS)

#################################################################################
#################################################################################
##     Assess features discrimination power by Ensemble Machine Learning       ##
#################################################################################
# set.seed to reproduce result:
set.seed(12345)
# perform prediction:
Classification.res <- DaMiR.EnsembleLearning(as.data.frame(fs_GARS),
                                             as.factor(metadata$class),
                                             fSample.tr = 0.7,
                                             fSample.tr.w = 0.7, 
                                             iter = 500, 
                                             cl_type = c("RF", "SVM", "PLS")
)

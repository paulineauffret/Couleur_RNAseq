###############################################################################
### R script to compare several conditions with the SARTools and DESeq2 packages
### Hugo Varet
### Dec 11th, 2017
### designed to be executed with SARTools 1.6.0
################################################################################
#Last modified 22 may 2018

#Chargement des librairies
source("http://bioconductor.org/biocLite.R")
library(RSvgDevice)
library(DESeq2)
library(tximport)
library(ggplot2)
library("pheatmap")
library('RColorBrewer')
library("assertthat")
library("scales")
library(genefilter)
library(ggrepel)

library(devtools)
#biocLite(c("limma", "edgeR", "genefilter"))
#install_github("PF2-pasteur-fr/SARTools", build_vignettes=TRUE)
vignette("SARTools")

###################################################################################################
#
#     Set parameters
#
###################################################################################################
curdir=""   #working directory
sampleSheet=paste(curdir,"/SampleSheet.txt",sep="")

#Imports samples infos
samplesInfo <- read.table(sampleSheet, h=T) ; head(samplesInfo) ; dim(samplesInfo)
target_file <- samplesInfo[,c(3,1,5)]
colnames(target_file) <- c("label","files","color") ; head(target_file) ; dim(target_file)
write.table(target_file,paste(curdir,"/targetfile.txt",sep=""),row.name=F,sep="\t",quote=F)

################################################################################
###                parameters: to be modified by the user                    ###
################################################################################
workDir <- ""       #output directory

projectName <- "Pinctada margaritifera Flesh Color"          # name of the project

author <- "Pauline Auffret"                                             # author of the statistical analysis/report

targetFile <- paste(curdir,"/targetfile_f1.txt",sep="")                       # path to the design/target file

rawDir <- ""    # path to the directory containing raw counts files


#featuresToRemove <- c("alignment_not_unique",        # names of the features to be removed
#"ambiguous", "no_feature",     # (specific HTSeq-count information and rRNA for example)
#"not_aligned", "too_low_aQual")# NULL if no feature to remove
featuresToRemove <- NULL

varInt <- "color"                                    # factor of interest
condRef <- "Black"                                   # reference biological condition
batch <- NULL                                        # blocking factor: NULL (default) or "batch" for example

fitType <- "parametric"                              # mean-variance relationship: "parametric" (default) or "local"
cooksCutoff <- TRUE                                  # TRUE/FALSE to perform the outliers detection (default is TRUE)
independentFiltering <- TRUE                         # TRUE/FALSE to perform independent filtering (default is TRUE)
alpha <- 0.01                                        # threshold of statistical significance
pAdjustMethod <- "BH"                                # p-value adjustment method: "BH" (default) or "BY"

typeTrans <- "rlog"                                  # transformation for PCA/clustering: "VST" or "rlog"
locfunc <- "median"                                  # "median" (default) or "shorth" to estimate the size factors

colors <- c("grey20","antiquewhite3","darkorange2")  # vector of colors of each biological condition on the plots

#Function to show colors
showCols <- function(cl=colors(), bg = "grey",
                     cex = 0.75, rot = 30) {
  m <- ceiling(sqrt(n <-length(cl)))
  length(cl) <- m*m; cm <- matrix(cl, m)
  require("grid")
  grid.newpage(); vp <- viewport(w = .92, h = .92)
  grid.rect(gp=gpar(fill=bg))
  grid.text(cm, x = col(cm)/m, y = rev(row(cm))/m, rot = rot,
            vp=vp, gp=gpar(cex = cex, col = cm))
}
showCols(bg="gray20",cl=colors()[1:100], rot=30, cex=0.9)


################################################################################
###                             running script                               ###
################################################################################
setwd(workDir)
library(SARTools)

# checking parameters
checkParameters.DESeq2(projectName=projectName,author=author,targetFile=targetFile,
                       rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                       condRef=condRef,batch=batch,fitType=fitType,cooksCutoff=cooksCutoff,
                       independentFiltering=independentFiltering,alpha=alpha,pAdjustMethod=pAdjustMethod,
                       typeTrans=typeTrans,locfunc=locfunc,colors=colors)

# loading target file
target <- loadTargetFile(targetFile=targetFile, varInt=varInt, condRef=condRef, batch=batch)

# loading counts
counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove)

# description plots
majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

# analysis with DESeq2
out.DESeq2 <- run.DESeq2(counts=counts, target=target, varInt=varInt, batch=batch,
                         locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                         cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)

# PCA + clustering
exploreCounts(object=out.DESeq2$dds, group=target[,varInt], typeTrans=typeTrans, col=colors)

# summary of the analysis (boxplots, dispersions, diag size factors, export table, nDiffTotal, histograms, MA plot)
summaryResults <- summarizeResults.DESeq2(out.DESeq2, group=target[,varInt], col=colors,
                                          independentFiltering=independentFiltering,
                                          cooksCutoff=cooksCutoff, alpha=alpha)

# save image of the R session
save.image(file=paste0(projectName, ".RData"))

# generating HTML report
writeReport.DESeq2(target=target, counts=counts, out.DESeq2=out.DESeq2, summaryResults=summaryResults,
                   majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
                   targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
                   condRef=condRef, batch=batch, fitType=fitType, cooksCutoff=cooksCutoff,
                   independentFiltering=independentFiltering, alpha=alpha, pAdjustMethod=pAdjustMethod,
                   typeTrans=typeTrans, locfunc=locfunc, colors=colors)


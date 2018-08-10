#!/usr/bin/Rscript
#Differential Analysis Flesh RNA-seq 14 nov 2017 : functions file
#see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf


#Custom plotPCA function to plot PC1 et PC3
plotPCA.san <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  
  d <- data.frame(PC1 = pca$x[, 1], PC3 = pca$x[, 3], group = group, 
                  intgroup.df, name = colData(rld)[,1])
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:3]
    return(d)
  }
  ggplot(data = d, aes_string(x = "PC1", y = "PC3", color = "group", label = "name")) + geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + ylab(paste0("PC3: ", round(percentVar[3] * 100), "% variance")) + coord_fixed() + geom_text_repel(size=3) 
  
}

#Heatmap
#Little trick to modify xlab orientation
draw_colnames_45 <- function (coln, ...) {
  m = length(coln)
  x = (1:m)/m - 1/2/m
  grid.text(coln, x = x, y = unit(0.96, "npc"), vjust = .5, 
            hjust = 1, rot = 45, gp = gpar(...)) ## Was 'hjust=0' and 'rot=270'
}

## For pheatmap_1.0.8 and later:
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}

#MAplots
#defining a new function to plot all genes and not only log2FC > 2
plotMA.DESeqResults <- function(object, alpha, main="", xlab="mean of normalized counts", ylim, MLE=FALSE, ...) {
  if (missing(alpha)) {
    alpha <- if (is.null(metadata(object)$alpha)) {
      0.1
    } else {
      metadata(object)$alpha
    }
  }
  df <- if (MLE) {
    # test if MLE is there
    if (is.null(object$lfcMLE)) {
      stop("lfcMLE column is not present: you should first run results() with addMLE=TRUE")
    }
    data.frame(mean = object$baseMean,
               lfc = object$lfcMLE,
               isDE = ifelse(is.na(object$padj), FALSE, object$padj < alpha))
  } else {
    data.frame(mean = object$baseMean,
               lfc = object$log2FoldChange,
               isDE = ifelse(is.na(object$padj), FALSE, object$padj < alpha))
    #isDE = ifelse(object$padj < alpha & abs(object$log2FoldChange) >= 2, TRUE, FALSE))
    
  }
  print(dim(df[df$isDE==TRUE,]))
  if (missing(ylim)) {
    plotMA(df, main=main, xlab=xlab, ...)
  } else {
    plotMA(df, main=main, xlab=xlab, ylim=ylim, ...)
  }  
}
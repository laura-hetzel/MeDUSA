#' Volcano plot
#'
#' @param data dataframe result from statistical tests within the sumR package with significant column (true or false)
#' @param xvalues log2 folchange between two groups column
#' @param yvalues nominal p value column
#' @param title title of the plot
#'
#' @return volcano plot
#' @import ggplot2
#' @importFrom ggpubr theme_pubr labs_pubr
#' @export
volcanoPlot <- function(exp, test, title = "") {
  if (!validateExperiment(exp)) return(NULL)

  data <- as.data.frame(rowData(exp)[[test]])
  foldChanges <- rowData(exp)$foldChange[,4]
  pvalues <- data$p.value
  ggplot(data) +
    geom_point(aes(x = foldChanges, y = -log10(pvalues), colour = significant)) + ## color by significant of fdr <0.1
    ggtitle(title) +
    xlab("log2 fold change") +
    ylab("-log10 nominal p.value") +
    # scale_y_continuous(limits = c(0,50)) +
    theme(
      legend.position = "right",
      plot.title = element_text(size = rel(1.5), hjust = 0.5),
      axis.title = element_text(size = rel(1.25))
    ) +
    theme_pubr() +
    labs_pubr()
}

#' @title Plot UMAP
#' @param exp
#' @param assay
#' @param doPCA
#' @param components
#' @importFrom umap umap
#' @importFrom stats prcomp
#' @importFrom ggrepel geom_text_repel
#' @export
plotUMAP <- function(exp, assay = 1, components = 20){
  if (!validateExperiment(exp)) return(NULL)

  data <- as.matrix(t(assay(exp, assay)))
  data <- prcomp(data, scale. = T, center = T)$x[,1:components]

  um <- umap(data)$layout
  colnames(um) <- c("UMAP_1", "UMAP_2")
  um <- as.data.frame(cbind(um, colData(exp)))
  suppressWarnings(
    ggplot(um, aes(x = UMAP_1, y = UMAP_2,
                      label = rownames(colData(exp)),
                      color = .data[[metadata(exp)$phenotype]])) +
    geom_point() +
    geom_text_repel() +
    ggtitle("Umap") +
    theme_bw()
  )
}

plotCellSds <- function(exp, assay = 1, top = nrow(exp)){
  if (!validateExperiment(exp)) return(NULL)

  vars <- colSds(as.matrix(assay(exp, assay)), na.rm = T)
  df <- data.frame(Cell = colnames(exp), Variation = vars)
  df <- df[order(df$Variation, decreasing = T), ]
  df <- df[1:top, ]
  ggplot(df, aes(x = reorder(Cell, -Variation), y = Variation)) +
    geom_bar(stat = "identity") +
    xlab("Cell") +
    ylab("Standard Deviation") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

#' @title Plot the standard deviations
#' @export
plotFeatureSds <- function(exp, assay = 1, top = nrow(exp)){
  if (!validateExperiment(exp)) return(NULL)

  vars <- rowSds(as.matrix(assay(exp, assay)), na.rm = T)
  df <- data.frame(Feature = rownames(exp), Variation = vars)
  df <- df[order(df$Variation, decreasing = T), ]
  df <- df[1:top, ]
  ggplot(df, aes(x = reorder(Feature, -Variation), y = Variation)) +
    geom_bar(stat = "identity") +
    xlab("Feature") +
    ylab("Standard Deviation") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

#' @title Plot Accuracy of Cross Validation
#' @param exp SummarizedExperiment object
#' @param modelName Either name of number of the model
#' @export
plotCrossValidation <- function(exp, modelName = 1){
  if (!validateExperiment(exp)) return(NULL)

  plot(model(exp, modelName)$model)
}

#' #' Plot PCA
#' #'
#' #' @description visualization using either factoMineR and factorextra package
#' #' or stats and ggplot2 .. click next to show the next plot
#' #'
#' #' @param data transposed dataframe with m/z as columns
#' #' @param method Which method of calculating the PCA should be used.
#' #' "facto" uses the FactoMineR and factoextra packages. "stats" uses prcomp from
#' #' the basic stats package and ggpubr.
#' #' @param classifiers sample groups as factor
#' #' @export
#' #' @importFrom FactoMineR PCA
#' #' @importFrom factoextra get_pca_var fviz_eig fviz_pca_var fviz_pca_ind
#' #' @importFrom ggpubr ggscatter
#' #' @importFrom stats prcomp
#' #' @importFrom ggplot2 geom_hline geom_vline
#' #' @importFrom dplyr %>%
#' #' @importFrom graphics par
#' plot_PCA <- function(data, method = c("facto", "stats"), classifiers) {
#'   if (method[1] == "facto") {
#'     res_pca <- PCA(data, graph = FALSE)
#'
#'     ## Visualisation of variance explained (plot the variance against the no of dimension)
#'     print(fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 100), main = "PCA - scree plot"))
#'     par(ask = TRUE)
#'
#'     ## Extracting results of variables
#'     var <- get_pca_var(res_pca)
#'     print(fviz_pca_var(res_pca, col.var = "grey", col.circle = "grey", title = "variables-PCA"))
#'
#'     ## Plotting the individuals
#'     par(ask = TRUE)
#'     print(fviz_pca_ind(res_pca,
#'       col.ind = "cos2",
#'       gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#'       repel = TRUE, # Avoid text overlapping (slow if many points)
#'       title = "individuals-PCA - names of the sample"
#'     ))
#'
#'     ## plotting the ellipses
#'     par(ask = TRUE)
#'     fviz_pca_ind(res_pca,
#'       geom.ind = "point", # show points only (nbut not "text")
#'       col.ind = classifiers, # color by groups
#'       palette = "viridis",
#'       addEllipses = TRUE, # Concentration ellipses
#'       legend.title = "Sample type",
#'       title = "PCA samples "
#'     )
#'   } else if (method[1] == "stats") {
#'     ## different way to plot the PCA
#'     ## PCA from basic stats
#'
#'     PCA2 <- prcomp(as.matrix(data), scale. = F) # PCA model using transposed df
#'     PCA_scores <- as.data.frame(PCA2$x) %>% dplyr::select(PC1, PC2)
#'     PCA_scores$Sample <- classifiers ## we add our classifiers here
#'
#'     ## plotting the samples and ellipses using ggplot
#'     print(fviz_eig(PCA2, addlabels = TRUE, ylim = c(0, 100), main = "PCA -scree plot"))
#'     par(ask = TRUE)
#'     ggscatter(PCA_scores,
#'       x = "PC1", y = "PC2",
#'       color = "Sample", shape = "Sample", palette = "aaas",
#'       mean.point = TRUE, ellipse = TRUE, title = "PCA Samples"
#'     ) +
#'       geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
#'       geom_vline(xintercept = 0, linetype = "dashed", color = "black")
#'   }
#' }


#' @title Plot Principle Component Analysis
#' @param se SummarizedExperiment object
#' @param colors optional vector of colors
#' @param assay Assay to use
#' @importFrom ggplot2 ggplot aes scale_shape_manual stat_ellipse
#' @importFrom stats prcomp
#' @export
plotPCA <- function(se, assay = 1, group = metadata(se)$phenotype){
  ratio <- log2(assay(se, assay))
  aliquots <- apply(ratio, 2, function(x) !all(is.na(x)))
  ratio <- ratio[, aliquots]
  se <- se[, aliquots]
  ratio <- as.data.frame(apply(ratio, 2, function(x) {
    x[!is.finite(x)] <- min(x[is.finite(x)]) * 0.01
    x
  }))
  aliquots <- vapply(colnames(ratio), function(x) all(!is.infinite(ratio[, x])), logical(1))
  ratio <- ratio[, aliquots]
  se <- se[, aliquots]

  pc <- prcomp(t(ratio), retx = T,
               scale = T, center = T)

  expvar <- (pc$sdev)^2 / sum(pc$sdev^2)
  pc1name <- sprintf("PC1 (%.1f%%)", expvar[1] * 100)
  pc2name <- sprintf("PC2 (%.1f%%)", expvar[2] * 100)

  df_pc <- data.frame(
    PC1 = pc$x[, 1],
    PC2 = pc$x[, 2],
    Group = as.factor(se[[group]])
  )

  ggplot(df_pc, aes(x = .data$PC1, y = .data$PC2, fill = .data$Group,
                    color = .data$Group)) +
    geom_point(size = 3, stroke = .2) +
    xlab(pc1name) +
    ylab(pc2name) +
    theme_bw() +
    ggtitle("Principle Component Analysis")
}

# df <- data.frame(pvalue = rowData(expModel)$welchTest$p.value,
#            logfc = rowData(expModel)$foldChange[,4])
# df$pvalue <- -log10(df$pvalue)
#
# library(ggplot2)
# ggplot(df, aes(x = logfc, y = pvalue)) + geom_point()


#' @title PCA scree plot
#' @importFrom FactoMineR PCA
#' @importFrom factoextra fviz_eig
#' @param exp
#' @param assay
#' @export
screePCA <- function(exp, assay = 1) {
  if (!validateExperiment(exp)) return(NULL)

  data <- t(assay(exp, assay))
  res_pca <- PCA(data, graph = FALSE)
  fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 100), main = "PCA - scree plot")
}

compoundPCA2 <- function(exp, assay = 1) {
  pca <- prcomp(t(assay(exp, assay)))
  sds <- pca$sdev

}

#' @title PCA variables plot
#' @importFrom FactoMineR PCA
#' @importFrom factoextra fviz_pca_var
#' @param exp
#' @param assay
#' @export
compoundPCA <- function(exp, assay = 1) {
  if (!validateExperiment(exp)) return(NULL)

  data <- t(assay(exp, assay))
  res_pca <- PCA(data, graph = FALSE)
  fviz_pca_var(res_pca, col.var = "grey", col.circle = "grey", title = "variables-PCA")
}


#' @title PCA individuals plot
#' @param exp
#' @param assay
#' @importFrom FactoMineR PCA
#' @importFrom factoextra fviz_pca_ind
#' @export
samplePCA <- function(exp, assay = 1) {
  if (!validateExperiment(exp)) return(NULL)
  data <- t(assay(exp, assay))
  res_pca <- PCA(data, graph = FALSE)
  suppressWarnings(fviz_pca_ind(res_pca,
    col.ind = "cos2",
    gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
    repel = TRUE, # Avoid text overlapping (slow if many points)
    title = "individuals-PCA - names of the sample"
  ))
}

#' @title PCA ellipse plot
#' @importFrom FactoMineR PCA
#' @importFrom factoextra fviz_pca_ind
#' @param data transposed dataframe with m/z as columns
#' @param classifiers sample groups as factor
#' @export
PCA_ellipse <- function(exp, classifiers, assay = 1) {
  if (!validateExperiment(exp)) return(NULL)
  data <- t(assay(exp, assay))
  classifiers <- as.factor(exp[[classifiers]])
  res_pca <- PCA(data, graph = FALSE)
  fviz_pca_ind(res_pca,
    geom.ind = "point", # show points only (nbut not "text")
    col.ind = classifiers, # color by groups
    palette = "viridis",
    addEllipses = TRUE, # Concentration ellipses
    legend.title = "Sample type",
    title = "PCA samples "
  )
}


#' @title PCA ellipse plot from stats
#' @importFrom ggpubr ggscatter
#' @importFrom stats prcomp
#' @importFrom ggplot2 geom_hline geom_vline
#' @importFrom dplyr %>%
#' @param classifiers sample groups as factor
#' @param data transposed dataframe with m/z as columns
#' @export
PCA_ellipse_stats <- function(exp, classifiers, assay = 1) {
  if (!validateExperiment(exp)) return(NULL)
  data <- t(assay(exp, assay))
  classifiers <- as.factor(exp[[classifiers]])
  PCA2 <- prcomp(as.matrix(data), scale. = F) # PCA model using transposed df
  PCA_scores <- as.data.frame(PCA2$x) %>% dplyr::select(PC1, PC2)
  PCA_scores$Sample <- classifiers ## we add our classifiers here
  ggscatter(PCA_scores,
    x = "PC1", y = "PC2",
    color = "Sample", shape = "Sample", palette = "aaas",
    mean.point = TRUE, ellipse = TRUE, title = "PCA Samples"
  ) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black")
}

#' Plot PLS-DA
#'
#' @param data transposed dataframe with m/z as columns
#' @param classifiers samples groups as factor
#' @param comp integer value indicating the component of interest from the object
#' @param method method for contribution "median", "mean"
#' @importFrom mixOmics plsda plotIndiv auroc plotLoadings
#' @importFrom graphics par
plotPLSDA <- function(data, classifiers, comp, method) {
  plsda <- plsda(data, classifiers)
  ## plotting the samples classifiers with ellipses
  plotIndiv(plsda, ind.names = FALSE, star = TRUE, ellipse = TRUE, legend = TRUE, title = "PLS-DA samples")
  par(ask = TRUE)
  ## plotting the ROC curve and calculating the auc of the plsda model
  auroc(plsda, roc.comp = comp, title = "PLS-DA ROC Curve ")
  par(ask = TRUE)
  plotLoadings(plsda, contrib = "max", method = method, comp = comp, title = "PLS-DA contribution")
}


#' @title PLSDA individual plot
#' @importFrom mixOmics plsda plotIndiv
#' @param data transposed dataframe with m/z as columns
#' @param classifiers sample groups as factor
#' @export
PLSDA_ind <- function(exp, classifiers, assay = 1) {
  if (!validateExperiment(exp)) return(NULL)
  data <- t(assay(exp, assay))
  classifiers <- as.factor(exp[[classifiers]])
  plsda <- plsda(data, classifiers)
  ## plotting the samples classifiers with ellipses
  plotIndiv(plsda, ind.names = FALSE, star = TRUE, ellipse = TRUE, legend = TRUE, title = "PLS-DA samples")
}


#' @title PLSDA ROC plot
#' @importFrom mixOmics plsda auroc
#' @param data transposed dataframe with m/z as columns
#' @param classifiers sample groups as factor
#' @param comp integer value indicating the component of interest from the object (default=1)
#' @export
PLSDA_ROC <- function(exp, classifiers, assay = 1, comp = 1) {
  if (!validateExperiment(exp)) return(NULL)
  data <- t(assay(exp, assay))
  classifiers <- as.factor(exp[[classifiers]])
  plsda <- plsda(data, classifiers)
  auroc(plsda, roc.comp = comp, title = "PLS-DA ROC Curve ")
}


#' @title PLSDA ROC plot
#' @importFrom mixOmics plsda plotLoadings
#' @param data transposed dataframe with m/z as columns
#' @param classifiers sample groups as factor
#' @param comp integer value indicating the component of interest from the object (default=1)
#' @param method method for contribution "median" or "mean", default is "median"
#' @export
PLSDA_loadings <- function(exp, classifiers, assay = 1, comp = 1, method = "median") {
  if (!validateExperiment(exp)) return(NULL)
  data <- t(assay(exp, assay))
  classifiers <- as.factor(exp[[classifiers]])
  plsda <- plsda(data, classifiers)
  plotLoadings(plsda, contrib = "max", method = method, comp = comp, title = "PLS-DA contribution")
}


#' heatmap
#' @description heatmap plot of m/z intensities of grouped samples
#' @param data dataframe with m/z as columns
#' @param classifiers samples group as factor
#' @param pretty.order.rows logical vector the default is TRUE
#' @param pretty.order.cols logical vector the default is TRUE
#' @importFrom superheat superheat
#' @importFrom viridis mako
#' @importFrom dplyr %>%
#' @export
heatMap <- function(exp, assay, classifiers, pretty.order.rows = T, pretty.order.cols = T) {
  if (!validateExperiment(exp)) return(NULL)
  data <- t(assay(exp, assay))
  classifiers <- as.factor(exp[[classifiers]])
  superheat <- data %>% superheat(
    left.label.size = 0.05,
    left.label.text.size = 7,
    left.label.col = "white",
    left.label.text.angle = 90, ## the classifiers names angle
    grid.hline.col = "white",
    grid.vline.col = "white",
    grid.hline.size = 5,
    grid.vline.size = 3,
    membership.rows = classifiers,
    heat.pal = viridis::mako(100),
    pretty.order.rows = pretty.order.rows,
    pretty.order.cols = pretty.order.cols,
    scale = F, # scales columns
    padding = 0.1,
    legend.width = 4,
    title = "heatmap of m/z intensity", title.size = 10
  )
}

#' parallel coordinates plot
#'
#' @description parallel coordinates plot with either highlighting a certain value or not
#'
#' @param data dataframe
#' @param tag logical vector TRUE or FALSE if highlighting a certain value is needed or not
#' @param columns the columns index for plotting either intensity or median intensity .. etc
#' @param groupColumn the column name of the groups or the index of the group column
#' @param scale character vector specifying the type of scaling "globalminmax","uniminmax","std","robust"
#' @param boxplot logical vector to add a boxplot to the plot, the default is FALSE
#' @param highlight if tag=TRUE, add highlighting ifelse statement:
#' example : highlight=if_else(data$ppm_difference < 5, "ppm difference less than 5", "ppm difference more than 5")
#' example : highlight=if_else(data$class1 > data$class2 ,"class1","class2")
#' @param highlight_color the highlight color, the default is c("maroon","navy")
#' @importFrom GGally ggparcoord
#' @importFrom ggpubr theme_pubr
#' @importFrom ggplot2 geom_line geom_point scale_color_manual
parallel_coord <- function(data, tag = c(TRUE, FALSE), columns, groupColumn, scale, boxplot = FALSE,
                           highlight = highlight, highlight_color = c("maroon", "navy")) {
  if (tag[1] == FALSE) {
    ggparcoord(data,
      columns = columns, groupColumn = groupColumn,
      showPoints = TRUE,
      boxplot = FALSE,
      scale = scale
    ) +
      labs(x = "Groups", y = "m/z intensity") +
      theme_pubr() +
      theme(
        panel.grid = element_blank(), panel.grid.major.x = element_line(colour = "black"),
        text = element_text(size = 20), legend.position = "right"
      ) +
      geom_line(size = 1.5) +
      geom_point(size = 4) +
      font("xlab", size = 20) +
      font("ylab", size = 20) +
      font("legend.title", face = "bold", size = 15) +
      font("legend.text", size = 15)
  } else if (tag[1] == TRUE) {
    ## tagging a certain value to highlight it in the parallel coords
    tag <- within(data, Highlight <- highlight)
    ggparcoord(tag[order(tag$Highlight), ],
      columns = columns, groupColumn = "Highlight",
      ## we use the extra column we created as groupcolumn to show this highlight
      showPoints = TRUE,
      boxplot = FALSE,
      scale = scale, title = "m/z intensity of Groups - Highlight"
    ) +
      scale_color_manual(values = highlight_color) +
      labs(x = "Groups", y = "m/z intensity") +
      theme(
        plot.title = element_text(size = 20),
        panel.grid = element_blank()
      ) + geom_line(size = 1)
  }
}

# random forest visualization ---------------------------------------------

# random forest visualization ---------------------------------------------

#' Plot random forest cross validation ROC
#' @description plotting ROC curves for different 5 cross-validation sub models vs ROC curve of the final model
#' Note : for this plot cross-validation from sumR package needs to be done first
#' @param roc_all ROC list from randomForest_CV function output
#' @param auc_all AUC list from randomForest_CV function output
#' @param rocfinal rocfinal from RF_model function output
#' @param aucfinal aucfinal from RF_model function output
#' @importFrom graphics par
RF_CV_plots <- function(roc_all, auc_all, rocfinal, aucfinal) {
  # ROC curve of each fold
  par(mfrow = c(2, 3))
  plot(roc_all[[1]], col = "Red", main = paste("Fold_1, AUC:", as.character(round(auc_all[[1]], 3))))
  plot(roc_all[[2]], col = "Red", main = paste("Fold_2, AUC:", as.character(round(auc_all[[2]], 3))))
  plot(roc_all[[3]], col = "Red", main = paste("Fold_3, AUC:", as.character(round(auc_all[[3]], 3))))
  plot(roc_all[[4]], col = "Red", main = paste("Fold_4, AUC:", as.character(round(auc_all[[4]], 3))))
  plot(roc_all[[5]], col = "Red", main = paste("Fold_5, AUC:", as.character(round(auc_all[[5]], 3))))
  plot(rocfinal, col = "blue", main = paste("final model, AUC:", as.character(round(aucfinal, 3))))
}


#' @description plotting ROC curves for any choosen model (final (rocfinal and aucfinal) or
#' from cross validation model(sepcifiy model number roc_all[[1]] and auc_all[[1]]))
#' @param roc ROC of desired model
#' @param auc AUC of desired model
ROC_plot <- function(roc, auc) {
  plot(roc, col = "Red", main = paste("AUC:", as.character(round(auc, 3))))
}

#' random forest plots
#' @description random forest plots .. click next to view the next plots
#' 1. OOB plot 2.important variables plot 3.MDS plot
#' @param model random forest model
#' @param training_set training set that was used for the random forest model with samples column
#' @importFrom randomForest varImpPlot
RandomForestPlots <- function(model, training_set) {
  plot(model)
  par(ask = TRUE)
  # Plot the important variables
  varImpPlot(model, col = "blue", pch = 2)
  par(ask = TRUE)
  MDS_plot(model, training_set)
}



#' random forest plots
#' @title MDS plot
#' @param model random forest model
#' @param training_set training set that was used for the random forest model with samples column
#' @importFrom ggplot2 ggplot geom_label theme_bw ggtitle
#' @importFrom stats cmdscale
#' @importFrom graphics par
MDS_plot <- function(model, training_set) {
  distance_matrix <- as.dist(1 - model$proximity)
  mds <- cmdscale(distance_matrix, eig = TRUE, x.ret = TRUE)
  ## calculate the percentage of variation that each MDS axis accounts for...
  mds_var <- round(mds$eig / sum(mds$eig) * 100, 1)
  ## now make a fancy looking plot that shows the MDS axes and the variation:
  mds_values <- mds$points
  mds_data <- data.frame(
    Sample = rownames(mds_values),
    X = mds_values[, 1],
    Y = mds_values[, 2],
    Status = training_set$samples
  )
  ggplot(data = mds_data, aes(x = X, y = Y, label = Sample)) +
    geom_label(aes(color = Status)) +
    theme_bw() +
    xlab(paste("MDS1 - ", mds_var[1], "%", sep = "")) +
    ylab(paste("MDS2 - ", mds_var[2], "%", sep = "")) +
    ggtitle("MDS plot using (1 - Random Forest Proximities)")
}

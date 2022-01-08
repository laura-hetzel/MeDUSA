#' Volcano plot
#'
#' @param data dataframe result from statistical tests within the sumR package with significant column (true or false)
#' @param xvalues log2 folchange between two groups column
#' @param yvalues nominal p value column
#' @param title title of the plot
#'
#' @return volcano plot
#' @import ggplot2
#' @importFrom ggsci scale_color_aaas
#' @importFrom ggpubr theme_pubr
#' @importFrom ggpubr labs_pubr
volcanoPlot <- function(data, xvalues, yvalues, title){
  ggplot(data) +
    geom_point(aes(x=xvalues, y=-log10(yvalues), colour=significant)) + ## color by significant of fdr <0.1
    ggtitle(title) +
    xlab("log2 fold change") +
    ylab("-log10 nominal p.value") +
    #scale_y_continuous(limits = c(0,50)) +
    theme(legend.position = "right",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25))) +
    scale_color_aaas() +
    theme_pubr() +
    labs_pubr()
}


#' Plot PCA
#'
#' @description visualization using either factoMineR and factorextra package
#' or stats and ggplot2 .. click next to show the next plot
#'
#' @param data transposed dataframe with m/z as columns
#' @param method Which method of calculating the PCA should be used.
#' "facto" uses the FactoMineR and factoextra packages. "stats" uses prcomp from
#' the basic stats package and ggpubr.
#' @param classifiers sample groups as factor
#' @return
#' @export
#' @importFrom FactoMineR PCA
#' @importFrom factoextra get_pca_var
#' @importFrom factoextra fviz_eig
#' @importFrom factoextra fviz_pca_var
#' @importFrom factoextra fviz_pca_ind
#' @importFrom ggpubr ggscatter
#' @importFrom stats prcomp
#' @importFrom ggplot2 geom_hline
#' @importFrom ggplot2 geom_vline
#'
#' @examples
plot_PCA <- function(data, method = c("facto", "stats"),classifiers){
  if(method[1] == "facto"){
    res_pca <- PCA(data, graph = FALSE)

    ## Visualisation of variance explained (plot the variance against the no of dimension)
    print(fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 100),main="PCA - scree plot" ))
    par(ask=TRUE)

    ## Extracting results of variables
    var <- get_pca_var(res_pca)
   print(fviz_pca_var(res_pca,col.var = "grey",col.circle = "grey",title="variables-PCA"))

    ##Plotting the individuals
    par(ask=TRUE)
    print(fviz_pca_ind(res_pca, col.ind = "cos2",
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE,# Avoid text overlapping (slow if many points)
                 title="individuals-PCA - names of the sample"))

    ## plotting the ellipses
    par(ask=TRUE)
    fviz_pca_ind(res_pca,
                 geom.ind = "point", # show points only (nbut not "text")
                 col.ind = classifiers, # color by groups
                 palette = "viridis",
                 addEllipses = TRUE, # Concentration ellipses
                 legend.title = "Sample type",
                 title = "PCA samples "
    )
  }else if(method[1] == "stats"){
    ## different way to plot the PCA
    ## PCA from basic stats

    PCA2 <- prcomp(as.matrix(data), scale. = F) #PCA model using transposed df
    PCA_scores <- as.data.frame(PCA2$x) %>% dplyr::select(PC1, PC2)
    PCA_scores$Sample <- classifiers ## we add our classifiers here

    ##plotting the samples and ellipses using ggplot
    print(fviz_eig(PCA2, addlabels = TRUE, ylim = c(0, 100),main="PCA -scree plot" ))
    par(ask=TRUE)
    ggscatter(PCA_scores, x = "PC1", y = "PC2",
              color = "Sample", shape = "Sample", palette = "aaas",
              mean.point = TRUE, ellipse = TRUE, title = "PCA Samples") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black")

  }

}
#' Plot PLS-DA
#'
#' @param data transposed dataframe with m/z as columns
#' @param classifiers samples groups as factor
#' @param comp integer value indicating the component of interest from the object
#' @param method method for contribution "median", "mean"
#' @importFrom mixOmics plsda
#' @importFrom mixOmics plotIndiv
#' @importFrom mixOmics auroc
#' @importFrom mixOmics plotLoadings
#'
#' @examples

plotPLSDA <- function(data, classifiers,comp,method){
  plsda <-  plsda(data, classifiers)

  ## plotting the samples classifiers with ellipses
  plotIndiv(plsda, ind.names = FALSE, star = TRUE, ellipse = TRUE, legend = TRUE,title = "PLS-DA samples")
   par(ask=TRUE)

  ## plotting the ROC curve and calculating the auc of the plsda model
   auroc(plsda,roc.comp = comp, title = "PLS-DA ROC Curve ")
 par(ask=TRUE)

   plotLoadings(plsda, contrib = 'max', method = method, comp = comp,title="PLS-DA contribution")

}

#' heatmap
#' @description heatmap plot of m/z intensities of grouped samples
#' @param data dataframe with m/z as columns
#' @param classifiers samples group as factor
#' @param pretty.order.rows logical vector the default is TRUE
#' @param pretty.order.cols logical vector the default is TRUE
#' @importFrom superheat superheat
#' @importFrom viridis mako
heatMap<-function(data,classifiers,pretty.order.rows=T,pretty.order.cols =T){
 superheat <- data %>% superheat(left.label.size = 0.05,
            left.label.text.size =7,
            left.label.col = "white",
            left.label.text.angle = 90, ## the classifiers names angle
            grid.hline.col = "white",
            grid.vline.col = "white",
            grid.hline.size = 5,
            grid.vline.size = 3,
            membership.rows = classifiers,
            heat.pal = viridis::mako(100),
            pretty.order.rows = pretty.order.rows,
            pretty.order.cols =pretty.order.cols,
            scale = F, #scales columns
            padding = 0.1,
            legend.width=4,
            title="heatmap of m/z intensity",title.size=10)
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
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_color_manual
parallel_coord <- function(data, tag = c(TRUE,FALSE),columns,groupColumn,scale,boxplot = FALSE,
                           highlight=highlight,highlight_color=c("maroon","navy")){
  if(tag[1] == FALSE){

    ggparcoord(data,
               columns = columns, groupColumn = groupColumn,
               showPoints = TRUE,
               boxplot = FALSE,
               scale=scale) +
      labs(x = "Groups",y= "m/z intensity") +
      theme_pubr()+
      theme(panel.grid = element_blank(),panel.grid.major.x=element_line(colour="black"),
            text = element_text(size=20),  legend.position="right")+
      geom_line(size=1.5)+
      geom_point(size=4)+
      font("xlab", size = 20)+
      font("ylab", size = 20)+
      font("legend.title", face = "bold",size=15)+
      font("legend.text",size=15)

  }
  else if(tag[1] == TRUE){

    ## tagging a certain value to highlight it in the parallel coords
    tag<-within(data, Highlight<-highlight)

    ggparcoord(tag[order(tag$Highlight),], columns = columns, groupColumn = "Highlight",
               ## we use the extra column we created as groupcolumn to show this highlight
               showPoints = TRUE,
               boxplot = FALSE,
               scale=scale,title="m/z intensity of Groups - Highlight" ) +
      scale_color_manual(values=highlight_color)+
      labs(x = "Groups",y= "m/z intensity") +
      theme(plot.title = element_text(size=20),
            panel.grid = element_blank())+geom_line(size=1)


  }

  }

# random forest visualization ---------------------------------------------

#' Plot random forest cross validation ROC
#' @description plotting ROC curves for different 5 cross-validation sub models vs ROC curve of the final model
#' Note : for this plot cross-validation from sumR package needs to be done first
#' @param model the final model
#' @param test_set the final test set with samples column
#' @importFrom pROC roc
#' @importFrom pROC auc
RF_CV_plots<-function(model,test_set){
  test_set$samples <- as.factor(test_set$samples)
  results <- predict(model,newdata=test_set,type="class")
  results<-as.data.frame(cbind("actual"=test_set$samples,"prediction"=results))
  rocfinal<- roc(results$actual,results$prediction )
  aucfinal<- auc(rocfinal)

  # ROC curve of each fold
  par(mfrow=c(2,3))
  plot(roc_all[[1]],col = "Red", main = paste("Fold_1, AUC:", as.character(round(auc_all[[1]], 3))))
  plot(roc_all[[2]],col = "Red", main = paste("Fold_2, AUC:", as.character(round(auc_all[[2]], 3))))
  plot(roc_all[[3]],col = "Red", main = paste("Fold_3, AUC:", as.character(round(auc_all[[3]], 3))))
  plot(roc_all[[4]],col = "Red", main = paste("Fold_4, AUC:", as.character(round(auc_all[[4]], 3))))
  plot(roc_all[[5]],col = "Red", main = paste("Fold_5, AUC:", as.character(round(auc_all[[5]], 3))))
  plot(rocfinal,col = "blue", main = paste("final model, AUC:", as.character(round(aucfinal, 3))))

}


#' random forest plots
#' @description random forest plots .. click next to view the next plots
#' 1. OOB plot 2.important variables plot 3.MDS plot
#' @param model random forest model
#' @param training_set training set that was used for the random forest model with samples column
#' @importFrom randomForest varImpPlot
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_label
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 ggtitle
RandomForestPlots<-function(model,training_set){

  plot(model)
  par(ask=TRUE)

  # Plot the important variables
  varImpPlot(model,col="blue",pch= 2)
  par(ask=TRUE)


  # MDS plot
  distance_matrix <- as.dist(1-model$proximity)

  mds <- cmdscale(distance_matrix, eig=TRUE, x.ret=TRUE)

  ## calculate the percentage of variation that each MDS axis accounts for...
  mds_var <- round(mds$eig/sum(mds$eig)*100, 1)

  ## now make a fancy looking plot that shows the MDS axes and the variation:
  mds_values <- mds$points
  mds_data <- data.frame(Sample=rownames(mds_values),
                         X=mds_values[,1],
                         Y=mds_values[,2],
                         Status=training_set$samples)

  ggplot(data=mds_data, aes(x=X, y=Y, label=Sample)) +
    geom_label(aes(color=Status))+
    theme_bw() +
    xlab(paste("MDS1 - ", mds_var[1], "%", sep="")) +
    ylab(paste("MDS2 - ", mds_var[2], "%", sep="")) +
    ggtitle("MDS plot using (1 - Random Forest Proximities)")
}

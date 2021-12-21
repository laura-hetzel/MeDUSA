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
plotPCA <- function(data, method = c("facto", "stats"),classifiers){
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


#' Volcano plot
#'
#' @param data 
#' @param xvalues 
#' @param yvalues 
#' @param title 
#'
#' @return
#' @export
#' @import ggplot2
#'
#' @examples
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
#' or stats and ggplot2
#'
#' @param master_df_pos_t transposed dataframe
#' @param method Which method of calculating the PCA should be used. 
#' "facto" uses the FactoMineR and factoextra packages. "stats" uses prcomp from
#' the basic stats package and ggpubr.
#'
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
plotPCA <- function(master_df_pos_t, method = c("facto", "stats")){
  if(method[1] == "facto"){
    res_pca <- PCA(master_df_pos_t, graph = FALSE)
    
    ## Visualisation of variance explained (plot the variance against the no of dimension)
    fviz_eig(res_pca, addlabels = TRUE, ylim = c(0, 50),main="PCA - scree plot" )
    
    
    ## Extracting results of variables
    var <- get_pca_var(res_pca)
    fviz_pca_var(res_pca,col.var = "grey",col.circle = "grey",title="variables-PCA")
    
    ##Plotting the individuals
    
    fviz_pca_ind(res_pca, col.ind = "cos2",
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE,# Avoid text overlapping (slow if many points)
                 title="individuals-PCA - names of the sample"
    )
    
    ## plotting the ellipses
    
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
    
    PCA2 <- prcomp(as.matrix(master_df_pos_t), scale. = F) #PCA model using transposed df
    fviz_eig(PCA2, addlabels = TRUE, ylim = c(0, 100),main="PCA -scree plot" )
    PCA_scores <- as.data.frame(PCA2$x) %>% dplyr::select(PC1, PC2)
    PCA_scores$Sample <- classifiers ## we add our classifiers here
    
    ##plotting the samples and ellipses using ggplot
    ggscatter(PCA_scores, x = "PC1", y = "PC2",
              color = "Sample", shape = "Sample", palette = "aaas",
              mean.point = TRUE, ellipse = TRUE, title = "PCA ", subtitle = "PC1(25.4%) - PC2(11.2%)") + ## add PC 1 and PC 2 percentage from scree plot
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "black")
    
  }
  
}

#' Plot PLS-DA
#'
#' @param master_df_pos_t 
#' @param classifiers 
#'
#' @return
#' @export
#' @importFrom mixOmics plsda
#' @importFrom mixOmics plotIndiv
#' @importFrom mixOmics auroc
#' @importFrom mixOmics plotLoadings
#'
#' @examples
plotPLSDA <- function(master_df_pos_t, classifiers){
  plsda <-  plsda(master_df_pos_t, classifiers)
  
  ## plotting the samples classifiers with ellipses
  plotIndiv(plsda, ind.names = FALSE, star = TRUE, ellipse = TRUE, legend = TRUE,title = "PLS-DA samples")
  
  ## plotting the ROC curve and calculating the auc of the plsda model
  auc.plsda <- auroc(plsda,roc.comp = 1, title = "PLS-DA ROC Curve ")
  
  ## plotting the pls-da contribution to median
  plotLoadings(plsda, contrib = 'max', method = 'median', comp = 1,title="pls-da contribution by median")
  
}
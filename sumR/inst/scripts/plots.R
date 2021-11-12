# Visualization -----------------------------------------------------------

## heatmap
## we need transposed df
superheat <- master_df_t %>%
  superheat(left.label.size = 0.05,
            left.label.text.size = 10,
            left.label.col = "white",
            left.label.text.angle = 90, ## the classifiers names angle
            grid.hline.col = "white",
            grid.vline.col = "white",
            grid.hline.size = 5,
            grid.vline.size = 3,
            membership.rows = classifiers,
            heat.pal = viridis::mako(100),
            pretty.order.rows = T,
            pretty.order.cols =T,
            scale = F, #scales columns
            padding = 0.1,
            legend.width=4)


## but we plot the points of log2 folchange between two groups and the values of nominal p value but it will be colored
# by significance according to adjusted p value
## we can use here the results of wilcox or welch t test
volcanoPlot(wilcox,wilcox$log2fc_M1_M2,wilcox$p.value, "Volcanoplo P-values M1 (welch t test and mann whiteny u)
            significance of adjusted p values with FDR < 0.1")


# Random forest visualization  ------------------------------------------------------------------
##AUC and roc
## if we did a random forest cross validation using 5 models to obtain the final model

## and then checking the final model against a permuted model we can plot the ROC Curve of all these model to compare them against
## each other , ideally the final model will have the highest accuracy and the permuted model will have the lowest accuracy

roc_1<- roc(results_1$actual,results_1$prediction ) ## the actual vs the predictions of the first model .. etc
auc_1<- auc(roc_1)

roc_2<- roc(results_2$actual,results_2$prediction )
auc_2<- auc(roc_2)

roc_3<- roc(results_3$actual,results_3$prediction )
auc_3<- auc(roc_3)

roc_4<- roc(results_4$actual,results_4$prediction )
auc_4<- auc(roc_4)

roc_5<- roc(results_5$actual,results_5$prediction )
auc_5<- auc(roc_5)

roc_final<- roc(results$actual,results$prediction )
auc_final<- auc(roc_final)

roc_permuted<- roc(results_permuted$actual,results_permuted$prediction )
auc_permuted<- auc(roc_permuted)


# ROC curve of each fold
par(mfrow=c(2,4))
plot(roc_1,col = "Red", main = paste("Fold_1, AUC:", as.character(round(auc_1, 3))))
plot(roc_2,col = "Red", main = paste("Fold_2, AUC:", as.character(round(auc_2, 3))))
plot(roc_3,col = "Red", main = paste("Fold_3, AUC:", as.character(round(auc_3, 3))))
plot(roc_4,col = "Red", main = paste("Fold_4, AUC:", as.character(round(auc_4, 3))))
plot(roc_5,col = "Red", main = paste("Fold_5, AUC:", as.character(round(auc_5, 3))))
plot(roc_final,col = "blue", main = paste("final model, AUC:", as.character(round(auc_final, 3))))
plot(roc_permuted,col = "green", main = paste("permuted model, AUC:", as.character(round(auc_permuted, 3))))


# Plot the important variables of a random forest model , importance = True in the model to plot this
varImpPlot(model,col="blue",pch= 2)

# model evaluation of random forest
## we still have to know more about this package library(rfUtilities)
## it is used from random forest cross validation



rf_cv <- rf.crossValidation(model, training_set,seed=555, p=0.1, n=99, ntree=500)
## use the same seed,training set , n_tree as the final model , n= number of cross validation the default is 99
#p ratio of splitting the training set for cross validation aka Proportion data withhold (default p=0.10)

# Plot cross validation versus model producers accuracy
par(mfrow=c(1,2))
plot(rf_cv, type = "cv", main = "CV producers accuracy")
plot(rf_cv, type = "model", main = "Model producers accuracy")

# Plot cross validation versus model oob
par(mfrow=c(1,2))
plot(rf_cv, type = "cv", stat = "oob", main = "CV oob error")
plot(rf_cv, type = "model", stat = "oob", main = "Model oob error")


## decision boundary from library(ElemStatLearn)
# you can do this for training set or test set ..

set <- test_set  # or training set

## i made a model specifically for this visualization because I can't find a way around it .
# this plot uses two features only for the model or we are trying to search more about that because Ahmed wants to :D
model_vis<-randomForest(x=cbind(training_set[,2:3]),y=training_set$samples,data=training_set,ntree=500)


X1 = seq(min(set[,2]) - 1, max(set[,2]) + 1, by = 0.01)
X2 = seq(min(set[,3]) - 1, max(set[,3]) + 1, by = 0.01)
grid_set = expand.grid(X1, X2) ## we make a grid using the sequence from min and max values of the two features we used in the model
colnames(grid_set) = colnames(set[,2:3]) ## setting the same name of the feature to our grid set
y_grid = predict(model_vis, newdata = grid_set, type = 'class') ## predict our y_grid which we will use it to color
## the background for this plot creating a decision boundary of predicition between two groups
## using the grid set we created and our model (2 features only !)
# NOTE we need class here because we have a y_grid is a matrix!
plot(set[,2:3],
     main = 'Random Forest classification (test set)',
     xlab = 'intensity', ylab = 'intensity',
     xlim = range(X1), ylim = range(X2)) # this bit creates the limits to the values plotted this is also a part of the MAGIC as it creates the line between green and red
contour(X1, X2, matrix(as.numeric(y_grid), length(X1), length(X2)), add = TRUE)
# here we run through all the y_pred data and use if else to color the dots
# note the dots are the real data, the background is the pixel by pixel determination of prediction
# graph the dots on top of the background give you the image
points(grid_set, pch = '.', col = ifelse(y_grid == "GroupA", 'springgreen3', 'tomato')) ## plotinng background
points(set[,2:3], pch = 21, bg = ifelse(set[, 1] == "GroupA", 'green4', 'red3')) ## plotting the real data of either training or test set


# parallel coordinates ----------------------------------------------------



## plotting the parallel coordinates from GGally package
## we use a dataset which have the group column = lipids ,, median scaled transformed intensity of each group or class as columns
ggparcoord(median_data,
           columns = 2:5, groupColumn = "Lipids",  ## columns = here we specify the range of columns of median data for each group
           showPoints = TRUE,                       ##groupColumn = either significant peaks (m/z) or lipids name
           #boxplot = TRUE,                     ## we can add boxplot to the plot
           scale="globalminmax") +              ## scale globalminmax plot the data as it is without scaling on the plot
  ## there is uniminmax scaling the data between 0-1
  labs(x = "groups",y= "Log2 pareto scaled intensity") +
  theme_pubr()+
  theme(panel.grid = element_blank(),panel.grid.major.x=element_line(colour="black"),
        text = element_text(size=20),  legend.position="right")+
  geom_line(size=1.5)+
  geom_point(size=4)+
  font("xlab", size = 20)+
  font("ylab", size = 20)+
  font("legend.title", face = "bold",size=15)+
  font("legend.text",size=15)


## if we want to add p values to the plot and we have more than two groups so we can specify our comparison group
my_comparisons <- list( c("A", "B"), c("A", "C"), c("A", "D"),
                        c("B", "C"),c("B", "D"),c("C", "D"))
plot+stat_compare_means(comparisons = my_comparisons)+  ## we add this line to the plot
  stat_compare_means(label.y = 10)  ## computing global p value , label.y= the position of p value printed in the plot



## tagging a certain value to highlight it in the parallel coords
## here we were tagging a slope if it is positive or negative using a slope column we already had in our data
## you can tag anything in your data set
tag<-within(median_data, highlight<-if_else(slope>0, "positive", "negative"))

ggparcoord(slope_tag[order(tag$highlight),], columns = 2:5, groupColumn = "highlight",
           ## we use the extra column we created as groupcolumn to show this highlight
           showPoints = TRUE,
           #boxplot = TRUE,
           scale="globalminmax",title="Median of logged and scaled data(pareto_scaling) - Slope highlight" ) +
  scale_color_manual(values=c("red","navy"))+
  labs(x = "groups",y= "Log2 pareto scaled intensity") +
  theme(plot.title = element_text(size=20),
        panel.grid = element_blank())+geom_line(size=1)

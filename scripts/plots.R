#' @title Plot feature coverage when adding cells
#' @param cells
#' @param assay
#' @param by
#' @param seed
#' @importFrom stats loess
#' @export
featureCovPlot <- function(cells, assay = 1, by = 5, seed = 42){
  set.seed(seed)
  cellNames <- sample(colnames(cells))
  is <- seq(by, ncol(cells), by)
  feats <- vapply(is, function(i){
    m <- assay(cells[, cellNames[1:i]], assay)
    nrow(m) - sum(rowSums(is.na(m)) == ncol(m))
  }, double(1))
  df <- data.frame(Cells = is, Features = feats)
  loess_fit <- loess(Features ~ Cells, df)
  df2 <- as.data.frame(predict(loess_fit))
  df2$Cells <- df$Cells
  colnames(df2) <- c("Features", "Cells")

  ggplot(df, aes(x = Cells, y = Features)) +
    geom_line(aes(colour = "Features")) +
    geom_line(data = df2, aes(colour = "Trend")) +
    labs(colour = "Line Type")
}

#' @title Plot features per cell
#' @param cells
#' @param assay
#' @export
featureCellPlot <- function(cells, assay = "Area"){
  df <- as.data.frame(stack(colSums(!is.na(assay(cells, assay)))))
  df <- df[order(df$values, decreasing = T), ]
  df$Cell <- 1:nrow(df)
  colnames(df)[1] <- "Features"
  ggplot(df, aes(y = Features, x = Cell)) + geom_bar(stat = "identity")
}

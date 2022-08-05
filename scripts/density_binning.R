library(sumR)

# Peak Picking (centroiding)
files <- list.files(file.path(r"(F:\Myrthe\mzml)"), full.names = T)
fileList <- prepareFiles(files[1:2], massWindow = c(100, 200), polarity = "+", centroid = T)

file = 1
scan = 1
peaks <- fileList





grid <- expand.grid(method = c("binning"),# "clustering", "density"),
            ppm = rev(c(5, 10, 15, 20)))
grid
for (i in 1:nrow(grid)) {
  method = as.vector(grid$method[i])
  ppm = as.vector(grid$ppm[i])
  cat(method, ppm, "\n")
  groups <- spectraAlignment(fileList, method = method, cores = 1, ppm = ppm)
  cells <- cellAlignment(groups[groups$npeaks >= 5, ], method = method, ppm = ppm)
  saveRDS(cells, file = sprintf("SE_method_%s__ppm_%s__npeaks_%s.RDS", method, ppm, 5))


  cells <- cellAlignment(groups[groups$npeaks >= 10, ], method = method, ppm = ppm)
  saveRDS(cells, file = sprintf("SE_method_%s__ppm_%s__npeaks_%s.RDS", method, ppm, 10))

  cells <- cellAlignment(groups[groups$npeaks >= 20, ], method = method, ppm = ppm)
  saveRDS(cells, file = sprintf("SE_method_%s__ppm_%s__npeaks_%s.RDS", method, ppm, 20))
}


barplot_counts <- function(minSamples = 1){
  files <- sort(list.files(pattern = ".*binning.*RDS"))

  exps <- lapply(files, readRDS)
  exps <- setNames(exps, basename(files))
  exps <- lapply(exps, function(exp){
    exp[rowSums(!is.na(assay(exp))) >= minSamples, ]
  })

  df <- data.frame(comps = unlist(lapply(exps, nrow)))
  df$ppm <- c(rep(10, 4), rep(15, 4), rep(20, 4), rep(5, 4))
  df$nscans <- as.factor(rep(c(1, 10, 20, 5), 4))
  library(ggplot2)
  ggplot(df, aes(y = comps, x = ppm, fill = nscans)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_y_continuous(labels = scales::comma) +
    ylab("Compounds") +
    ggtitle("Compounds per binning ppm",
            "With minimum occurences (scans)")
}

p <- cowplot::plot_grid(
  plotlist = list(
    barplot_counts(),
    barplot_counts(5),
    barplot_counts(10),
    barplot_counts(20)
  ), labels = c("1", "5", "10", "20")
)

library(SummarizedExperiment)

files <- sort(list.files(pattern = ".*binning.*RDS"))
exps <- lapply(files, readRDS)
exps <- setNames(exps, basename(files))

names(exps)
exp <- exps[[12]]
exp1 <- exp[rowSums(!is.na(assay(exp))) >= 6, ]
rowData(exp[rowData(exp)$mzmed > 372.2 & rowData(exp)$mzmed < 372.4, ])
assay(exp[17147, ])

exp <- exps[[12]]

exp2 <- exp[, grep("BLANK", colnames(exp), invert = T)]
exp2 <- exp2[rowSums(!is.na(assay(exp2))) >= 6, ]
rowData(exp2[rowData(exp2)$mzmin > 372.1 & rowData(exp2)$mzmax < 372.4, ])


rowData(exp1)
rowData(exp2)

df <- do.call(rbind, lapply(1:20, function(minSamples){
  exps <- lapply(exps, function(exp){
    exp <- exp[, grep("BLANK", colnames(exp), invert = T)]
    exp[rowSums(!is.na(assay(exp))) >= minSamples, ]
  })
  data.frame(comps = unlist(lapply(exps, nrow)), minSamples = minSamples)
}))


df

df$ppm <- factor(rep(c(rep(10, 4), rep(15, 4), rep(20, 4), rep(5, 4)), 20))
df$nscans <- as.factor(rep(rep(c(1, 10, 20, 5), 4), 20))
library(ggplot2)

p2 <- ggplot(df, aes(y = comps, x = minSamples, group = ppm, color = ppm)) +
  geom_line(lwd = 1, lty = 1) +
#  geom_bar(stat = "identity", position = "dodge") +
#  scale_y_continuous(labels = scales::comma) +
  scale_y_log10() +
  ylab("Compounds") +
  ggtitle("minimum occurrences scans versus samples") +
  facet_wrap(~ nscans, scales = "free") +
  theme_minimal()

p1 <- ggplot(df, aes(y = comps, x = minSamples, group = ppm, color = ppm)) +
  geom_line(lwd = 1, lty = 1) +
  #  geom_bar(stat = "identity", position = "dodge") +
  scale_y_continuous(labels = scales::comma) +
  #scale_y_log10() +
  ylab("Compounds") +
  ggtitle("minimum occurrences scans versus samples") +
  facet_wrap(~ nscans, scales = "free") +
  theme_minimal()

p

p1

p2



#fileList <- peakFilter(fileList, peakValue = 5000)
#lapply(fileList[[2]], nrow)[1:5]

groups <- spectraGrouping(fileList, method = "clustering", cores = 8, ppm = 30)
cells <- sumR:::cellClustering(groups, ppm = 30)


groups[groups$mz > 372.22 & groups$mz < 372.25, ]
groups <- groups[groups$npeaks >= 5, ]

cells <- sumR:::cellBinning(groups, tolerance = 0.2)

rowData(cells[rowData(cells)$mzmed > 372.2 & rowData(cells)$mzmed < 372.3, ])

cells <- sumR:::cellClustering(groups, ppm = 30)

cells


a <- sumR:::spectraClustering(fileList[1:3])
b <- sumR:::spectraDensity(fileList)
c <- sumR:::spectraBinning(fileList)



which.min(a$mz)
c2 <- sumR:::cellBinning(a, tolerance = 0.002)

sum(is.na(assay(c2))) / length(c(t(assay(c2))))

c3 <- sumR:::cellDensity(a)




c[c2$peakidx[[3]], ]

spectra <- c
n <- length(unique(spectra$sample))
groups <- rep(1, n)
res <- xcms::do_groupChromPeaks_density(as.data.frame(spectra), groups,
                                        binSize = 0.05, minFraction = 0)
res <- res[, grep("mz|npeaks|peakidx", colnames(res))]

samples <- as.character(unique(spectra$sample))

m <- matrix(NA, nrow = length(bins), ncol = length(samples),
            dimnames = list(1:length(bins), 1:length(samples)))

noise <- matrix(NA, nrow = length(bins), ncol = length(samples),
                dimnames = list(1:length(bins), 1:length(samples)))

print(head(spectra))

print(head(binned))

print(m[1:5, ])

pbapply::pblapply(1:nrow(binned), function(i) {
  idx <- as.integer(binned$peakidx[[i]])

  cells <- spectra$sample[idx]

  noise[i, cells] <<- spectra$noise[idx]
  m[i, cells] <<- spectra$i[idx]
})

print(m[1:5, ])

se <- SummarizedExperiment(
  assays = list(Area = m, Noise = noise),
  rowData = DataFrame(binned)
)
if (!is.null(cellData)) colData(se) <- DataFrame(cellData[colnames(m), ])
se


df_list <- sumR:::doBinning(as.data.frame(c), split = "sample", tolerance = 5e-6)

df <- do.call(rbind, df_list)
rownames(df) <- rownames(cells)

head(df)
head(cells[,1:3])

df$mzdiff <- df$mz - cells$mz
df
df <- df[order(df$mz), ]
bins <- split.data.frame(df, df$mz)

res <- DataFrame(
  mzmed = unique(df$mz),
  mzmin = sapply(bins, function(x) min(x$mz - x$mzdiff)),
  mzmax = sapply(bins, function(x) max(x$mz - x$mzdiff)),
  npeaks = sapply(bins, nrow)
)
res

res$peakidx <- lapply(bins, function(x) as.integer(rownames(x)))
res$i <- vapply(res$peakidx, function(idx) sum(cells$i[idx]), double(1))
res$noise <- vapply(res$peakidx, function(idx) sum(cells$noise[idx]), double(1))

head(res)

cells[res$peakidx[[3]], ]


targets <- c(675.483905, 702.567575, 728.583225, 730.598875,
             757.562155, 785.593455, 814.692776, 842.724076,
             843.671706)

res <- lapply(targets, function(target){
  lapply(fileList, function(spectra){
    sum(sapply(spectra, function(spectrum){
      s1 <- spectrum[,1] + 0.003
      s2 <- spectrum[,1] - 0.003
      any(s1 > target & s2 < target)
    }))
  })
})

df <- as.data.frame(do.call(cbind, res))
colnames(df) <- as.character(round(targets, 3))
df

target <- targets[1]
spectra <- fileList[[2]]



base <- 1
res <- pbapply::pblapply(fileList, function(l) {
  non_nulls <- !vapply(l, is.null, logical(1))
  l <- l[non_nulls]
  if (length(l) == 0) return(NULL)
  df <- data.frame(sample = rep(which(non_nulls), sapply(l, nrow)), do.call(rbind, l)) #
  if (nrow(df) == 0)
    return(NULL)
  rownames(df) <- base:(base + nrow(df) - 1)
  base <<- base + nrow(df)
  as.data.frame(df)
})

res2 <- do.call(rbind, res[[1]])

peaks <- res[[1]]
peaks$rt <- -1
sampleGroups <- unique(peaks$sample)
groups <- xcms::do_groupPeaks_mzClust(peaks, minFraction = 0, sampleGroups = rep(1, length(unique(peaks$sample))))

df <- DataFrame(groups$featureDefinitions)
df$peakidx <- groups$peakIndex

df

peaks[df[df$mzmed > 675.4 & df$mzmed < 675.6, ]$peakidx[[8]], ]

head(res[1])

a <- spectraBinning(res[2:length(res)], tolerance = 5e-6)

df_list <- sumR:::doBinning(a, split = NULL,
                     tolerance = 5e-6)

b <- do.call(rbind, lapply(1:length(df_list), function(i){
  a[[i]]$mzAlign <- df_list[[i]]$mz
  a[[i]]$sample <- i
  a[[i]]
}))
head(b[order(b$mzAlign), ], 10)


max(a[[1]]$mz - df_list[[1]]$mz)
head(a[[1]])
head(df_list[[1]])


a <- do.call(rbind, a)
rsd <- function(x) sd(x) / mean(x) * 100

b <- vapply(seq_len(nrow(a)), function(x){
  rsd(res2[a$peakidx[[x]], ]$i)
}, double(1))

a$peakidx[[4001]]



rsd(res2[a$peakidx[[99]], ]$mz)

mean(x)
sd(x)

res2[a[[2]]$peakidx[[1]], ]

res2[as.integer(unlist(a[[2]]$peakidx[2])), ]



res[[1]][a[[1]]$peakidx[[5]], ]


a[[1]]

res[[1]][ unlist(a[[1]]$peakidx[12]), ]


a[[1]]
names(a) <- names(fileList)


df <- do.call(rbind, lapply(1:length(a), function(i) {
  df <- a[[i]]

  if (is.null(df)) return(NULL)
  if (nrow(df) == 0) return(NULL)
  df$cell <- names(a)[i]
  df
}))
df

a

peakidx <- df$peakidx
df_list <- sumR:::doBinning(as.data.frame(df[, c("mz", "i", "cell")]), split = "cell",
                            tolerance = 5e-6)

df2 <- do.call(rbind, df_list)
df2 <- df2[order(df2$mz), ]

rownames(df2) <- rownames(df)



#df2 <- df2[order(df2$mz), ]
#df2$oldmz <- df$mz[order(df$mz)]
#df2$mzdiff <- df2$mz - df2$oldmz
#df2 <- df2[order(rownames(df2)), ]

df2 <- data.frame(df2, cell = rep(names(df_list), sapply(df_list, nrow)))

bins <- split.data.frame(df2, df2$mz)

length(bins)

res2[unlist(df[as.integer(rownames(bins[[13]])), ]$peakidx), ]



df[rownames(bins[[11]]), ]

df[41995, ]

rownames(df) <- as.character(rownames(df))


bins[[11]]
df
df[as.integer(rownames(bins[[11]])), ]


res[[1]][434381, ]

a[1]
unlist(a[[1]]$peakidx[1])

peakidx[as.integer(rownames(bins[[11]]))]

nrow(df)

res2[434381, ]

res2[unlist(peakidx[as.integer(rownames(bins[[11]]))]), ]


unlist(peakidx[as.integer(rownames(bins[[3]]))])

peakidx[[]]

tail(subdf)

res2[as.integer(unlist(peakidx[as.integer(rownames(bins[[11]]))])), ]

spectra <- res[[1]]
sumR:::doBinning(list(subdf), tolerance = 5e-6)



spectra[unlist(df[3, "peakidx"]), ]

















a <- spectraBinning(res, tolerance = 5e-6)
res[[1]]

a[[1]][5, ]
res[[1]][unlist(a[[1]][5, ]$peakidx), ]

res[[1]][ unlist(a[[1]][5, ]$peakidx), ]



spectra <- res[[1]]
df_list <- sumR:::doBinning(spectra, split = "scan", tolerance = 5e-6)

df <- do.call(rbind, df_list)

df <- df[order(df$mz), ]
df$oldmz <- spectra$mz[order(spectra$mz)]

df$mzdiff <- df$mz - df$oldmz
rownames(df) <- 1:nrow(df)
#df <- df[order(rownames(df)), ]
bins <- split.data.frame(df, df$mz)

head(bins[[5]])
head(df)

head(res[[1]])

res <- res[!vapply(res, is.null, logical(1))]
res <- do.call(rbind, res)
table(res$sample)

bins <- sumR:::doBinning(res, split = "sample", tolerance = 5e-6)
bins[[2]]

a <- split.data.frame(res, res$sample)
head(a[[1]])
head(bins[[1]], 30)

df <- do.call(rbind, bins)
df <- df[order(df$mz), ]
head(df, 200)

split.data.frame(df, df$mz)[[1]]


length(bins)
head(bins[[1]])
head(res[[10]])


if (typeof(res[[1]]) != "list") spectra <- split.data.frame(res[[1]], res[[1]][, "sample"])

class(res[1])


df <- do.call(rbind, res2)

df <- df[order(df$mz), ]
df$oldmz <- res$mz[order(res$mz)]

df$mzdiff <- df$mz - df$oldmz
df <- df[order(rownames(df)), ]
bins <- split.data.frame(df, df$mz)
df <- data.frame(
  mz = unique(df$mz), npeaks = sapply(bins, nrow),
  i = sapply(bins, function(x) switch(method, "sum" = sum(x$i), "mean" = mean(x$i))),
  SNR = sapply(bins, function(x) mean(x$snr)),
  mzmin = sapply(bins, function(x) min(x$mzdiff)),
  mzmax = sapply(bins, function(x) max(x$mzdiff))
)

head(res2[[1]])
head(res[[1]])




res <- sumR:::doBinning(res)

nonNull <- which(!vapply(res, is.null, logical(1)))
groups <- rep(1, length(nonNull))
res <- do.call(rbind, res[nonNull])
res$rt <- -1
binned <- xcms::do_groupPeaks_mzClust(res, groups, minFraction = 0.2, ppm = 5)

binned <- xcms::do_groupChromPeaks_density(res, groups, binSize = 0.02, minFraction = 0.2)
peakIdx <- binned$peakidx
binned <- binned[, -which(colnames(binned) == "peakidx")]
binned[binned$mzmed > 372 & binned$mzmed < 373, ]

res[peakIdx[[11650]], ]
df <- res[peakIdx[[4782]], ]
df[order(-df$mz), ]
names(fileList)[12]

bin <- 0.05
spectraList <- sumR:::xcmsSpectraBinning(fileList, binSize = bin)
se <- sumR:::xcmsCellBinning(spectraList, binSize = bin)



binsizes <- c(0.1, 0.05, 0.02, 0.01, 0.005)
for (bin in binsizes) {
  spectraList <- sumR:::xcmsSpectraBinning(fileList, binSize = bin)
  se <- sumR:::xcmsCellBinning(spectraList, binSize = bin)
  se$Datetime <- 1:ncol(se)
  se$Type <- "SAMPLE"
  se$Type[grep("BLANK", colnames(se))] <- "BLANK"
  assay(se, "Value") <- assay(se, "Area") - assay(se, "Noise")

  df <- expToCombined(se, assay = "Value")
  p <- ggplot(df, aes(x = Aliquot, y = Value, fill = Type)) +
    ggplot2::scale_y_log10() +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ggsave(sprintf("heatmap_Value_binsize_%s.png", bin))

  df <- expToCombined(se, assay = "Area")
  p <- ggplot(df, aes(x = Aliquot, y = Area, fill = Type)) +
    ggplot2::scale_y_log10() +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  ggsave(sprintf("heatmap_Area_binsize_%s.png", bin))

  sub <- se[rowData(se)$mzmin > 372 & rowData(se)$mzmax < 373, ]
  rownames(sub) <- round(rowData(sub)$mzmed, 4)
  p <- plotHeatmap(sub, assay = "Area", main = sprintf("Peak Areas, BinSize: %s", bin), xlab = "Sample", ylab = "m/z bin median")
  htmlwidgets::saveWidget(p, file =  sprintf("Area_BinSize_%s.html", bin))

  p <- plotHeatmap(sub, assay = "Value", main = sprintf("Peak Values, BinSize: %s", bin), xlab = "Sample", ylab = "m/z bin median")
  htmlwidgets::saveWidget(p, file =  sprintf("Value_BinSize_%s.html", bin))
}



names(fileList)[29]
l <- fileList[[29]]
non_nulls <- !vapply(l, is.null, logical(1))
l <- l[non_nulls]
df <- data.frame(sample = rep(which(non_nulls), sapply(l,
                                                       nrow)), do.call(rbind, l))
df$rt <- -1
n <- length(which(non_nulls))
groups <- rep(1, n)
res <- suppressMessages(xcms::do_groupChromPeaks_density(df,
                                                         groups, binSize = binSize, minFraction = 0))
res$i <- vapply(res$peakidx, function(idx) sum(df$i[idx]),
                double(1))
res$noise <- vapply(res$peakidx, function(idx) sum(df$Noise[idx]),
                    double(1))
res


df$silhouette <- 0
for (n in 2:(nrow(res) - 1)){
  bin1 <- df[unlist(res[n - 1,"peakidx"]), ]

  ids <- unlist(res[n,"peakidx"])
  bin2 <- df[ids, ]

  bin2 <- df[unlist(res[n + 1,"peakidx"]), ]

  df$silhouette[ids] <- vapply(1:nrow(bin2), function(i){
    a <- mean(abs(bin2$mz[i] - bin2$mz[-i]))
    b <- mean(c(abs(bin2$mz[i] - bin1$mz), abs(bin2$mz[i] - bin3$mz)))
    (b - a) / max(a, b)
  }, double(1))
}

n <- 3342
bin1 <- df[unlist(res[n - 1,"peakidx"]), ]

ids <- unlist(res[n,"peakidx"])
bin2 <- df[ids, ]

bin2 <- df[unlist(res[n + 1,"peakidx"]), ]

silhouettes <- vapply(1:nrow(bin2), function(i){
  a <- mean(abs(bin2$mz[i] - bin2$mz[-i]))
  b <- mean(c(abs(bin2$mz[i] - bin1$mz), abs(bin2$mz[i] - bin3$mz)))
  (b - a) / max(a, b)
}, double(1))

silhouettes
res[res$mzmed > 372 & res$mzmed < 373, ][4,]

n <- 4
ids <- unlist(res[n,"peakidx"])
bin1 <- df[ids, ]

ids <- unlist(res[n + 1,"peakidx"])
bin2 <- df[ids, ]

ids <- unlist(res[n + 2,"peakidx"])
bin3 <- df[ids, ]


df$silhoutte[ids] <- vapply(1:nrow(bin2), function(i){
  a <- mean(abs(bin2$mz[i] - bin2$mz[-i]))
  b <- mean(c(abs(bin2$mz[i] - bin1$mz), abs(bin2$mz[i] - bin3$mz)))
  (b - a) / max(a, b)
}, double(1))




lapply(spectraList, function(df){
  df[df[, 1] < 728.6 & df[, 1] > 728.53, ]
})

samplePCA(sub)

sub <- se[rowData(se)$mzmin > 372.2 & rowData(se)$mzmax < 372.25, ]
sub <- se[1,]

t(assay(sub, "Area") / assay(sub, "Noise"))

sum(is.na(assay(se, "Area"))) / length(c(t(assay(se))))




idx <- assay(se, "Value") <= 0
assay(se, "Value")[idx] <- NA




library(ggplot2)





imputation(se, useAssay = "Area", noise = 1000)

data <- as.data.frame(assay(se, "Value"))

data[is.na(data)] <- runif(sum(is.na(data)), min = 1, max = 100)
assay(se, "Imputed") <- log2(data)

sub <- se[rowData(se)$mzmin > 372 & rowData(se)$mzmax < 373, ]
rownames(sub) <- round(rowData(sub)$mzmed, 4)
p <- plotHeatmap(sub, assay = "Area", main = sprintf("BinSize: %s", 0.05), xlab = "Sample", ylab = "m/z bin median")

ggsave(p)
a <- plotly::to_basic(p)
htmlwidgets::saveWidget(p, file =  sprintf("BinSize_%s.html", 0.05))

df <- reshape2::melt(assay(sub, "Area"))



# PC(14:0/14:1) seems good
assay(se, "Value")[rowData(se)$mzmax < 675.5 & rowData(se)$mzmin > 675.4, ]
assay(se[5618, ], "Value")

# SM(d18:1/16:0) Looks good, but quantification might be a bit off
rowData(se)[rowData(se)$mzmed < 702.58 & rowData(se)$mzmed > 702.54, ]
assay(se[6077, ], "Value")

# SM(d18:1/18:1)
rowData(se)[rowData(se)$mzmed < 728.6 & rowData(se)$mzmed > 728.5, ]
assay(se[6077, ], "Value")




head(df)
write.table(df, file = "sumR_test.tsv", sep = "\t", row.names = F)

#' @title Convert to a data frame in the long format
#' @param rowIndex Name of the rownames
#' @param colIndex Name of the columnnames
#' @export
expToCombined <- function(se, rowIndex = "Compound", colIndex = "Aliquot", assay = "Area"){
  mat <- as.matrix(assay(se, assay))
  x <- reshape2::melt(mat)
  colnames(x) <- c(rowIndex, colIndex, assay)

  df <- as.data.frame(cbind(x, colData(se)[x[, colIndex], ], rowData(se)[x[, rowIndex], ]))
  df <- df[order(df[, rowIndex], df[, colIndex]), ]
  rownames(df) <- 1:nrow(df)
  df
}

library(sumR)
files <- list.files(file.path(r"(F:\Myrthe\mzml)"), full.names = T)[1:34]
df <- DataFrame(Type = rep("SAMPLE", length(files)),
                row.names = tools::file_path_sans_ext(basename(files)))
df$Type[grep("BLANK", basename(files))] <- "BLANK"

df$Treatment <- c(
  rep("DMSO", 7),
  rep("TMX", 6),
  "DMSO",
  rep("TMX", 3),
  "DMSO",
  "TMX",
  rep("DMSO", 9),
  rep("TMX", 6)
)

exp2 <- extractPeaks(files, massWindow = c(100, 200), polarity = "-", centroid = T, cores = 16) %>%
  spectraAlignment(method = "binning", ppm = 5, nPeaks = 5) %>%
  cellAlignment(method = "binning", ppm = 5, nCells = 10, cellData = df, phenotype = "Treatment")

saveRDS(exp2, file = 'Aligned_negative.RDS')

saveRDS(massDefectFilter(exp), file = "Mass_defect_positve.RDS")
saveRDS(massDefectFilter(exp2), file = "Mass_defect_negative.RDS")
saveRDS(massDefectFilter(exp[, exp$Treatment == "DMSO"]), file = "Mass_defect_positive_DMSO.RDS")
saveRDS(massDefectFilter(exp[, exp$Treatment == "TMX"]), file = "Mass_defect_positive_TMX.RDS")

saveRDS(massDefectFilter(exp2[, exp2$Treatment == "DMSO"]), file = "Mass_defect_negative_DMSO.RDS")
saveRDS(massDefectFilter(exp2[, exp2$Treatment == "TMX"]), file = "Mass_defect_negative_TMX.RDS")

saveRDS(massDefectFilter(exp[, exp$Type == "BLANK"]), file = "Mass_defect_positive_BLANK.RDS")
saveRDS(massDefectFilter(exp[, exp$Type == "SAMPLE"]), file = "Mass_defect_positive_SAMPLE.RDS")

saveRDS(massDefectFilter(exp2[, exp2$Type == "BLANK"]), file = "Mass_defect_negative_BLANK.RDS")
saveRDS(massDefectFilter(exp2[, exp2$Type == "SAMPLE"]), file = "Mass_defect_negative_SAMPLE.RDS")

library(sumR)
exp <- readRDS("../Aligned_positive.RDS")
exp <- exp[, exp$Type == "BLANK"]
print(exp)
exp <- exp[rowSums(is.na(assay(exp))) != ncol(exp), ]
print(exp)



filterCells(exp)


exp <- extractPeaks(files, massWindow = c(100, 200), polarity = "+", centroid = T, cores = 16) %>%
  spectraAlignment(method = "binning", ppm = 5, nPeaks = 5) %>%
  cellAlignment(method = "binning", ppm = 5, nCells = 10, cellData = df, phenotype = "Treatment") %>%
  blankSubstraction(blankThresh = 5, nSamples = 10, removeBlanks = TRUE) %>%
  massDefectFilter() %>%
  imputation(method = "noise", noise = 100, seed = 42) %>%
  isotopeTagging(corr = 0.8) %>%
  fragmentFilter(corr = 0.95, method = "spearman")




exp2 <- extractPeaks(files, massWindow = c(100, 200), polarity = "-", centroid = T, cores = 16) %>%
  spectraAlignment(method = "binning", ppm = 5, npeaks = 5) %>%
  cellAlignment(method = "binning", ppm = 5, ncells = 10, cellData = df, phenotype = "Treatment")  %>%
  blankSubstraction(blankThresh = 5, nSamples = 10, removeBlanks = TRUE) %>%
  massDefectFilter() %>%
  imputation(method = "noise", noise = 100, seed = 42) %>%
  isotopeTagging(corr = 0.8) %>%
  fragmentFilter(corr = 0.95, method = "spearman")


exp3 <- combineExperiments(exp, exp2)

exp3 <- exp3 %>%
  shapiroTest() %>%
  leveneTest() %>%
  welchTest() %>%
  foldChange()


metadata(exp3)$phenotype <- "Treatment"

exp3 <- leveneTest(exp3)





exp3 <- exp3[rowData(exp3)$leveneTest$unequal_variance, ]
exp3

exp3 <- keepVariableFeatures(exp3)

samplePCA(exp3)

compoundPCA(exp3)
screePCA(exp3)
plotUMAP(exp3, components = 4)


exp3 <- exp3 %>%
  generateModel(modelName = "rf", seed = 42, cv = 5, ratio = 0.632) %>%
  generateModel(modelName = "glmnet", seed = 42)

rowData(exp[rowData(exp)$mz > 372 & rowData(exp)$mz < 373, ])

model(exp3, "rf")
as.data.frame(rowData(exp3[varImportance(exp3, 2)$Compound, ]))

head(rowData(exp3[varImportance(exp3)$Compound, ]))
head(assay(exp3[varImportance(exp3)$Compound, ]))

as.data.frame(df)
plot(model(exp3)$varImp)


as.data.frame()

head(assay(exp3)[order(rowData(exp3)$foldChange$log2fc_DMSO_TMX), ])

df[df$mz > 327 & df$mz < 328, ]
rownames(df) == "1385"

saveRDS(model(exp3), file = "model.RDS")

Rdisop::decomposeIsotopes(c(147.0529,148.0563), c(100.0,5.561173), maxElements = "C99")


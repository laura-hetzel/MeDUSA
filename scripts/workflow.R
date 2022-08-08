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

exp <- prepareFiles(files, massWindow = c(100, 200), polarity = "+", centroid = T, cores = 16) %>%
  spectraAlignment(method = "binning", ppm = 5, npeaks = 5) %>%
  cellAlignment(method = "binning", ppm = 5, nCells = 10, cellData = df, phenotype = "Treatment") %>%
  blankSubstraction(blankThresh = 5, nSamples = 10, removeBlanks = TRUE) %>%
  massDefectFilter(mz_MD_plot = F) %>%
  imputation(method = "noise", noise = 100, seed = 42) %>%
  fragmentFilter(corr = 0.95, method = "spearman") #%>%
  #isotopeTagging()

exp

exp2 <- prepareFiles(files, massWindow = c(100, 200), polarity = "-", centroid = T, cores = 16) %>%
  spectraAlignment(method = "binning", ppm = 5, npeaks = 5) %>%
  cellAlignment(method = "binning", ppm = 5, ncells = 10, cellData = df, phenotype = "Treatment")  %>%
  blankSubstraction(blankThresh = 5, nSamples = 10, removeBlanks = TRUE) %>%
  massDefectFilter(mz_MD_plot = F) %>%
  imputation(method = "noise", noise = 100, seed = 42) %>%
  fragmentFilter(corr = 0.95, method = "spearman") #%>%
  #isotopeTagging()

exp3 <- combineExperiments(exp, exp2)

exp3 <- exp3 %>%
  shapiroTest() %>%
  leveneTest() %>%
  welchTest() %>%
  foldChange()



exp2
exp

m <- as.matrix(assay(exp2))
n <- 6
cors <- sapply(1:nrow(m)[-n], function(i){
  cor(m[n,], m[i,])
})
rowData(exp2[cors > 0.8,])

assay(exp2[cors > 0.8,])
cors <- cor(t(m), y = t(m), method = "spearman")
df <- stack(cors)
df <- df[df$value > 0.8, ]
df <- df[as.integer(df$row) < as.integer(df$col), ]

mzs <- rowData(exp2)$mz
diff <- mzs[as.integer(df$col)] - mzs[as.integer(df$row)]
df <- df[diff < 1.1 & diff > 1, ]
df
i <- 4
assay(exp2[c(as.integer(df$row[i]), as.integer(df$col[i])), ])

exp2

l <- Rdisop::decomposeMass(x, ppm = 5, maxisotopes = 5, maxElements = "C10")
l

rowData(exp)


metadata(exp3)$phenotype <- "Treatment"

exp3 <- leveneTest(exp3)





exp3 <- exp3[rowData(exp3)$leveneTest$unequal_variance, ]
exp3

exp3 <- keepVariableFeatures(exp3, top = 100)

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


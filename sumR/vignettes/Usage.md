---
title: "Usage"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Usage}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
# bibliography: references.bib
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7, 
  fig.height = 4,
  eval = F
)

```

*This package is under heavy development and therefore is subject to change.*

## Introduction in sumR

![](images/sumr.png)

`sumR` (**s**ingle cell **u**ntargeted **m**etabolomics in **R**) is an R package designed for the analysis of single cell untargeted metabolomics datasets. Due to the small volume of cells and other challenges, conventional LC-MS strategies aren't suitable to analyse individual cells. This can be solved by using DI-MS, but brings its own challenges. `sumR` has been designed as a suite of functions that can deal with DI-MS data across many cells. This document functions as a showcase of a pipeline constructed with `sumR`.

## Pre-processing

sumR starts by pre-processing the raw data to features by extracting peaks and aligning these across spectra and cells. Unlike LC-MS, here we define a peak as a m/z -- intensity combination, and features as peaks that have been aligned across cells. This allows us to compare cells by e.g. filtering possible contaminants and/or performing statistical modelling. Consequently, features are defined using a min -- max m/z range across cells. Be aware that due to the acquisition method, sumR does not consider retention time of samples.

To start a sumR pipeline, we load the package.

```{r setupEcho, echo = T, eval = F}
library(sumR)
```

```{r setup, echo = F, include=FALSE}
suppressWarnings(suppressPackageStartupMessages(library(sumR)))
```

### Converting vendor format
SumR requires that the data is in `mzML` file format. In case you have vendor-format files, you can use the `rawToMzml` function to convert the files to mzML files. This function requires Proteowizard to be installed on your system, which you can find [here](https://proteowizard.sourceforge.io/). Alternatively, you can use Proteowizard without sumR to perform the conversion and skip this section.

After obtaining mzML files, we set the directory where our data is located. The path to the directory needs to be the full path, which can be assured with the `file.path()` function.

The rawToMzml function requires the defined input directory and an output folder to store the mzML files. It is recommended to use an empty or non-existing folder to ensure the output folder only contains mzML files. The function returns the file paths of the newly created mzML files.

```{r convert}
dir <- file.path(r"(F:\Myrthe\raw)")
files <- rawToMzml(dir, output = file.path(r"(F:\Myrthe\mzml)"))
```

### Defining Metadata

To make use of the full pipeline, it is mandatory to supply a `data.frame` of metadata. To ensure proper integration, ensure that the row names are equal to the base name of the file, without file extension. E.g. a file called `sample001.mzML` would translate into `sample001` as a proper row name. This name is matched to the data file names to match data with the supplied metadata.

```{r metadata2, eval = T, echo = F}
# knitr::kable(head(df), caption = "Start of the metadata file") 
# knitr::kable(tail(df), caption = "End of the metadata file")
```

Columns may contain metadata about each of the samples/cells. Here we prepared metadata with three additional columns: `Sample.number`, `phenotype`, and `Type`. The first column contains the number of each cell, the second contains the cell differentiation status and the third columns contains the sample type, here either `SAMPLE` or `BLANK`. We will use this phenotype during modelling after post-processing. The table below shows the contents of the metadata.



Now that we've extracted features from our data, we need to post-process the data in order to remove any contaminants, isotopes, adducts, and/or other artifacts. Due to the low cell volume, we are restricted from using pooled QCs. However, we can use Blanks and/or Lab QCs.

```{r}
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
knitr::kable(head(df), caption = "Start of the cellData metadata") 
```

### Extracting peaks
Peaks can be extracted from the mzML file using the `extractPeaks()` function. This function converts profile-mode data to centroid-mode by doing the following steps:

1. Extract scans within a massWindow (e.g. a scan of 200-350 dalton has a massWindow of 150 dalton). Defaults to `c(0, Inf)`
2. Extract scans of either positive (+) or negative (-) polarity. 
3. Smooth these scans using a Savitz-Golay filter
4. Centroid by picking local maxima

By setting `cores` to any value higher than 1, each scan is centroided in parallel using that amount of cores. The function returns a list of data frames, one for each cell. 

```{r prepareScans}
peaks <- extractPeaks(
  files, 
  massWindow = c(100, 200), 
  polarity = "+", 
  centroid = T, 
  cores = 8
)
```

We can plot the result of this function using a `spectrumPlot`. This plots the centroided peaks of a single scan by providing the filename and a scan number:

```{r spectrumPlot}
spectrumPlot(peaks, file = 1, scan = 1)
```

### Filtering Peaks

After inspection, we can filter noisy peaks by a minimum intensity value using the `peakFilter()` function. It returns a similar list of data frames, but only peaks with at least the given intensity are kept. We can view the difference by plotting the results in a `spectrumPlot`:

```{r peakPlot}
peakFilter(peaks, intensity = 1e3) %>%
  spectrumPlot(file = 1, scan = 1)
```

To save the peaks with our chosen intensity filter, we can override the results in the `peaks` variable.
```{r peakFilter}
peaks <- peakFilter(peaks, intensity = 1e3)
```

### Spectra Alignment

After we've detected our peaks, we need to identify which peaks should be considered equal across spectra within a cell. We use the function `spectraAlignment()` to align masses with a given `ppm` (Defaults to `5`). Alignment can be done using either the `binning` (default), `clustering`, or `density` method. While chromatography cannot be used to identify compounds, we do expect to find a compound multiple times across scans. This is adjusted using the `nPeaks` parameter. Here, we set the value of nPeaks to `5` to discard peaks that are likely to be false positives and should not be considered novel compounds:

```{r spectraAlignment}
spectra <- spectraAlignment(
   peaks, 
   method = "binning", 
   ppm = 5, 
   nPeaks = 5
)
knitr::kable(spectra, caption = "Compounds aligned across spectra") 

```

#### Inspecting Spectra Shifts
To inspect if the alignment did not over-correct, we can use the function `spectraShiftPlot()`. It shows the difference between the original mass and the mass shift after alignment. Changing the parameters `nPeaks` and `ppm` of the `spectraAlignment()` function will affect this mass shift. Since the alignment is based on ppm, we expect this mass shift to increase with increasing masses. The `cell` parameter either takes the cell name as defined in the metadata or the cell number. 

```{r spectraShift}
spectraShiftPlot(spectra, cell = 1)
```


### Cell Alignment

After the spectra within cells are aligned, we repeat the alignment between cells. Our arguments are similar to `spectraAlignment()`, however, the function `cellAlignment` returns a `SummarizedExperiment` object instead. Since this object is used throughout the sumR workflow, you can read more about interacting with a SummarizedExperiment object [here](https://bioconductor.org/help/course-materials/2019/BSS2019/04_Practical_CoreApproachesInBioconductor.html).

Other differences include the `nCells` parameter instead of `nPeaks`, which sets the threshold for the minimum number of cells a compound should be found. the nCells parameter is highly dependent on the experimental design, as it affects the biological variation between cells. The metadata of each cell, `cellData`, can be set here as well. Omitting this will result in a warning, as it is necessary to continue the post-processing pipeline. Lastly, the `phenotype` parameter can be set and refers to the column in the cellData that will be used for statistical tests and modelling.

```{r cellBinning}
exp <- cellAlignment(
  spectra, 
  method = "binning", 
  ppm = 5, 
  nCells = 10,
  cellData = df,
  phenotype = "Treatment"
)
```

#### Inspecting Cell Shifts

We can inspect the mass shift between cells that occured during alignment. 

```{r cellShift}
cellShiftPlot(exp)
```

## Post-processing
With aligned features across cells, the pre-processing is complete. However, not all aligned features are actually novel compounds. This is due to artifacts, fragments, isotopes, and other non-biological compounds. The following functions are part of the post-processing pipeline in sumR to remove these unwanted features. All these functions are optional (and can be executed in any order), but are recommended for most uses.

### Blank substracton

When Blank samples are measured, these can be used to filter any features. `blankThresh` indicates the fold-change needed to exceed for a sample to be kept. By default, samples must have an intensity of at least 5x higher than the blank samples in order to be kept. This ensures that noisy peaks are removed. The `nSamples` parameter determines the number of samples that the blank threshold should exceed in order for a compound to be removed. The default for nSamples is infinite, thus no features will be removed. To remove the blank samples after filtering compounds, set the `removeBlanks` parameter to TRUE (the default)

```{r Blanks}
exp <- blankSubstraction(exp, blankThresh = 5, nSamples = 10, removeBlanks = TRUE)
```

### Mass Defect Filter
As described by [McMillan et al.](10.1186/s13321-016-0156-0), calculating the mass defect can be used to identify potential salt clusters and other non-biological compounds. This is implemented in the `massDefectFilter` which can plot the mass versus mass defect.
```{r}
exp <- massDefectFilter(exp, plot = TRUE)
```

### Missing value imputation

Due to the dropout effect, single cell datasets can contain over 50% of NA values which will need to be imputed for statistical analysis. sumR currently has two methods for imputation, determined by the parameter `method`. This parameter can either be set to `"noise"` or `"saver"`. The former is a noise-based random imputation method. It will generate a random number between 1 and the value of `noise`, which defaults to 100. The saver method is model-based that assumes a poisson distribution and will generate values that are similar to highly correlated cells.

```{r imoputation}
exp <- imputation(exp, method = "noise", noise = 100, seed = 42)
```

### Isotope identification

When interested in novel compounds, isotopes and adducts are a source of unwanted features. The function `isotopeTagging()` can detect isotopes and adducts using known ratios of mass and intensities of atoms. Calling this function will store detailed results about isotopes and adducts in the `rowData` slot of the object. 

```{r isotopes}
exp <- isotopeTagging(exp)
```

### Fragment filtering

Another source of unwanted features is fragmented compounds. These are features that show highly correlated intensity values and are generally removed during post-processing. Here we use the function `fragmentFilter()` to remove highly correlated features. The threshold of correlation can be determined by the parameter `corr`, which is set here at 95%. Available methods for fragment filtering are `pearson` and `spearman`. Here spearman correlation is used, as for datasets with a lot of imputed values pearson correlation might overestimate the correlation.

```{r fragments}
exp <- fragmentFilter(exp, method = "spearman", corr = 0.95)
```

## Processing workflow

We can combine the functions of the pre- and post-processing in a pipeline using the 'pipe' notation of the `dplyr` package. Here we are using the same settings as the previous steps except for positive polarity instead of negative polarity.

```{r pipeline}
# Pre-processing

exp2 <- extractPeaks(files, massWindow = c(100, 200), polarity = "+", centroid = T, cores = 16) %>%
  spectraAlignment(method = "binning", ppm = 5, nPeaks = 5) %>%
  cellAlignment(method = "binning", ppm = 5, nCells = 10, cellData = df, phenotype = "Treatment")

# Post-processing
exp2 <- exp2 %>%
  blankSubstraction(blankThresh = 5, nSamples = 10, removeBlanks = TRUE) %>%
  massDefectFilter() %>%
  imputation(method = "noise", noise = 100, seed = 42) %>%
  isotopeTagging() %>%
  fragmentFilter(method = "spearman", corr = 0.95)
```

### Combining experiments

We can combine the features found in both polarities with the `combineExperiments()` function. This function takes in any number of `SummarizedExperiment` objects with the assumption of equal 'colData'. The function combines these objects row-wise, essentially adding features together and returns a new SummarizedExperiment object.

```{r comb}
expModel <- combineExperiments(exp, exp2)
expModel
```


## Modelling & Statistics

`sumR` supports building models to determine phenotypes that might be discriminatory. The package includes several univariate tests, as well as multivariate tests using models in the [caret package](https://topepo.github.io/caret/). Finally we can inspect the results using PCA, UMAP and other plots.

### Univariate statistical tests

sumR supports the shapiro-wilk test for normality, levene test for difference in variance between groups, and welch T-test for significant difference testing between groups. These tests will test the features (row-wise), not samples. Finally, while not a statistical test, the fold change between the given phenotype can be calculated. 

```{r statTests}
expModel <- expModel %>%
  shapiroTest() %>%
  leveneTest() %>% 
  welchTest() %>%
  foldChange()

  # TODO: Plot each of these in appropriate plots
```

### Filter non-variable features

Due to the drop-out effect, a large amount of values need to be imputed. Generally, this causes these features to have a low variance and consequently a low impact on models. With the function `keepVariableFeatures()` we can subset the found features based on their variance, determined by the `leveneTest()` function. Using the `plotFeatureSds()` function we can inspect the standard deviation of each feature in a barplot.

```{r variableFeatures}
expModel <- keepVariableFeatures(expModel)
plotFeatureSds(expModel)
```


### Data inspection using PCA & UMAP

The results can be inspected using several plot types. `samplePCA()` is a function that plots a Principle Components Analysis (PCA) on the samples. In contrast, the `compoundPCA()` function plots a PCA for the compounds instead. Next, the `screePCA()` function plots a barplot with the variance explained for each Principle Component (PC). Finally, we can use the `plotUMAP()` function to plot a UMAP after doing a PCA first. This function takes in the number of PCs to construct the UMAP.

```{r pcaPlots}
samplePCA(expModel)
compoundPCA(expModel)
screePCA(expModel)
plotUMAP(expModel, components = 10)
```

### Model generation

As stated before, sumR utilizes the caret package to generate models. The 'Available Models' section in the caret documentation indicates which models are available (you may need to install dependencies for some models). The `generateModel()` function has the following workflow:

1.   Data partitioning in a train/test set
2.  Set the train control with scaling and centering
3.  Perform training with cross-validation
4.  Predict on the partitioned test set
5.  Produce a confusion matrix and determine important variables

the generateModel function takes in the following parameters:

- exp: The SummarizedExperiment object
- assay: The assay to use (should be an assay without missing values)
- modelName: Name of the model to use as defined by `caret`
- folds: Number of folds to use for cross-validation
- ratio: Ratio of train-test data split
- seed: Used for reproducible results

```{r models}
expModel <- generateModel(
    exp = expModel,
    assay = "Imputed",
    modelName = "rf", 
    folds = 5, 
    ratio = 0.632,
    seed = 42
) 
```

All generated models are stored in the metadata of the SummarizedExperiment object under `model`. 

### Model Assessment

Using the `model()` function, information about a given model can be accessed. It will return a list with the following entries: `train`, `test`, `model`, `prediction`, `varImp`, and `confMatrix`. These can be used to assess the quality of the model for the given phenotype. Here we are printing the confusion matrix.

```{r assessment}
print(model(expModel, "rf")$confMatrix)

# TODO: Plot the ROC, Varimportance, scatterplot (PCA?) with predictions, etc. 
```

### Cross validation

We can also inspect the performance of each model by plotting the accuracy during cross validation. This is done by the `plotCrossValidation()` function with a given model name.

```{r CV}
plotCrossValidation(expModel, "rf")
```

# SessionInfo

```{r sessionInfo}
sessionInfo()
```

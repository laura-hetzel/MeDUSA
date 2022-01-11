# README

This is the README for the `sumR` package.

## Install R
1. Download and install R and Rstudio Desktop from:

https://cran.r-project.org/bin/windows/base/ and https://www.rstudio.com/products/rstudio/download/

2. Start Rstudio in administrator mode (right-click, run as Administrator). This ensures that the packages in the next step can be installed.

## Install package

1. Create a new project. This can be done by clicking on File -> New Project -> New Directory -> New Project. Pick an appropriate name for your project so you can find it on your computer, e.g. '`sumR`'.

2. Create a new file by clicking on File -> New File -> R Script

3. Install required packages using the following commands.

```{r eval = F, echo = T}
if (startsWith(.libPaths()[1], "\\")) {
    .libPaths(sort(.libPaths(), decreasing = T))
}

install.packages("remotes")
```

4. Request a personal access token in Gitlab. Go to: https://gitlab.com/-/profile/personal_access_tokens and add a Token name. Ensure that 'read_api' scope is selected. Click on 'Create personal access token' to retrieve the token. Copy this code, you will need it in the next step. If you already have a token from another project, you can use that token instead.

5. Use the following code to install the `sumR` package. Replace the value "abc123" with your personal access token. The `sumR` package should install automatically. 


```{r eval = F, echo = T}
if (startsWith(.libPaths()[1], "\\")) {
    .libPaths(sort(.libPaths(), decreasing = T))
}

my_token <- "abc123"
remotes::install_gitlab(
  repo = "lacdr-abs/sum-r",
  subdir = "sumR",
  auth_token = my_token,
  dependencies = TRUE,
  force = TRUE,
  type = "binary",
  upgrade = FALSE,
  build_manual = TRUE,
  build_vignettes = TRUE
)
```

# Basic usage

For basic usage of `sumR`, you can execute the following code:

```{r eval = F, echo = T}
vignette(topic = "Usage", package = "sumR")
```

# Troubleshooting

If you encounter issues during installation or execution, please submit an issue in the issue tracker here: https://gitlab.com/lacdr-abs/sum-r/-/issues

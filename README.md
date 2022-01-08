# README

This is the README for the `sumR` package.

## Install on Windows
1. Download and install R and Rstudio Desktop from:

https://cran.r-project.org/bin/windows/base/ and https://www.rstudio.com/products/rstudio/download/

2. Start Rstudio in administrator mode (right-click, run as Administrator). This ensures that the packages in step 3 can be installed.

3. Create a new project. This can be done by clicking on File -> New Project -> New Directory -> New Project. Pick an appropriate name for your project so you can find it on your computer, e.g. '`sumR`'.

4. Create a new file by clicking on File -> New File -> R Script

5. Install required packages using the following commands. This will create a project-specific library where packages are installed. 

```{r eval = F, echo = T}
if (startsWith(.libPaths()[1], "\\")) {
    .libPaths(sort(.libPaths(), decreasing = T))
}

install.packages("remotes")
```

6. Request a personal access token in Gitlab. Go to: https://gitlab.com/-/profile/personal_access_tokens and add a Token name. Ensure that 'read_api' scope is selected. Click on 'Create personal access token' to retrieve the token. Copy this code, you will need it in the next step

7. Use the following code to install the `sumR` package. Replace the value "abc123" with your personal access token.  


```{r eval = F, echo = T}
if (startsWith(.libPaths()[1], "\\")) {
    .libPaths(sort(.libPaths(), decreasing = T))
}
my_token <- "abc123"
remotes::install_gitlab(
  repo = "lacdr-abs/sum-r",
  subdir = "package",
  auth_token = my_token,
  dependencies = TRUE,
  force = TRUE,
  type = "binary",
  upgrade = FALSE
)
```

`sumR` should install automatically. If you encounter issues during installation, please submit an issue in the issue tracker here: https://gitlab.com/lacdr-abs/sum-r/-/issues

# Basic usage

# Troubleshooting

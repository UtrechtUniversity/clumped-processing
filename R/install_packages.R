# typically, you do NOT need to run this!
# all packages are installe in C:/R/win-library/4.1/

# install the base packages
install.packages(c(
  # meta-package that has most of what we need
  "tidyverse",
  # work with strings
  "glue",
  # show informative messages when using tidyverse functions
  "tidylog",
  # the workflow/pipeline package
  "targets",
  # for rolling offset correction and empirical transfer function
  "slider", "zoo",
  # read in excel and csv files etc.
  "readxl", "writexl", "readr",
  # efficient filetypes
  "qs", "fst",
  # plotting helpers
  "plotly", "gghighlight", "ggforce", "patchwork",
  "ggnewscale",
  # to install development packages
  "devtools", "remotes",
  # optional extensions to isoreader
  "feather", "openxlsx", "xml2", "BiocManager"
))

# install the development packages
BiocManager::install("rhdf5")
devtools::install_github("leeper/unf")
devtools::install_github("isoverse/isoreader")
devtools::install_github("isoverse/clumpedr")

library(tidyverse) # collection of packages dplyr, purrr, readr,
library(glue)
library(readxl)
# remotes::install_github("isoverse/isoreader")
library(isoreader)
# remotes::install_github("isoverse/clumpedr")
library(clumpedr)
library(targets)
library(slider)
library(plotly) # not needed for workflow
library(patchwork)
library(gghighlight)
library(ggforce)
# I use tidylog for interactive sessions for nice logs, but I call it explicitly when desired
library(gghighlight) # not needed for workflow

# this used to be nice, but may be confusing
sam <- scale_alpha_manual(values=c("TRUE" = 0.2, "FALSE" = 1))

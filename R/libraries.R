library(tidyverse) # collection of packages dplyr, purrr, readr,
library(glue)
library(readxl)
library(isoreader)
library(clumpedr)
library(targets)
library(slider)
library(plotly) # not needed for workflow
# I use tidylog for interactive sessions for nice logs, but I call it explicitly when desired
library(gghighlight) # not needed for workflow
sam <- scale_alpha_manual(values=c("TRUE" = 0.2, "FALSE" = 1))

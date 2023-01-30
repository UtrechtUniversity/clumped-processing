# you can run this script as a Background Job so you can continue using an
# interactive R session.
library(targets)
library(tidyverse)
targets::tar_make(
  names = c(
    
  # motu_scn_meta_update, 
  #  motu_meta_update,
  motu_export, motu_export_csv, motu_out, 
  
  # pacman_meta_update,
  # pacman_scn_meta_update,
  # pacman_export, pacman_export_csv, pacman_out,
  ),
)
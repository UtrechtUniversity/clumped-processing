# this just runs all targets
targets::tar_make()

# use the shortcut if you don't care about checking every dependency
targets::tar_make(motu_sample_level, shortcut = TRUE)
targets::tar_make(motu_out, shortcut = TRUE)


# all of the below should be removed once we solve #3

# our workflow currently requires manual intervention after this step
targets::tar_make(names=c(
  # the background scans
  motu_scn_meta_update, pacman_scn_meta_update,
  # the new rows for the metadata files
  motu_meta_update, pacman_meta_update, pacman_caf_meta_update))
# now update your metadata files by copying out/..._update.xslx to the metadata
# files on OneDrive

# now we can build the rest of the target
targets::tar_make(names=c(
  # to excel
  motu_export, pacman_export, pacman_caf_export,
  # to CSV
  motu_export_csv, pacman_export_csv, pacman_caf_export_csv,
  # rds copies
  motu_out, pacman_out, pacman_caf_out))
# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### TODO

- fix pacman workflow #6
  - etf_group is empty, causing slider::hop to have an NA for the .starts argument
  - d13C_PDB_mean is NA for more than half and rubbish for the rest
  - a lot of scans don't have the proper scan_group etc.: need to tweak ranges
    further back in time in pacman_scn_metadata_parameters file
  - the export_metadata_update function doesn't work for pacman, because it 
    already has all the Analysis numbers, but not all the file_id's in it...

### Added 

- 2022-12-19: gave a live demo to UU clumpy group, showing how to use the project.
- running_clumped-processing.Rmd notebook that walks through how everything works
- inspecting_clumped-processing.Rmd notebook to make figures and check output
- this changelog file
- 2023-01-23: runsimple.R file to run workflow as Background Job in RStudio

### Fixed

- 2023-01-16: Pacman
  - updated pacman_did_metadata_parameters with the newest version, made sure 
    the etf_grp column is formatted the same as in motu_metadata_parameters.
  - changed the root folder for pacman scans so that it finds all scan files.
  
- 2023-01-23: s44_init and r44_init were always NA for motu and pacman #8
  - Changed add_inits and extract_file_info functions (re-order, fix logic)


## [v1.0.0] - 2022-04-07

### Added
- copying notice
- GPL-3 Licence
- pacman processing code

### Changed
- Updated background correction from linear to 3rd order polynomial
# Clumped Isotope R Workflow
# written by Ilja J. Kocken https://orcid.org/0000-0003-2196-8718
# Copyright 2022 © Ilja Kocken

# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.

# libraries and options ---------------------------------------------------
library(targets)
library(tarchetypes)

# Note that we're using the development package ~clumpedr~, which I'm writing. Install it with:
# remotes::install_github("isoverse/clumpedr")

tar_option_set(
  packages = c("tidyverse",
               ## "readr",
               "readxl",
               "isoreader",
               "clumpedr",
               "slider"),
  workspace_on_error = TRUE
)

# options(crayon.enabled = FALSE)
options(clustermq.scheduler = "multicore")

source("R/functions.R")

# after how many days should we guarantee a re-run?
our_update_interval <- as.difftime(50, units = "days") 

# general and logbooks ---------------------------------------------------------
# These general targets contain accepted standard values and excel logbooks.
# Currently, the latter are not used in the pipeline.
list(
  tar_target(
    accepted_standard_values_file,
    "C:/Archive/OneDrive - Universiteit Utrecht/Archive/accepted_standard_values.csv",
    format = "file"
  ),
  
  tar_target(
    accepted_standard_values,
    read_csv(accepted_standard_values_file)
  ),
  
  tar_target(
    stdnames,
    c(paste0("ETH-", 1:4), paste0("IAEA-C", 1:2), "Merck")
  ),
  
  # logfiles currently not used
  tar_target(
    motu_log_file,
    "//ad.geo.uu.nl/GML/Rawdata/253pluskiel/logbook_253plus.xlsx",
    format = "file"
  ),
  
  tar_target(
    motu_log,
    readxl::read_excel(
      motu_log_file,
      sheet = "logbook  253plus",
      range = "A1:AB1000",
      col_types = c("date",
                    "date",
                    "text",
                    ## "text", "text", "text",
                    rep("guess", 25))
    ) |>
      mutate(
        datetime = paste(
          as.character(Date),
          as.character(`Time start prep (heat PP from May 2019, unless otherwise stated)`) |>
            str_replace("^1899-12-31 ", "")
        ) |>
          lubridate::as_datetime()
      )
  ),
  
  tar_target(
    motu_maintenance,
    readxl::read_excel(
      motu_log_file,
      sheet = "Maintenance  253plus",
      range = "A1:D1000",
      col_types = c("date", rep("guess", 3))
    )
  ),
  
  
  # motu pipeline -----------------------------------------------------------
  # MotU stands for Master of the Universe, and is our fanciest newest mass
  # spectrometer, the 253 plus with a Kiel-IV device.
  
  # This uses dynamic targets for all the specific files. This allows us to
  # process files independently and only combine them at the ETF level. We use
  # ~iteration = "list"~ to make dynamic targets per directory, so that
  # preparations only need to be read into R once.
  
  # These are the measurement files for the standards and the samples. That's 46
  # measurements per run/preparation/sequence.
  tar_target(
    motu_dids_paths_all,
    list_files("Kiel Raw Data",
               wd = "//ad.geo.uu.nl/GML/Rawdata/253pluskiel/Raw Data") |>
      file_info() |>
      remove_copies() |>
      batch_files(),
    # it now iterates over the directories
    iteration = "list",
    # run if it hasn't been run in 3 days, for now
    cue = tarchetypes::tar_cue_age(name = motu_dids_paths_all,
                                   age = our_update_interval)
  ),
  
  tar_target(
    motu_dids_paths,
    motu_dids_paths_all, # |>
    # this is to quickly play around with a subset
    ## vctrs::vec_c() |>
    ## vctrs::vec_slice(c(1:3, floor(length(.)/2) + c(-1,0,1), length(.) + c(-2, -1, 0))),
    iteration = "list"
  ),
  
  tar_target(
    motu_did_files,
    motu_dids_paths,
    format = "file",
    pattern = map(motu_dids_paths)
  ),
  
  # scn files
  # These are the background scans. We create 5 files per run, and they are used
  # to correct all the measurements that follow it until the next scans.
  tar_target(
    motu_scn_paths_all,
    list_files("253pluskiel", ".scn$",
               # unfortunately the BGs are kept in many folders at the root level
               wd = "//ad.geo.uu.nl/GML/Rawdata") |>
      file_info() |>
      remove_copies() |>
      batch_month(),
    iteration = "list",
    cue = tarchetypes::tar_cue_age(name = motu_scn_paths_all,
                                   age = our_update_interval)
  ),
  
  tar_target(
    motu_scn_paths,
    motu_scn_paths_all, ##  |>
    # small subset!
    ## vctrs::vec_c() |>
    ## vctrs::vec_slice(c(1:3, floor(length(.)/2) + c(-1,0,1), length(.) + c(-2, -1, 0))),
    iteration = "list"
  ),
  
  tar_target(
    motu_scn_files,
    motu_scn_paths,
    format = "file",
    pattern = map(motu_scn_paths)
  ),
  
  ## read in as isoreader files
  # The above only listed the files and cut them up into list chunks per run. Here
  # we read in the data in the files. This is quite slow and usually only needs to
  # happen once, unless we have an update in the ~isoreader~ package.
  tar_target(
    motu_dids,
    read_di(motu_did_files),
    pattern = map(motu_did_files),
    iteration = "list",
    format = "qs"
    #cue = tar_cue(command = FALSE)
  ),
  
  tar_target(
    motu_scn,
    read_scn(motu_scn_files),
    pattern = map(motu_scn_files),
    iteration = "list",
    format = "qs"
    #cue = tar_cue(command = FALSE)
  ),
  
  ## extract raw data
  ## This gets the raw data, i.e. individual cycles of intensities per mass, from
  ## the above files.
  
  tar_target(
    motu_raw,
    iso_get_raw_data(motu_dids, include_file_info = Analysis),
    #|>
    # this now iterates over the folders, so it won't have to re-run this
    # expensive function
    pattern = map(motu_dids),
    iteration = "list",
    format = "qs"
  ),
  
  tar_target(
    motu_scn_raw,
    iso_get_raw_data(motu_scn,
                     include_file_info = c(file_root, file_datetime)),
    pattern = map(motu_scn),
    iteration = "list",
    format = "qs"
  ),
  
  ## read in metadata
  ## These files hold the current metadata fixes with desired parameters for data
  ## processing.
  tar_target(
    motu_meta_file,
    "C:/Archive/OneDrive - Universiteit Utrecht/Archive/motu_metadata_parameters.xlsx",
    format = "file"
  ),
  
  tar_target(
    motu_metadata,
    readxl::read_excel(motu_meta_file, guess_max = 1e5) |>
      meta_fix_types() |> # TODO: switch to parse_col_types?
      tidylog::distinct(Analysis,
                        # there are some with unique file_id's but the same file contents
                        # file_id,
                        file_size, file_datetime, .keep_all = TRUE),
    format = "fst_tbl"
  ),
  
  tar_target(
    motu_scn_meta_file,
    "C:/Archive/OneDrive - Universiteit Utrecht/Archive/motu_scn_metadata_parameters.xlsx",
    format = "file"
  ),
  
  tar_target(
    motu_scn_meta,
    read_xlsx(
      motu_scn_meta_file,
      sheet = "data",
      guess_max = 2e3,
      col_types = c(
        "text",
        "text",
        "date",
        rep("numeric", 16),
        "text",
        "date",
        "text",
        "text",
        "numeric",
        "text",
        rep("numeric", 18),
        "logical",
        "logical",
        "logical",
        "text",
        "numeric",
        "text",
        "date",
        "text"
      ),
      na = c("", "NA")
    )
  ),
  
  ## process scans
  
  # TODO: import/export motu_scn_metadata so that I output all parameter columns
  tar_target(
    motu_scn_fix,
    motu_scn_raw |>
      nest_by(file_id, file_root, file_datetime) |>
      # this gets some metadata from the raw scan
      file_name_scn() |>
      # this is a way to create the metadata file for the first time:
      ## mutate(min = -500, max = 50000,
      ##        min_start_44 = 9.392386, min_end_44 = 9.395270,
      ##        min_start_45_49 = 9.424277, min_end_45_49 = 9.429723,
      ##        max_start = 9.464633, max_end = 9.468291) |>
      add_scan_info(
        motu_scn_meta,
        c(
          "min",
          "max",
          "min_start_44",
          "min_end_44",
          "min_start_45_49",
          "min_end_45_49",
          "max_start",
          "max_end",
          "manual_outlier",
          "fix_software",
          "scan_group_overwrite",
          "voltage_overwrite",
          "checked_by",
          "checked_date",
          "checked_comment"
        )
      )  |>
      fix_scan_meta() |>
      unnest(cols = c(data)) |>
      fix_motu_scans(),
    pattern = map(motu_scn_raw),
    iteration = "list",
    format = "qs"
  ),
  
  tar_target(
    motu_scn_mod,
    motu_scn_fix |>
      tidy_scans() |>
      flag_scan_ranges() |> # also gets rid of manual outliers
      flag_scan_capped() |>
      calculate_min_max() |>
      # this combines the scans of the same scan_group into one row
      calculate_scan_models(),
    ## unnest(data) |>
    pattern = map(motu_scn_fix),
    iteration = "list",
    format = "qs"
  ),
  
  tar_target(
    motu_scn_meta_update,
    export_scan_metadata(
      data = motu_scn_mod |>
        bind_rows() |>
        unnest(c(data, scan_files)),
      meta = motu_scn_meta,
      file = "out/motu_scn_metadata_update.xlsx"
    ),
    format = "file"
  ),
  
  ## clean up metadata, make file info
  
  # extracted because it's slow and never changes after reading it once
  tar_target(
    motu_file_info_raw,
    extract_file_info(motu_dids),
    pattern = map(motu_dids),
    iteration = "list",
    ## cue = tar_cue(command = FALSE),
    format = "qs"
  ),
  
  ## # quickly subset to date range for experimenting with bg factor
  ## tar_target(my_filter, motu_file_info_raw |>
  ##                    bind_rows() |>
  ##                    tidylog::filter(file_datetime > lubridate::ymd("2020-01-01"),
  ##                                    file_datetime < lubridate::ymd("2020-11-01"))
  ##            ),
  
  tar_target(
    motu_file_info,
    motu_file_info_raw |>
      rowwise() |>
      # this adds all the _overwrite columns and manual_outlier etc.
      # it also tries to get the Preparation number from the filename if it doesn't exist
      fix_metadata(motu_metadata, irms = "MotU-KielIV") |>
      # this then applies them to calculate identifier_1 etc.
      overwrite_meta(stdnames = stdnames) |>
      # this creates bg_group based on the file_datetime and the scan_datetime
      add_scan_group(motu_scn_mod |> bind_rows()) |>
      # this adds the parameters that are now in motu_metadata in stead of parms
      add_parameters(motu_metadata) |>
      rename(c("outlier_manual" = "manual_outlier")),
    pattern = map(motu_file_info_raw),
    iteration = "list",
    format = "qs"
  ),
  
  # this is a subset target so that the raw part only needs to be run when these
  # specific metadata are updated
  tar_target(
    motu_raw_file_info,
    motu_file_info |>
      bind_rows() |>
      select(
        Analysis,
        file_id,
        dis_min,
        dis_max,
        dis_fac,
        dis_rel,
        # cycle_filter
        bg_group,
        starts_with("scan_"),
        #starts_with("lm_"), #starts_with("slope_"),
        bg_fac,
        d13C_PDB_wg,
        d18O_PDBCO2_wg,
        Mineralogy
      ),
    pattern = map(motu_file_info),
    iteration = "list"
  ),
  
  tar_target(
    motu_badruns,
    motu_file_info |>
      bind_rows() |>
      find_bad_runs()
  ),
  
  tar_target(
    motu_meta_update,
    export_metadata(
      data = motu_file_info |>
        bind_rows(),
      meta = motu_metadata,
      file = "out/motu_metadata_update.xlsx"
    ),
    format = "file"
  ),
  
  ## raw deltas
  
  # Most of the computations have already landed in my clumpedr package
  # https://github.com/isoverse/clumpedr/, but we do have some tricks here that
  # I've found not to be general enough for sharing with the wider community, such
  # as offset correction. I've made the calls to ~clumpedr~ explicit with ~::~ so
  # that it is clear which functions are mainained in this repository and which
  # ones are in the other package.
  tar_target(
    motu_raw_deltas,
    motu_raw |>
      # somehow it's become a character
      # there are some without Analysis
      tidylog::mutate(
        Analysis = ifelse(
          "Analysis" %in% colnames(motu_raw),
          readr::parse_double(Analysis),
          NA_real_
        )
      ) |>
      # write a wrapper function for this so that the targets are simpler
      # TODO figure out how to loop over two separate lists of both raw and meta info
      add_info(
        motu_raw_file_info,
        c("dis_min", "dis_max", "dis_fac", "dis_rel")
      ) |>
      clumpedr::find_bad_cycles(
        min = dis_min,
        max = dis_max,
        fac = dis_fac,
        # TODO: get relative_to parms call to work based on dataframe itself
        relative_to = "init"
      ) |>
      filter_duplicated_raw_cycles() |>
      clumpedr::spread_match() |>
      add_info(.info = motu_raw_file_info,
               cols = c("bg_group", "bg_fac")) |>
      add_background_info(motu_scn_mod |>
                            bind_rows()) |>
      # TODO: use neighbouring scans before and after sample to get rid of scan noise?
      correct_backgrounds_scn(fac = .data$bg_fac) |>
      # remove the scan models because they take up a lot of memory as list columns
      select(-starts_with("lm_")) |>
      add_info(.info = motu_raw_file_info,
               c("d13C_PDB_wg", "d18O_PDBCO2_wg")) |>
      clumpedr::abundance_ratios(s44, s45_bg, s46_bg, s47_bg, s48_bg, s49_bg) |>
      clumpedr::abundance_ratios(
        r44,
        r45_bg,
        r46_bg,
        r47_bg,
        r48_bg,
        r49_bg,
        R45_wg,
        R46_wg,
        R47_wg,
        R48_wg,
        R49_wg
      ) |>
      clumpedr::little_deltas() |>
      add_info(motu_raw_file_info, c("Mineralogy")) |>
      add_R18() |>
      # TODO check if this works for dolomite samples, not sure if vectorized
      clumpedr::bulk_and_clumping_deltas(R18_PDB = .data$R18_PDB) |>
      # outlier on the cycle level now contains all the reasons for cycle outliers
      clumpedr::summarise_outlier(quiet = TRUE),
    # TODO: exclude values mass 54/48/49 < -490
    # TODO: decide whether p49 can be ignored here? I think so because we're doing it at sample level now
    ## add_info(motu_file_info |> bind_rows(), c("Analysis", "p49_crit")) |>
    ## clumpedr::find_R_flags() |>  # TODO: get rid of R_flags? do they find anything of value?
    pattern = map(motu_raw, motu_raw_file_info),
    iteration = "list",
    format = "qs"
  ),
  
  
  # nesting and summarising still happens within each folder, because this is
  # slow for the big db
  tar_target(
    motu_nested,
    motu_raw_deltas |>
      clumpedr::nest_cycle_data() |>
      summarize_d13C_d18O_D47(),
    pattern = map(motu_raw_deltas),
    iteration = "list",
    format = "qs"
  ),
  
  ## sample level summaries
  tar_target(
    motu_sample_level,
    motu_nested |>
      bind_rows() |>  # finally the data are rbinded into one big df!
      string_scan_files() |>
      add_remaining_meta(motu_file_info |> bind_rows()) |>
      # there are 24 measurements for motu that are so messed up
      # that they cannot be matched by both file_id and Analysis somehow
      tidylog::filter(!is.na(init_low)) |>
      clumpedr::find_init_outliers(
        init_low = .data$init_low,
        init_high = .data$init_high,
        init_diff = .data$init_diff
      ) |>
      summarise_cycle_outliers() |>
      mutate(outlier_param49 = param_49_mean > p49_crit |
               param_49_mean < -p49_crit) |>
      ## summarize_outlier() |>
      # try out conservative outlier selection
      mutate(
        outlier = outlier_noscan | outlier_nodelta |
          (!is.na(outlier_cycles) & outlier_cycles) |
          (!is.na(outlier_init) & outlier_init) |
          (!is.na(outlier_manual) & outlier_manual)
      ) |>
      # get rid of raw cycle data
      ## select(-cycle_data) |>
      tidylog::select(-where(is.list)) |>
      arrange(file_datetime) |>
      # get rid of duplicated rows
      tidylog::distinct(Analysis, file_id, file_size, s44_init, r44_init, .keep_all = TRUE) |>
      offset_correction_wrapper(acc = accepted_standard_values),
    format = "fst_tbl"
  ),
  
  tar_target(
    motu_temperature,
    motu_sample_level |>
      # there are many ways of calculating the ETF
      ## raw session
      clumpedr:::calculate_etf(
        raw = D47_raw_mean,
        exp = expected_D47,
        session = etf_grp,
        etf = etf,
        etf_coefs = etf_coefs,
        slope = etf_slope_grp,
        intercept = etf_intercept_grp
      ) |>
      ## offset corrected session
      clumpedr:::calculate_etf(
        raw = D47_offset_corrected,
        exp = expected_D47,
        session = etf_grp,
        etf = etf_off,
        etf_coefs = etf_coefs_off,
        slope = etf_slope_grp_off,
        intercept = etf_intercept_grp_off
      ) |>
      ## raw rolling, 201
      rolling_etf(
        x = expected_D47,
        y = D47_raw_mean,
        std = etf_stds,
        width = etf_width,
        slope = etf_slope_raw,
        intercept = etf_intercept_raw
      ) |>
      ## offset rolling, 201
      rolling_etf(
        x = expected_D47,
        y = D47_offset_corrected,
        std = etf_stds,
        width = etf_width,
        slope = etf_slope,
        intercept = etf_intercept
      ) |>
      apply_etf(
        intercept = etf_intercept_raw,
        slope = etf_slope_raw,
        raw = D47_raw_mean,
        out = D47_final_roll_no_offset
      ) |>
      apply_etf(
        intercept = etf_intercept,
        slope = etf_slope,
        raw = D47_offset_corrected,
        out = D47_final_roll
      ) |>
      apply_etf(
        intercept = etf_intercept_grp,
        slope = etf_slope_grp,
        raw = D47_raw_mean,
        out = D47_final_no_offset
      ) |>
      apply_etf(
        intercept = etf_intercept_grp_off,
        slope = etf_slope_grp_off,
        raw = D47_offset_corrected,
        out = D47_final
      ) %>%
      temperature_calculation(
        D47 = D47_final,
        slope = .data$temperature_slope,
        intercept = .data$temperature_intercept
      ) %>%
      temperature_calculation(
        D47 = D47_final_no_offset,
        temp = temperature_no_offset,
        slope = .data$temperature_slope,
        intercept = .data$temperature_intercept
      ) |>
      create_reason_for_outlier() |>
      select(-where(is.list)) |> # this might solve hanging?
      order_columns() |>
      arrange(file_datetime),
    format = "qs"
  ),
  
  ## export
  tar_target(
    motu_export,
    tar_excel(motu_temperature, "C:/Archive/OneDrive - Universiteit Utrecht/Archive/motu_all_data_RAW.xlsx"),
    format = "file"
  ),
  
  tar_target(
    motu_out,
    tar_write(motu_temperature,  "C:/Archive/OneDrive - Universiteit Utrecht/Archive/motu_cycle_level_summaries.rds"),
    format = "file"
  ),
  
  tar_target(
    motu_export_csv,
    tar_csv(motu_temperature, "C:/Archive/OneDrive - Universiteit Utrecht/Archive/motu_all_data_RAW.csv"),
    format = "file"
  ),
  
  
  # pacman pipeline ---------------------------------------------------------
  
  # This is our older mass spectrometer. It is a MAT 253 with a Kiel IV, but it
  # had a Kiel III attached earlier with a different software version. The newer
  # measurement files are the ~.did~ files and the older files are the ~.caf~
  # files.
  tar_target(
    pacman_dids_paths_all,
    list_files("Kiel IV data",
               wd = "//ad.geo.uu.nl/GML/Rawdata/Kiel 253") |>
      file_info() |>
      remove_copies() |>
      batch_files(),
    # it now iterates over the directories
    iteration = "list",
    cue = tarchetypes::tar_cue_age(name = pacman_dids_paths_all,
                                   age = our_update_interval)
  ),
  
  tar_target(
    pacman_dids_paths,
    pacman_dids_paths_all, # |>
    # this is to quickly play around with a subset
    ## vctrs::vec_c() |>
    ## vctrs::vec_slice(c(1:3, floor(length(.)/2) + c(-1,0,1), length(.) + c(-2, -1, 0))),
    iteration = "list"
  ),
  
  tar_target(
    pacman_did_files,
    pacman_dids_paths,
    format = "file",
    pattern = map(pacman_dids_paths)
  ),
  
  ## caf files
  
  tar_target(
    pacman_caf_paths_all,
    list_files("clumped/Results", ".caf$",
               wd = "//ad.geo.uu.nl/GML/Rawdata/Kiel 253") |>
      file_info() |>
      remove_copies() |>
      batch_files(),
    # it now iterates over the directories
    iteration = "list",
    cue = tarchetypes::tar_cue_age(name = pacman_caf_paths_all,
                                   age = our_update_interval)
  ),
  
  tar_target(
    pacman_caf_paths,
    pacman_caf_paths_all, # |>
    # this is to quickly play around with a subset
    ## vctrs::vec_c() |>
    ## vctrs::vec_slice(c(1:3, floor(length(.)/2) + c(-1,0,1), length(.) + c(-2, -1, 0))),
    iteration = "list"
  ),
  
  tar_target(
    pacman_caf_files,
    pacman_caf_paths,
    format = "file",
    pattern = map(pacman_caf_paths)
  ),
  
  ## scn files
  #"Background Scans"
  #wd = "//ad.geo.uu.nl/GML/Rawdata/Kiel 253"
  # scn files
  tar_target(
    pacman_scn_paths_all,
    list_files("Scans", ".scn$",
               wd = "//ad.geo.uu.nl/GML/Rawdata/Kiel 253/clumped") |>
      file_info() |>
      remove_copies() |>
      batch_month(),
    iteration = "list",
    cue = tarchetypes::tar_cue_age(name = pacman_scn_paths_all,
                                   age = our_update_interval))
  ),
  
  tar_target(
    pacman_scn_paths,
    pacman_scn_paths_all, ##  |>
    # small subset!
    ## vctrs::vec_c() |>
    ## vctrs::vec_slice(c(1:3, floor(length(.)/2) + c(-1,0,1), length(.) + c(-2, -1, 0))),
    iteration = "list"
  ),
  
  tar_target(
    pacman_scn_files,
    pacman_scn_paths,
    format = "file",
    pattern = map(pacman_scn_paths)
  ),
  
  ## read in as isoreader files
  tar_target(
    pacman_cafs,
    read_di(pacman_caf_files),
    pattern = map(pacman_caf_files),
    iteration = "list",
    format = "qs"
    #cue = tar_cue(command = FALSE)
  ),
  
  tar_target(
    pacman_dids,
    read_di(pacman_did_files),
    pattern = map(pacman_did_files),
    iteration = "list",
    format = "qs"
    #cue = tar_cue(command = FALSE)
  ),
  
  tar_target(
    pacman_scn,
    read_scn(pacman_scn_files),
    pattern = map(pacman_scn_files),
    iteration = "list",
    format = "qs"
    #cue = tar_cue(command = FALSE)
  ),
  
  ## extract raw data
  tar_target(
    pacman_caf_raw,
    iso_get_raw_data(pacman_cafs, include_file_info = Analysis),
    # this now iterates over the folders, so it won't have to re-run this expensive function
    pattern = map(pacman_cafs),
    iteration = "list",
    format = "qs"
  ),
  
  tar_target(
    pacman_raw,
    iso_get_raw_data(pacman_dids, include_file_info = Analysis),
    #|>
    # this now iterates over the folders, so it won't have to re-run this expensive function
    pattern = map(pacman_dids),
    iteration = "list",
    format = "qs"
  ),
  
  tar_target(
    pacman_scn_raw,
    iso_get_raw_data(pacman_scn, include_file_info = c(file_root, file_datetime)),
    pattern = map(pacman_scn),
    iteration = "list",
    format = "qs"
  ),
  
  ## read in metadata
  ## TODO: figure out how to get this from onedrive automatically
  tar_target(
    pacman_did_meta_file,
    "C:/Archive/OneDrive - Universiteit Utrecht/Archive/pacman_did_metadata_parameters.xlsx",
    format = "file"
  ),
  
  tar_target(
    pacman_metadata,
    readxl::read_excel(pacman_did_meta_file, guess_max = 1e5) |>
      meta_fix_types() |> # TODO: switch to parse_col_types?
      tidylog::distinct(Analysis, ## file_id, # there are some with unique file_id's but the same file contents
                        file_size, file_datetime, .keep_all = TRUE),
    format = "fst_tbl"
  ),
  
  tar_target(
    pacman_caf_meta_file,
    "C:/Archive/OneDrive - Universiteit Utrecht/Archive/pacman_caf_metadata_parameters.xlsx",
    format = "file"
  ),
  
  tar_target(
    pacman_caf_metadata,
    readxl::read_excel(pacman_caf_meta_file, guess_max = 1e5) |>
      meta_fix_types() |> # TODO: switch to parse_col_types?
      # hardcoded hack to deal with weight inconsistencies
      ## mutate(`Weight [mg]` = as.character(`Weight [mg]`)) |>
      tidylog::distinct(Analysis, ## file_id, # there are some with unique file_id's but the same file contents
                        file_size, file_datetime, .keep_all = TRUE),
    format = "fst_tbl"
  ),
  
  tar_target(
    pacman_scn_meta_file,
    "C:/Archive/OneDrive - Universiteit Utrecht/Archive/pacman_scn_metadata_parameters.xlsx",
    format = "file"
  ),
  
  tar_target(
    pacman_scn_meta,
    read_xlsx(pacman_scn_meta_file, sheet = "data", guess_max = 1e5)
  ),
  
  ## process scans
  # TODO: import/export pacman_scn_metadata so that I output all parameter columns
  tar_target(
    pacman_scn_fix,
    pacman_scn_raw |>
      nest_by(file_id, file_root, file_datetime) |>
      # this gets some metadata from the raw scan
      file_name_scn() |>
      # this is a way to create the metadata file for the first time:
      ## mutate(min = -500, max = 50000,
      ##        min_start_44 = 9.392386, min_end_44 = 9.395270,
      ##        min_start_45_49 = 9.424277, min_end_45_49 = 9.429723,
      ##        max_start = 9.464633, max_end = 9.468291) |>
      add_scan_info(
        pacman_scn_meta,
        c(
          "min",
          "max",
          "min_start_44",
          "min_end_44",
          "min_start_45_49",
          "min_end_45_49",
          "max_start",
          "max_end",
          "manual_outlier",
          "fix_software",
          "scan_group_overwrite",
          "voltage_overwrite",
          "checked_by",
          "checked_date",
          "checked_comment"
        )
      )  |>
      fix_scan_meta() |>
      unnest(cols = c(data)) |>
      fix_motu_scans(),
    # this hopefully does nothing here!
    pattern = map(pacman_scn_raw),
    iteration = "list",
    format = "qs"
  ),
  
  tar_target(
    pacman_scn_mod,
    pacman_scn_fix |>
      tidy_scans() |>
      flag_scan_ranges() |>
      flag_scan_capped() |>
      calculate_min_max() |>
      calculate_scan_models(),
    ## unnest(data) |>
    pattern = map(pacman_scn_fix),
    iteration = "list",
    format = "qs"
  ),
  
  ## clean up metadata, make file info
  
  # extracted because it's slow and never changes after reading it once
  tar_target(
    pacman_file_info_raw,
    extract_file_info(pacman_dids),
    pattern = map(pacman_dids),
    iteration = "list",
    ## cue = tar_cue(command = FALSE),
    format = "qs"
  ),
  
  tar_target(
    pacman_caf_file_info_raw,
    extract_file_info(pacman_cafs),
    pattern = map(pacman_cafs),
    iteration = "list",
    format = "qs"
  ),
  
  tar_target(
    pacman_file_info,
    pacman_file_info_raw |>
      rowwise() |>
      # this adds all the _overwrite columns and manual_outlier etc.
      fix_metadata(pacman_metadata, irms = "Pacman-KielIV") |>
      # this then applies them to calculate identifier_1 etc.
      overwrite_meta(stdnames = stdnames) |>
      add_scan_group(
        pacman_scn_mod |>
          bind_rows() |>
          tidylog::distinct(scan_datetime, .keep_all = TRUE) |>
          tidylog::filter(!is.na(scan_group), !is.na(lm_47))
      ) |>
      # this adds the parameters that are now in pacman_metadata in stead of parms
      add_parameters(pacman_metadata) |>
      rename(c("outlier_manual" = "manual_outlier")),
    pattern = map(pacman_file_info_raw),
    iteration = "list",
    format = "qs"
  ),
  
  tar_target(
    pacman_caf_file_info,
    pacman_caf_file_info_raw |>
      # I don't fully understand why, but it does the mutate for whole preparations otherwise, resulting in duplicated identifier_1s
      rowwise() |>
      # this adds all the _overwrite columns and manual_outlier etc.
      fix_metadata(pacman_caf_metadata, irms = "Pacman-KielIII") |>
      # this then applies them to calculate identifier_1 etc.
      overwrite_meta(stdnames = stdnames) |>
      add_scan_group(
        pacman_scn_mod |>
          bind_rows() |>
          tidylog::distinct(scan_datetime, .keep_all = TRUE) |>
          tidylog::filter(!is.na(scan_group), !is.na(lm_47))
      ) |>
      # this adds the parameters that are now in pacman_metadata in stead of parms
      add_parameters(pacman_caf_metadata) |>
      rename(c("outlier_manual" = "manual_outlier")),
    pattern = map(pacman_caf_file_info_raw),
    iteration = "list",
    format = "qs"
  ),
  
  
  # this is a subset target so that the raw part only needs to be run when these
  # specific metadata are updated
  tar_target(
    pacman_raw_file_info,
    pacman_file_info |>
      bind_rows() |>
      select(
        Analysis,
        file_id,
        dis_min,
        dis_max,
        dis_fac,
        dis_rel,
        # cycle_filter
        bg_group,
        starts_with("scan_"),
        starts_with("lm_"),
        #starts_with("slope_"),
        bg_fac,
        d13C_PDB_wg,
        d18O_PDBCO2_wg,
        Mineralogy
      ),
    pattern = map(pacman_file_info),
    iteration = "list"
  ),
  
  tar_target(
    pacman_caf_raw_file_info,
    pacman_caf_file_info |>
      bind_rows() |>
      select(
        Analysis,
        file_id,
        dis_min,
        dis_max,
        dis_fac,
        dis_rel,
        # cycle_filter
        bg_group,
        starts_with("scan_"),
        starts_with("lm_"),
        #starts_with("slope_"),
        bg_fac,
        d13C_PDB_wg,
        d18O_PDBCO2_wg,
        Mineralogy
      ),
    pattern = map(pacman_caf_file_info),
    iteration = "list"
  ),
  
  tar_target(
    pacman_badruns,
    pacman_file_info |> 
      bind_rows() |> 
      find_bad_runs()
  ),
  
  # creating pacman_caf_badruns doesn't make sense as the caf files do not have
  # Preparation info
  
  tar_target(
    pacman_meta_update,
    export_metadata(
      data = pacman_file_info |>
        bind_rows(),
      meta = pacman_metadata,
      file = "out/pacman_metadata_update.xlsx"
    ),
    format = "file"
  ),
  
  tar_target(
    pacman_caf_meta_update,
    export_metadata(
      data = pacman_caf_file_info |>
        bind_rows(),
      meta = pacman_caf_metadata,
      file = "out/pacman_caf_metadata_update.xlsx"
    ),
    format = "file"
  ),
  
  tar_target(
    pacman_scn_meta_update,
    export_scan_metadata(
      data = pacman_scn_mod |>
        bind_rows() |>
        unnest(c(data, scan_files)),
      meta = pacman_scn_meta,
      file = "out/pacman_scn_metadata_update.xlsx"
    ),
    format = "file"
  ),
  
  ## pacman caf
  tar_target(
    pacman_caf_raw_deltas,
    pacman_caf_raw |>
      # somehow it's become a character
      # there are some without Analysis
      tidylog::mutate(
        Analysis = ifelse(
          "Analysis" %in% colnames(pacman_caf_raw),
          readr::parse_double(Analysis),
          NA_real_
        )
      ) |>
      # write a wrapper function for this so that the targets are simpler
      # TODO figure out how to loop over two separate lists of both raw and meta info
      add_info(
        pacman_caf_raw_file_info,
        c("dis_min", "dis_max", "dis_fac", "dis_rel")
      ) |>
      clumpedr::find_bad_cycles(
        min = dis_min,
        max = dis_max,
        fac = dis_fac,
        # TODO: get relative_to parms call to work based on dataframe itself
        relative_to = "init"
      ) |>
      filter_duplicated_raw_cycles() |>
      clumpedr::spread_match(masses = 44:49) |>
      add_info(.info = pacman_caf_raw_file_info, cols = c("bg_group", "bg_fac")) |>
      add_background_info(pacman_scn_mod |> bind_rows()) |>
      # TODO: use neighbouring scans before and after sample to get rid of scan noise?
      correct_backgrounds_scn(fac = .data$bg_fac) |>
      # remove the scan models because they take up a lot of memory as list columns
      select(-starts_with("lm_")) |>
      add_info(.info = pacman_caf_raw_file_info,
               c("d13C_PDB_wg", "d18O_PDBCO2_wg")) |>
      clumpedr::abundance_ratios(s44, s45_bg, s46_bg, s47_bg, s48_bg, s49_bg) |>
      clumpedr::abundance_ratios(
        r44,
        r45_bg,
        r46_bg,
        r47_bg,
        r48_bg,
        r49_bg,
        R45_wg,
        R46_wg,
        R47_wg,
        R48_wg,
        R49_wg
      ) |>
      clumpedr::little_deltas() |>
      add_info(pacman_caf_raw_file_info, c("Mineralogy")) |>
      add_R18() |>
      # TODO check if this works for dolomite samples, not sure if vectorized
      clumpedr::bulk_and_clumping_deltas(R18_PDB = .data$R18_PDB) |>
      clumpedr::summarise_outlier(quiet = TRUE),
    # TODO: exclude values mass 54/48/49 < -490
    # TODO: decide whether p49 can be ignored here? I think so because we're doing it at sample level now
    ## add_info(pacman_caf_file_info |> bind_rows(), c("Analysis", "p49_crit")) |>
    ## clumpedr::find_R_flags() |>  # TODO: get rid of R_flags? do they find anything of value?
    pattern = map(pacman_caf_raw, pacman_caf_raw_file_info),
    iteration = "list",
    format = "qs"
  ),
  
  # nesting and summarising still happens within each folder, because this is
  # slow for the big db
  tar_target(
    pacman_caf_nested,
    pacman_caf_raw_deltas |>
      clumpedr::nest_cycle_data(masses = 44:49) |>
      summarize_d13C_d18O_D47(),
    pattern = map(pacman_caf_raw_deltas),
    iteration = "list",
    format = "qs"
  ),
  
  ## sample level summaries
  tar_target(
    pacman_caf_sample_level,
    pacman_caf_nested |>
      bind_rows() |>  # finally the data are rbinded into one big df!
      string_scan_files() |>
      add_remaining_meta(pacman_caf_file_info |> bind_rows()) |>
      # there are 21 measurements for pacman_caf that are so messed up
      # that they cannot be matched by both file_id and Analysis somehow
      tidylog::filter(!is.na(init_low)) |>
      clumpedr::find_init_outliers(
        init_low = init_low,
        init_high = init_high,
        init_diff = init_diff
      ) |>
      summarise_cycle_outliers() |>
      mutate(outlier_param49 = param_49_mean > p49_crit |
               param_49_mean < -p49_crit) |>
      ## summarize_outlier() |>
      # try out conservative outlier selection
      mutate(
        outlier = outlier_noscan | outlier_nodelta |
          (!is.na(outlier_cycles) & outlier_cycles) |
          (!is.na(outlier_init) & outlier_init) |
          (!is.na(outlier_manual) & outlier_manual)
      ) |>
      # get rid of raw cycle data
      ## select(-cycle_data) |>
      select(-where(is.list)) |>
      arrange(file_datetime) |>
      # get rid of duplicated rows
      tidylog::distinct(Analysis, file_id, file_size, s44_init, r44_init, .keep_all = TRUE) |>
      tidylog::filter(!is.na(file_datetime)) |> # this is needed for pacman caf because there are 3 NA rows!
      offset_correction_wrapper(acc = accepted_standard_values),
    format = "fst_tbl"
  ),
  
  tar_target(
    pacman_caf_temperature,
    pacman_caf_sample_level |>
      # there are many ways of calculating the ETF
      ## raw session
      clumpedr:::calculate_etf(
        raw = D47_raw_mean,
        exp = expected_D47,
        session = etf_grp,
        etf = etf,
        etf_coefs = etf_coefs,
        slope = etf_slope_grp,
        intercept = etf_intercept_grp
      ) |>
      ## offset corrected session
      clumpedr:::calculate_etf(
        raw = D47_offset_corrected,
        exp = expected_D47,
        session = etf_grp,
        etf = etf_off,
        etf_coefs = etf_coefs_off,
        slope = etf_slope_grp_off,
        intercept = etf_intercept_grp_off
      ) |>
      ## raw rolling, 201
      rolling_etf(
        x = expected_D47,
        y = D47_raw_mean,
        std = etf_stds,
        width = etf_width,
        slope = etf_slope_raw,
        intercept = etf_intercept_raw
      ) |>
      ## offset rolling, 201
      rolling_etf(
        x = expected_D47,
        y = D47_offset_corrected,
        std = etf_stds,
        width = etf_width,
        slope = etf_slope,
        intercept = etf_intercept
      ) |>
      apply_etf(
        intercept = etf_intercept_raw,
        slope = etf_slope_raw,
        raw = D47_raw_mean,
        out = D47_final_roll_no_offset
      ) |>
      apply_etf(
        intercept = etf_intercept,
        slope = etf_slope,
        raw = D47_offset_corrected,
        out = D47_final_roll
      ) |>
      apply_etf(
        intercept = etf_intercept_grp,
        slope = etf_slope_grp,
        raw = D47_raw_mean,
        out = D47_final_no_offset
      ) |>
      apply_etf(
        intercept = etf_intercept_grp_off,
        slope = etf_slope_grp_off,
        raw = D47_offset_corrected,
        out = D47_final
      ) |>
      temperature_calculation(
        D47 = D47_final,
        slope = .data$temperature_slope,
        intercept = .data$temperature_intercept
      ) |>
      temperature_calculation(
        D47 = D47_final_no_offset,
        temp = temperature_no_offset,
        slope = .data$temperature_slope,
        intercept = .data$temperature_intercept
      ) |>
      create_reason_for_outlier() |>
      select(-where(is.list)) |> # this might solve hanging?
      order_columns() |>
      arrange(file_datetime),
    format = "qs"
  ),
  
  ## pacman
  tar_target(
    pacman_raw_deltas,
    pacman_raw |>
      rowwise() |>
      # somehow it's become a character
      # there are some without Analysis
      tidylog::mutate(
        Analysis = ifelse(
          "Analysis" %in% colnames(pacman_raw),
          readr::parse_double(Analysis),
          NA_real_
        )
      ) |>
      # write a wrapper function for this so that the targets are simpler
      # TODO figure out how to loop over two separate lists of both raw and meta info
      add_info(
        pacman_raw_file_info,
        c("dis_min", "dis_max", "dis_fac", "dis_rel")
      ) |>
      clumpedr::find_bad_cycles(
        min = dis_min,
        max = dis_max,
        fac = dis_fac,
        # TODO: get relative_to parms call to work based on dataframe itself
        relative_to = "init"
      ) |>
      filter_duplicated_raw_cycles() |>
      clumpedr::spread_match(masses = 44:49) |>
      add_info(.info = pacman_raw_file_info, cols = c("bg_group", "bg_fac")) |>
      add_background_info(pacman_scn_mod |> bind_rows()) |>
      # TODO: use neighbouring scans before and after sample to get rid of scan noise?
      correct_backgrounds_scn(fac = .data$bg_fac) |>
      # remove the scan models because they take up a lot of memory as list columns
      select(-starts_with("lm_")) |>
      add_info(.info = pacman_raw_file_info,
               c("d13C_PDB_wg", "d18O_PDBCO2_wg")) |>
      clumpedr::abundance_ratios(s44, s45_bg, s46_bg, s47_bg, s48_bg, s49_bg) |>
      clumpedr::abundance_ratios(
        r44,
        r45_bg,
        r46_bg,
        r47_bg,
        r48_bg,
        r49_bg,
        R45_wg,
        R46_wg,
        R47_wg,
        R48_wg,
        R49_wg
      ) |>
      clumpedr::little_deltas() |>
      add_info(pacman_raw_file_info, c("Mineralogy")) |>
      add_R18() |>
      # TODO check if this works for dolomite samples, not sure if vectorized
      clumpedr::bulk_and_clumping_deltas(R18_PDB = .data$R18_PDB) |>
      clumpedr::summarise_outlier(quiet = TRUE),
    # TODO: exclude values mass 54/48/49 < -490
    # TODO: decide whether p49 can be ignored here? I think so because we're doing it at sample level now
    ## add_info(pacman_file_info |> bind_rows(), c("Analysis", "p49_crit")) |>
    ## clumpedr::find_R_flags() |>  # TODO: get rid of R_flags? do they find anything of value?
    pattern = map(pacman_raw, pacman_raw_file_info),
    iteration = "list",
    format = "qs"
  ),
  
  # nesting and summarising still happens within each folder, because this is
  # slow for the big db
  tar_target(
    pacman_nested,
    pacman_raw_deltas |>
      clumpedr::nest_cycle_data(masses = 44:49) |>
      summarize_d13C_d18O_D47(),
    pattern = map(pacman_raw_deltas),
    iteration = "list",
    format = "qs"
  ),
  
  ## sample level summaries
  tar_target(
    pacman_sample_level,
    pacman_nested |>
      bind_rows() |>  # finally the data are rbinded into one big df!
      string_scan_files() |>
      # this shows that there are quite a few duplicates in pacman did!
      add_remaining_meta(pacman_file_info |> bind_rows()) |>
      clumpedr::find_init_outliers(
        init_low = init_low,
        init_high = init_high,
        init_diff = init_diff
      ) |>
      summarise_cycle_outliers() |>
      mutate(outlier_param49 = param_49_mean > p49_crit |
               param_49_mean < -p49_crit) |>
      ## summarize_outlier() |>
      # try out conservative outlier selection
      mutate(
        outlier = outlier_noscan | outlier_nodelta |
          (!is.na(outlier_cycles) & outlier_cycles) |
          (!is.na(outlier_init) & outlier_init) |
          (!is.na(outlier_manual) & outlier_manual)
      ) |>
      # get rid of raw cycle data
      ## select(-cycle_data) |>
      select(-where(is.list)) |>
      arrange(file_datetime) |>
      # get rid of duplicated rows
      tidylog::distinct(
        Analysis,
        file_id,
        file_datetime,
        file_size,
        s44_init,
        r44_init,
        .keep_all = TRUE
      ) |>
      offset_correction_wrapper(acc = accepted_standard_values),
    format = "fst_tbl"
  ),
  
  tar_target(
    pacman_temperature,
    pacman_sample_level |>
      # there are many ways of calculating the ETF
      ## raw session
      clumpedr:::calculate_etf(
        raw = D47_raw_mean,
        exp = expected_D47,
        session = etf_grp,
        etf = etf,
        etf_coefs = etf_coefs,
        slope = etf_slope_grp,
        intercept = etf_intercept_grp
      ) |>
      ## offset corrected session
      clumpedr:::calculate_etf(
        raw = D47_offset_corrected,
        exp = expected_D47,
        session = etf_grp,
        etf = etf_off,
        etf_coefs = etf_coefs_off,
        slope = etf_slope_grp_off,
        intercept = etf_intercept_grp_off
      ) |>
      ## raw rolling, 201
      rolling_etf(
        x = expected_D47,
        y = D47_raw_mean,
        std = etf_stds,
        width = etf_width,
        slope = etf_slope_raw,
        intercept = etf_intercept_raw
      ) |>
      ## offset rolling, 201
      rolling_etf(
        x = expected_D47,
        y = D47_offset_corrected,
        std = etf_stds,
        width = etf_width,
        slope = etf_slope,
        intercept = etf_intercept
      ) |>
      apply_etf(
        intercept = etf_intercept_raw,
        slope = etf_slope_raw,
        raw = D47_raw_mean,
        out = D47_final_roll_no_offset
      ) |>
      apply_etf(
        intercept = etf_intercept,
        slope = etf_slope,
        raw = D47_offset_corrected,
        out = D47_final_roll
      ) |>
      apply_etf(
        intercept = etf_intercept_grp,
        slope = etf_slope_grp,
        raw = D47_raw_mean,
        out = D47_final_no_offset
      ) |>
      apply_etf(
        intercept = etf_intercept_grp_off,
        slope = etf_slope_grp_off,
        raw = D47_offset_corrected,
        out = D47_final
      ) |>
      temperature_calculation(
        D47 = D47_final,
        slope = .data$temperature_slope,
        intercept = .data$temperature_intercept
      ) |>
      temperature_calculation(
        D47 = D47_final_no_offset,
        temp = temperature_no_offset,
        slope = .data$temperature_slope,
        intercept = .data$temperature_intercept
      ) |>
      create_reason_for_outlier() |>
      select(-where(is.list)) |> # this might solve hanging?
      order_columns() |>
      arrange(file_datetime),
    format = "qs"
  ),
  
  ## export
  tar_target(
    pacman_export,
    tar_excel(pacman_temperature,
              "C:/Archive/OneDrive - Universiteit Utrecht/Archive/pacman_all_data_RAW.xlsx"),
    format = "file"
  ),
  
  tar_target(
    pacman_out,
    tar_write(pacman_temperature,
              "C:/Archive/OneDrive - Universiteit Utrecht/Archive/pacman_cycle_level_summaries.rds"),
    format = "file"
  ),
  
  tar_target(
    pacman_caf_out,
    tar_write(
      pacman_caf_temperature,
      "C:/Archive/OneDrive - Universiteit Utrecht/Archive/pacman_caf_cycle_level_summaries.rds"
    ),
    format = "file"
  ),
  
  tar_target(
    pacman_export_csv,
    tar_csv(pacman_temperature,
            "C:/Archive/OneDrive - Universiteit Utrecht/Archive/pacman_all_data_RAW.csv"),
    format = "file"
  )
)

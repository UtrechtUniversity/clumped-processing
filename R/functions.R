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

# read in the raw files ---------------------------------------------------
list_files <- function(path = "motu/dids",
                       pattern = ".did$",
                       recursive = TRUE,
                       wd = "/home/japhir/Documents/archive") {
  list.files(path = paste(wd, path, sep = "/"),
             pattern = pattern,
             full.names = TRUE,
             recursive = recursive)
}

file_info <- function(files) {
  tibble(
    file_path = files,
    file_name = basename(file_path),
    file_dir = dirname(file_path),
    file_size = file.info(file_path)$size,
    file_datetime = file.info(file_path)$mtime,
    file_year = lubridate::year(file_datetime),
    file_month = lubridate::month(file_datetime))
    ## file_week = lubridate::week(file_datetime))
}

remove_copies <- function(data) {
  tidylog::distinct(data, file_name, file_size, .keep_all = TRUE)
}

# Batch reading in the files so that we have fewer dynamic targets. Do this per
# directory of results.
batch_files <- function(data) {
  tapply(data$file_path,
         INDEX = data$file_dir,  # also possible to batch by directory
         identity, simplify = FALSE) |>
    unname()
}

# The scans are not listed in separate directories, so we batch them by
# year+month.
batch_month <- function(data) {
  tapply(data$file_path,
         INDEX = data$file_year + 1/12 * data$file_month,
         identity, simplify = FALSE) |>
    unname()
}

read_di <- function(data, cache = FALSE, parallel = TRUE, quiet = FALSE) {
  iso_read_dual_inlet(data, cache = cache, parallel = parallel, quiet = quiet)
}

read_scn <- function(data, cache = FALSE, parallel = TRUE, quiet = FALSE) {
  iso_read_scan(data, cache = cache, parallel = parallel, quiet = quiet)
}

# clean up metadata -------------------------------------------------------
meta_fix_types <- function(data) {
  data |>
    # new format with parms included
    readr::type_convert(col_types = cols(Analysis = "i",
                                  file_id = "c",
                                  file_root = "c",
                                  file_subpath = "T",
                                  file_path = "c",
                                  file_datetime = "d",
                                  file_size = "i",
                                  Row = "i",
                                  `Peak Center` = "i",
                                  Background = "i",
                                  Pressadjust = "i",
                                  `Reference Refill` = "i",
                                  Line = "i",
                                  Sample = "i",
                                  `Weight [mg]` = "c",
                                  `Identifier 1` = "c",
                                  `Identifier 2` = "c",
                                  Comment = "c",
                                  Preparation = "c",
                                  Method = "c",
                                  # new columns!
                                  ref_mbar = "d",
                                  ref_pos = "d",
                                  bellow_pos_smp = "d",
                                  init_int = "d",
                                  background = "l",
                                  PC = "i",
                                  VM1_aftr_trfr = "i",
                                  CO2_after_exp = "i",
                                  no_exp = "i",
                                  total_CO2 = "i",
                                  p_gases = "i",
                                  p_no_acid = "i",
                                  extra_drops = "i",
                                  leak_rate = "i",
                                  acid_temperature = "d",
                                  MS_integration_time.s = "i",
                                  timeofday = "d",
                                  d13C_PDB_wg = "d",
                                  d18O_PDBCO2_wg = "d",
                                  # /new columns
                                  s44_init = "d",
                                  r44_init = "d",
                                  # more new parms columns
                                  ## bg_group = "c",
                                  scan_group = "c",
                                  scan_datetime = "T",
                                  scan_files = "c",
                                  scan_n = "i",
                                  bg_fac = "d",
                                  dis_min = "d",
                                  dis_max = "d",
                                  dis_fac = "d",
                                  dis_rel = "c",
                                  init_low = "d",
                                  init_high = "d",
                                  init_diff = "d",
                                  p49_crit = "d",
                                  prop_bad_param49 = "d",
                                  prop_bad_cyc = "d",
                                  sd_D47 = "d",
                                  sd_d13C = "d",
                                  sd_d18O = "d",
                                  off_D47_min = "d",
                                  off_D47_max = "d",
                                  off_D47_grp = "c",
                                  off_D47_width = "i",
                                  off_D47_stds = "c",
                                  off_d13C_min = "d",
                                  off_d13C_max = "d",
                                  off_d13C_grp = "c",
                                  off_d13C_width = "i",
                                  off_d13C_stds = "c",
                                  off_d18O_min = "d",
                                  off_d18O_max = "d",
                                  off_d18O_grp = "c",
                                  off_d18O_width = "i",
                                  off_d18O_stds = "c",
                                  etf_stds = "c",
                                  etf_width = "i",
                                  acid_fractionation_factor = "d",
                                  temperature_slope = "d",
                                  temperature_intercept = "d",
                                  # /parms columns
                                  manual_outlier = "l",
                                  Preparation_overwrite = "d",
                                  `Identifier 1_overwrite` = "c",
                                  `Identifier 2_overwrite` = "c",
                                  `Weight [mg]_overwrite` = "d",
                                  Comment_overwrite = "c",
                                  scan_group_overwrite = "c",
                                  Mineralogy = "c",
                                  checked_by = "c",
                                  checked_date = "T",
                                  checked_comment = "c")) #|>
     ## mutate(Preparation = as.double(Preparation))
}

filter_info_duplicates <- function(data) {
  data |>
    tidylog::distinct(file_id, file_datetime, file_size, .keep_all=TRUE)
}

add_timeofday <- function(data) {
  message("Info: adding timeofday")
  data |>
    mutate(timeofday = lubridate::hour(file_datetime) +
             lubridate::minute(file_datetime) / 60 +
             lubridate::second(file_datetime) / 60 / 60)
}

# This compares the preparation/run number inside the file with the one in the
# filename/filepath.
find_bad_runs <- function(data) {
  out <- data |>
    file_name_prep() |>
    tidylog::filter(Preparation != file_id_prep) |>
    select(file_id, Preparation, file_id_prep) |>
    tidylog::distinct(Preparation, file_id_prep, .keep_all = TRUE)
}

file_name_prep <- function(data) {
  data |>
    mutate(file_id_prep = str_extract(file_id, "_\\d{1,3}_?(restart_)?B?") |>
             str_replace_all("_", "") |> str_replace_all("restart", "") |>
             str_replace_all("B", "") |> parse_integer())
}

parse_col_types <- function(data) {
  data |>
    readr::type_convert(col_types = cols(file_id = "c",
                                  file_root = "c",
                                  file_path = "c",
                                  file_subpath = "c",
                                  file_datetime = "T",
                                  file_size = "i",
                                  Row = "i",
                                  `Peak Center` = "l",
                                  Background = "l",
                                  Pressadjust = "l",
                                  `Reference Refill` = "l",
                                  Line = "i",
                                  Sample = "i",
                                  `Weight [mg]` = "c",
                                  `Identifier 1` = "c",
                                  `Identifier 2` = "c",
                                  Analysis = "i",
                                  Comment = "c",
                                  Preparation = "c",
                                  Method = "c",
                                  measurement_info = "?",
                                  MS_integration_time.s = "d"))
}

# Here is an example of the info that we're trying to tease apart
# - Acid: 70.0 [°C]
# - LeakRate [µBar/Min]:  171
# - 0 xtra drops
# - P no Acid :    3
# - P gases:   27
# - Total CO2 :  550
# - # Exp.:  0
#   - CO2 after Exp.:  550
# - VM1 aftr Trfr.:    0
# - PC [62040]
# - Background: BGD 2018/Jan/23 03:15 -  (Administrator)
# - Init int: 18050.65
# - Bellow Pos: 100%
# - RefI: mBar r 67.1  pos r 33.7
split_meas_info <- function(data) {
    if (!"measurement_info" %in% colnames(data)) {
      warning("Column `measurement_info` not found in data.")
      return(data)
    }

    data |>
        tidyr::extract(measurement_info,
                into = "acid_temperature",
                regex = "Acid: *(-?\\d+\\.?\\d*) *\\[?°?C?\\]?",
                remove = FALSE,
                convert = TRUE) |>
        tidyr::extract(measurement_info,
                into = "leak_rate",
                regex =   "LeakRate *\\[µBar/Min\\]: *(-?\\d+\\.?\\d*)",
                remove = FALSE,
                convert = TRUE) |>
        tidyr::extract(measurement_info,
                into = "extra_drops",
                regex = "(\\d+) *xtra *drops",
                remove = FALSE,
                convert = TRUE) |>
        tidyr::extract(measurement_info,
                into = "p_no_acid",
                regex = "P no Acid : *(-?\\d+\\.?\\d*)",
                remove = FALSE,
                convert = TRUE) |>
        tidyr::extract(measurement_info,
                into = "p_gases",
                regex = "P gases: *(-?\\d+\\.?\\d*)",
                remove = FALSE,
                convert = TRUE) |>
        tidyr::extract(measurement_info,
                into = "total_CO2",
                regex = "Total CO2 *: *(-?\\d+\\.?\\d*)",
                remove = FALSE,
                convert = TRUE) |>
        tidyr::extract(measurement_info,
                into = "no_exp",
                regex = "# Exp\\.: *(-?\\d+\\.?\\d*)?",
                remove = FALSE,
                convert = TRUE) |>
        tidyr::extract(measurement_info,
                into = "CO2_after_exp",
                regex = "CO2 after Exp\\.: *(-?\\d+\\.?\\d*)",
                remove = FALSE,
                convert = TRUE) |>
        tidyr::extract(measurement_info,
                into = "VM1_aftr_trfr",
                regex = "VM1 *aftr *Trfr\\.: *(-?\\d+\\.?\\d*)",
                remove = FALSE,
                convert = TRUE) |>
        tidyr::extract(measurement_info,
                into = "PC",
                regex = "PC \\[(-?\\d+\\.?\\d*)\\]",
                remove = FALSE,
                convert = TRUE) |>
        tidyr::extract(measurement_info,
                into = "background",
                regex = "Background: (.*)\n",
                remove = FALSE) |>
        tidyr::extract(measurement_info,
                into = "init_int",
                regex =  "Init int: *(-?\\d+\\.?\\d*)",
                remove = FALSE,
                convert = TRUE) |>
        tidyr::extract(measurement_info,
                into = "bellow_pos_smp",
                regex = "Bellow Pos: *(-?\\d+\\.?\\d*)%",
                remove = FALSE,
                convert = TRUE) |>
        tidyr::extract(measurement_info,
                into = c("ref_mbar", "ref_pos"),
                regex = "RefI: *mBar *r *(-?\\d+\\.?\\d*) *pos *r *(-?\\d+\\.?\\d*)",
                remove = FALSE,
                convert = TRUE)
  }

#' this adds the initial intensities from dids to the metadata
add_inits <- function(data, dids) {
  inits <- dids |>
    iso_get_raw_data(select = c(cycle, type, v44.mV),
                     include_file_info = Analysis)

  # if it's a normal case with file info
  if (!is.null(nrow(inits))) {
    if (nrow(inits) > 0L) {
      inits <- inits |>
        get_inits()
    }
  } else {
    inits <- tibble(file_id = character(), 
                    Analysis = integer(), 
                    s44_init = double(), 
                    r44_init = double())
  }

  left_join(x = data, y = inits, by = c("Analysis", "file_id"))
}

fix_metadata <- function(data, meta, irms = "MotU-KielIV") {
  if (nrow(data) == 0L) {
    return(tibble(file_id = character()))
  }

  out <- data |>
    # add the metadata overwrite columns!
    tidylog::left_join(
               meta |>
               select(.data$Analysis, 
                      .data$file_id, #.data$`Identifier 1`,
                      ends_with("_overwrite"), 
                      .data$manual_outlier, 
                      .data$Mineralogy,
                      starts_with("checked_")), 
               by = c("Analysis", "file_id"))

  # make sure that weight exists and is a double
  if ("Weight [mg]" %in% colnames(out)) {
    out <- out |>
      # in case the weight is not a double, try to parse it automatically
      mutate(weight_double = parse_double(`Weight [mg]`)) |>
      tidylog::mutate(`Weight [mg]` = ifelse(!is.na(weight_double),
                                             weight_double,
                                             # parsing the weight simply failed, trying to extract it from the string
                                             str_extract(`Weight [mg]`, "\\d+.?\\d*") |> parse_double())) |>
      select(-weight_double)
  } else {
    out <- out |> mutate(`Weight [mg]` = NA_real_)
  }

  # this makes sure that Preparation exists and is an integer!
  if ("Preparation" %in% colnames(out)) {
    out <- out |>
      # we don't do anything with preparation_overwrite yet, that's for overwrite_meta
      mutate(Preparation_integer = parse_integer(Preparation)) |>
      # convert preparation to integer, if this didn't work first extract the
      # first number from the text then parse it.
      tidylog::mutate(Preparation = ifelse(!is.na(Preparation_integer),
                                           Preparation_integer,
                                           str_extract(Preparation, "\\d+") |> parse_integer())) |>
      select(-Preparation_integer)  # get rid of temporary column
  } else { # there's no preparation column
    out <- out |>
      mutate(Preparation = NA_integer_)
  }

  out |>
    # get the Preparation number from the directory name, if possible
    # NOTE: 2021-10-11 commented out because hopefully the metadata file has all this info now!
    ## tidylog::mutate(Preparation_overwrite =
    ##                             # Pacman caf naming convention (if adhered to) is YYMMDD_people (so we'll use the date)
    ##                   case_when(irms == "Pacman-KielIII" & is.na(Preparation_overwrite) ~
    ##                               str_extract(file_root, "cafs/\\d{6}") |>
    ##                               str_extract("\\d{6}") |>
    ##                               parse_integer(),
    ##                             # Pacman did naming convention (if adhered to) is _YYMMDD_prep number
    ##                             irms == "Pacman-KielIV" & is.na(Preparation_overwrite) ~
    ##                               str_extract(file_root, "\\d{6}_\\d+$") |>
    ##                               str_extract("\\d+$") |>
    ##                               parse_integer(),
    ##                             irms == "MotU-KielIV" & !is.na(Preparation_overwrite) ~
    ##                               Preparation_overwrite |> as.integer(),
    ##                             TRUE ~ NA_integer_)) |>
    mutate(masspec = irms)
}

add_parameters <- function(data, meta) {
  cd <- colnames(data)
  cm <- colnames(meta)
  cn <- cm[!cm %in% cd]

  data |>
    tidylog::left_join(
               meta |>
               select(.data$Analysis, 
                      .data$file_id, #.data$`Identifier 1`,
                      one_of(cn)),
               by = c("Analysis", "file_id"))
}

overwrite_meta <- function(meta, masspec = "MotU-KielIV", stdnames) {
  if (nrow(meta) == 0L) {
    return(tibble(file_id = character()))
  }

  desired_cols <- c("Preparation", "Identifier 1", "Identifier 2", "Weight [mg]", "Comment")
  cols_exist <- desired_cols %in% colnames(meta)
  if (!all(cols_exist)) {
    warning(glue::glue("Colname(s) '{glue::glue_collapse(desired_cols[!cols_exist], sep = ' ', width = 30L, last = ' and ')}' not found in meta"))
  }

  meta |>
    tidylog::mutate(
               preparation = ifelse("Preparation" %in% colnames(meta) &&
                                    is.na(.data$Preparation_overwrite),
                                    .data$Preparation,
                                    .data$Preparation_overwrite),
               identifier_1 = ifelse("Identifier 1" %in% colnames(meta) &&
                                     is.na(.data$`Identifier 1_overwrite`),
                                     .data$`Identifier 1`,
                                     .data$`Identifier 1_overwrite`),
               identifier_2 = ifelse("Identifier 2" %in% colnames(meta) &&
                                     is.na(.data$`Identifier 2_overwrite`),
                                     .data$`Identifier 2`, .data$`Identifier 2_overwrite`),
               weight = ifelse("Weight [mg]" %in% colnames(meta) &&
                               is.na(.data$`Weight [mg]_overwrite`),
                               .data$`Weight [mg]`,
                               .data$`Weight [mg]_overwrite` |> as.double()),
               comment = ifelse("Comment" %in% colnames(meta) &&
                                is.na(.data$Comment_overwrite),
                                .data$Comment, .data$Comment_overwrite),
               masspec = .data$masspec,
               ## scan_group = ifelse(is.na(scan_group_overwrite), scan_group, scan_group_overwrite),
               broadid = ifelse(.data$identifier_1 %in% stdnames, identifier_1, "other"))
}

filter_raw_duplicates <- function(data) {
  dups <- data |>
    filter(cycle==0, type=="standard") |>
    tidylog::distinct(Analysis, v44.mV, .keep_all = TRUE) # message tells us the number of dups

  data |>
    filter(file_id %in% dups$file_id & Analysis %in% dups$Analysis)
}

export_metadata <- function(data, meta, file) {
       # NOTE: this only filters out the last ones!
     # every now and then if you have issues, check if all
     # analyses are in the metadata file 
     # (e.g. see my 2022-11-04_pacman_extra_manual)
     # tidylog::filter(Analysis > max(meta$Analysis, na.rm = TRUE)) |>
     # NOTE: experimental fix: look for both Analysis and file_id 
     # (because pacman doesn't have unique file_id's)
     # and export everything that isn't in the metadata file yet
     #tidylog::filter(!(Analysis %in% meta$Analysis & 
     #                   file_id %in% meta$file_id)) |>
  if ("Analysis" %in% colnames(data)) {
    data <- data |>
      tidylog::filter(!Analysis %in% meta$Analysis)
  } else {
    warning("Column `Analysis` not found in data.", call. = FALSE)
  }
  
  if ("file_datetime" %in% colnames(data)) {
    data <- data |> 
      tidylog::filter(!file_datetime %in% meta$file_datetime)
  } else {
    warning("Column `file_datetime` not found in data.", call. = FALSE)
  }
  
  if ("file_id" %in% colnames(data)) {
    data <- data |>
      tidylog::filter(!file_id %in% meta$file_id)
  } else {
    warning("Column `file_id` not found in data.", call. = FALSE)
  }
  
  data |>
     rename(c("manual_outlier" = "outlier_manual")) |>
     tidylog::select(all_of(c("Analysis",
                              "file_id",
                              "file_root",
                              "file_subpath",
                              "file_path",
                              "file_datetime",
                              "file_size",
                              "Row",
                              "Peak Center",
                              "Background",
                              "Pressadjust",
                              "Reference Refill",
                              "Line",
                              "Sample",
                              "Weight [mg]",
                              "Identifier 1",
                              "Identifier 2",
                              "Comment",
                              "Preparation",
                              "Method",
                              # new columns!
                              "ref_mbar",
                              "ref_pos",
                              "bellow_pos_smp",
                              "init_int",
                              "background",
                              "PC",
                              "VM1_aftr_trfr",
                              "CO2_after_exp",
                              "no_exp",
                              "total_CO2",
                              "p_gases",
                              "p_no_acid",
                              "extra_drops",
                              "leak_rate",
                              "acid_temperature",
                              "MS_integration_time.s",
                              "timeofday",
                              "d13C_PDB_wg",
                              "d18O_PDBCO2_wg",
                              # /new columns
                              "s44_init",
                              "r44_init",
                              # more new parms columns
                              ## "bg_group",
                              "scan_group",
                              "scan_datetime",
                              "scan_files",
                              "scan_n",
                              "bg_fac",
                              "dis_min", "dis_max", "dis_fac", "dis_rel",
                              "init_low", "init_high", "init_diff",
                              "p49_crit",
                              "prop_bad_param49",
                              "prop_bad_cyc",
                              "sd_D47", "sd_d13C", "sd_d18O",
                              "off_D47_min", "off_D47_max", "off_D47_grp", "off_D47_width", "off_D47_stds",
                              "off_d13C_min", "off_d13C_max", "off_d13C_grp", "off_d13C_width", "off_d13C_stds",
                              "off_d18O_min", "off_d18O_max", "off_d18O_grp", "off_d18O_width", "off_d18O_stds",
                              "etf_stds", "etf_width",
                              "acid_fractionation_factor",
                              "temperature_slope", "temperature_intercept",
                              # /parms columns
                              "manual_outlier",
                              "Preparation_overwrite",
                              "Identifier 1_overwrite",
                              "Identifier 2_overwrite",
                              "Weight [mg]_overwrite",
                              "Comment_overwrite",
                              "scan_group_overwrite",
                              "Mineralogy",
                              "checked_by",
                              "checked_date",
                              "checked_comment"))) |>
     writexl::write_xlsx(file)
   file
}

extract_file_info <- function(did) {
  did |>
    iso_get_file_info() |>
    filter_info_duplicates() |>
    add_inits(did) |>
    parse_col_types() |>
    split_meas_info() |>
    select(-one_of("measurement_info")) |> # this is a list
    add_timeofday() |>
    clumpedr::append_ref_deltas(.did = did)
}

# function only used to create first set of metadata files
create_metadata <- function(meta, file) {
   meta |>
     rename(c("manual_outlier" = "outlier_manual")) |>
     arrange(file_datetime) |>
     tidylog::select(one_of(c("Analysis",
                              "file_id",
                              "file_root",
                              "file_subpath",
                              "file_path",
                              "file_datetime",
                              "file_size",
                              "Row",
                              "Peak Center",
                              "Background",
                              "Pressadjust",
                              "Reference Refill",
                              "Line",
                              "Sample",
                              "Weight [mg]",
                              "Identifier 1",
                              "Identifier 2",
                              "Comment",
                              "Preparation",
                              "Method",
                              # new columns!
                              "ref_mbar",
                              "ref_pos",
                              "bellow_pos_smp",
                              "init_int",
                              "background",
                              "PC",
                              "VM1_aftr_trfr",
                              "CO2_after_exp",
                              "no_exp",
                              "total_CO2",
                              "p_gases",
                              "p_no_acid",
                              "extra_drops",
                              "leak_rate",
                              "acid_temperature",
                              "MS_integration_time.s",
                              "timeofday",
                              "d13C_PDB_wg",
                              "d18O_PDBCO2_wg",
                              # /new columns
                              "s44_init",
                              "r44_init",
                              # more new parms columns
                              ## "bg_group",
                              "scan_group",
                              "scan_datetime",
                              "scan_files",
                              "scan_n",
                              "bg_fac",
                              "dis_min", "dis_max", "dis_fac", "dis_rel",
                              "init_low", "init_high", "init_diff",
                              "p49_crit",
                              "prop_bad_param49",
                              "prop_bad_cyc",
                              "sd_D47", "sd_d13C", "sd_d18O",
                              "off_D47_min", "off_D47_max", "off_D47_grp", "off_D47_width", "off_D47_stds",
                              "off_d13C_min", "off_d13C_max", "off_d13C_grp", "off_d13C_width", "off_d13C_stds",
                              "off_d18O_min", "off_d18O_max", "off_d18O_grp", "off_d18O_width", "off_d18O_stds",
                              "etf_stds", "etf_width",
                              "acid_fractionation_factor",
                              "temperature_slope", "temperature_intercept",
                              # /parms columns
                              "manual_outlier",
                              "Preparation_overwrite",
                              "Identifier 1_overwrite",
                              "Identifier 2_overwrite",
                              "Weight [mg]_overwrite",
                              "Comment_overwrite",
                              "scan_group_overwrite",
                              "Mineralogy",
                              "checked_by",
                              "checked_date",
                              "checked_comment"))) |>
     writexl::write_xlsx(file)
   file
}


# process background scans ------------------------------------------------
file_name_scn <- function(data) {
  if (nrow(data) == 0L) {
    return(tibble(file_id = character()))
  }

  data |>
    # we're searching for numbers/characters, then an underscore. Mostly we use
    # YYMMDD_#V.scn but sometimes something else
    tidylog::mutate(scan_group = str_extract(file_id, "^[\\dA-z]+?_") |>
                      # get rid of the underscore
                      str_replace_all("_", "") |>
                      # another format for 190215 :S
                      str_replace("BG\\d{1,2}V", ""),
                    # we look for the voltage in the filename, must be format NNV or NN.NNV
                    voltage = str_extract(file_id, "\\d+\\.?\\d*V") |>
                      str_replace("V", "") |>
                      parse_double()) |>
    group_by(scan_group) |>
    tidylog::mutate(scan_datetime = first(file_datetime)) |>
    group_by(file_id) |>
    tidylog::mutate(voltage_max = purrr::possibly(map_dbl, NA_real_)(
      data,
      ~ max(.$v44.mV, na.rm = TRUE))) |>
    ungroup(file_id)
}

fix_scan_meta <- function(data) {
  if (nrow(data) == 0L) {
    return(tibble(file_id = character()))
  }

  data |>
    tidylog::mutate(scan_group = ifelse(is.na(scan_group_overwrite),
                                        scan_group,
                                        scan_group_overwrite),
                    voltage = ifelse(is.na(voltage_overwrite),
                                     voltage,
                                     voltage_overwrite),
                    fix_software = ifelse(is.na(fix_software), FALSE, fix_software),
                    outlier_scan_manual = ifelse(is.na(manual_outlier), FALSE, manual_outlier)) |>
    select(-manual_outlier)
}

# We had a mistake in the software setting for some time. Here we undo that
# correction prior to analysis, based on the logical column ~fix_software~ in
# the metadata.
fix_motu_scans <- function(data) {
  if (nrow(data) == 0L) {
    return(tibble(file_id = character()))
  }

  if (!all(c("v47.mV", "v54.mV", "fix_software") %in% colnames(data))) {
    warning("Column names v47.mV, v54.mV and fix_software not found")
    return(data)
  }
  if (sum(data |> distinct(file_id, .keep_all = TRUE) |> pull(fix_software) > 0)) {
    glue::glue("Info: fixing software settings for {sum(data |> distinct(file_id, .keep_all = TRUE) |> pull(fix_software) > 0)} scans.") |>
      message()
  }
  data |>
    tidylog::mutate(v47.mV = ifelse(fix_software, v47.mV - v54.mV, v47.mV))
}

# Tidying is reshaping into long format https://r4ds.had.co.nz/tidy-data.html.
tidy_scans <- function(data) {
  if (!all(c("v44.mV", "v47.mV") %in% colnames(data)) | nrow(data) == 0) {
    return(tibble(file_id = character()))
  }

  data |>
    # there are a bunch of weird columns in Pacman scans that I get rid of here
    tidylog::select(-any_of(c("v17.6.mV", "v18.mV", "v18.4.mV", "v2.mV", "v3.mV")),
                    -matches("v\\d+\\.\\d+\\.mV"),
                    -matches("vC\\d+\\.mV")) |>
    tidylog::pivot_longer(cols = matches("v\\d+\\.mV"), names_pattern = "v(\\d+).mV") |>
    tidylog::mutate(name = parse_integer(name)) |>
    tidylog::rename("mass" = "name", "intensity" = "value")
}

# This creates logical columns to indicate whether a part of a scan should be
# used to calculate the minimum or maximum intensities. It does so based on the
# metadata columns. this one now uses columns!
flag_scan_ranges <- function(data) {
  if (nrow(data) == 0L) {
    return(tibble(file_id = character()))
  }

  if (! all(c("min", "max", "min_start_44", "min_end_44", "min_start_45_49", "min_end_45_49", "max_start", "max_end") %in% colnames(data))) {
    warning("Scan parameters not found, emptying this target!")
    return(tibble(file_id = character()))
  }

  data |>
    ## tidylog::filter(!outlier_scan_manual) |> # get rid of manually labelled failed scans
    tidylog::filter(intensity >= min | is.na(min)) |>
    tidylog::filter(intensity <= max | is.na(max)) |>
    tidylog::mutate(min_sub = ifelse(mass == 44,
                          x > min_start_44 & x < min_end_44,
                                 x > min_start_45_49 & x < min_end_45_49)) |>
    tidylog::mutate(max_sub = x > max_start & x < max_end)
}

# Some scans have values in the minimum range that are less than the sensor can
# actually record. We need to exclude those, so I mark them as outliers here.
# The capped minimum value differs per mass, so I've put the actual capped
# values in here.
flag_scan_capped <- function(data,
                             m44 = -499,
                             m45 = -499,
                             m46 = -499,
                             m47 = -499.0608,
                             m48 = -499.5371,
                             m49 = -498.8829,
                             m54 = -499.6343) {
  if (nrow(data) < 1) {
    return(tibble(file_id = character()))
  }

  crit <- tibble(mass = c(44, 45:49, 54), cap = c(m44, m45, m46, m47, m48, m49, m54))

  minrange <- data |>
    filter(min_sub) |>
    left_join(crit, by = "mass") |>
    group_by(file_id, mass) |>
    mutate(outlier_scan_minimumcap = any(intensity <= cap)) |> # low in the minimum range?
    ungroup(file_id, mass) |>
    distinct(file_id, mass, outlier_scan_minimumcap)

  data |>
    left_join(minrange, by = c("file_id", "mass"))
}

# This calculates the average minimum and maximum values in the flagged ranges.
calculate_min_max <- function(data) {
  if (nrow(data) == 0L) {
    return(tibble(scan_group = character())) # this one doesn't have file_id anymore!
  }

  # this makes sure I only add real metadata, not the min/max/model output
  meta <- data |>
    distinct(file_id,
             file_root,
             file_datetime,
             scan_datetime,
             voltage,
             voltage_max,
             scan_group, min, max,
             min_start_44,
             min_end_44,
             min_start_45_49,
             min_end_45_49,
             max_start,
             max_end,
             outlier_scan_manual,
             fix_software,
             scan_group_overwrite,
             voltage_overwrite,
             checked_by,
             checked_date,
             checked_comment)

  max_intensity <- data |>
    filter(max_sub | is.na(max_sub)) |>
    group_by(file_id, file_root, file_datetime, voltage, voltage_max, mass, scan_group, scan_datetime) |>
    summarise(measure = "max", value = mean(intensity))

  min_intensity <- data |>
    filter(min_sub | is.na(min_sub))  |>
    tidylog::filter(is.na(outlier_scan_minimumcap) | !outlier_scan_minimumcap) |>
    group_by(file_id, file_root, file_datetime, voltage, voltage_max, mass, scan_group, scan_datetime) |>
    summarise(measure = "min", value = mean(intensity))

  # SOME: how to make pivot_scans not remove all the stuff from before?
  bind_rows(min_intensity, max_intensity) |>
    pivot_scans()  |>
    left_join(meta,
              by = c("file_id",
                     "file_root",
                     "file_datetime",
                     "scan_datetime",
                     "voltage",
                     "voltage_max",
                     "scan_group"))
}

pivot_scans <- function(data) {
  data |>
    ungroup() |>
    tidylog::pivot_wider(names_from = c(measure, mass),
                         values_from = value)
}

# This fits linear models between the minima for the different masses and the
# maximum of mass 44.
calculate_scan_models <- function(data) {
  if (nrow(data) == 0L) {
    return(tibble(scan_group = character()))
  }

  data |>
    group_by(scan_group) |>
    tidyr::nest(data = c(starts_with("file_"), starts_with("voltage"),
                  starts_with("min_4"),
                  starts_with("min_54"),
                  starts_with("max_4"),
                  starts_with("max_54"),
                  # the min max min_start min_end max_start max_end columns
                  # SHOULD all be identical within a scan_group
                  # so we do not need to nest by those. But we do it anyway!
                  min, max, min_start_44, min_end_44, min_start_45_49, min_end_45_49, max_start, max_end,
                  fix_software, scan_group_overwrite,
                  outlier_scan_manual, checked_by, checked_date, checked_comment)) |>
    tidylog::mutate(scan_datetime = map_dbl(data, ~ min(.x$file_datetime)) |>
                      as.POSIXct(origin = "1970-01-01 00:00.00"),
                    scan_files = map(data, ~ paste(.x$file_id)),
                    scan_n = map_dbl(data, ~ nrow(.x)), ## 45 is not linear, but very minor
                    # first fit the mass 44 model to scale everything to 0 to max
                    ## lm_44 = map(data, purrr::possibly(~ lm(min_44 ~ max_44 - 1, data = .x), otherwise = em())),
                    # TODO: look into whether fitting a line through the origin works better? probably not, e.g. 45 behaves a bit non-linearly
                    ## max_44 = predict(lm_44, newdata = max_44),
                    lm_45 = map(data, purrr::possibly(~ lm(formula = min_45 ~ poly(max_44, 3, raw = TRUE) - 1, data = .x |> filter(!outlier_scan_manual)), otherwise = em())),
                    lm_46 = map(data, purrr::possibly(~ lm(formula = min_46 ~ poly(max_44, 3, raw = TRUE) - 1, data = .x |> filter(!outlier_scan_manual)), otherwise = em())),
                    lm_47 = map(data, purrr::possibly(~ lm(formula = min_47 ~ poly(max_44, 3, raw = TRUE) - 1, data = .x |> filter(!outlier_scan_manual)), otherwise = em())),
                    lm_48 = map(data, purrr::possibly(~ lm(formula = min_48 ~ poly(max_44, 3, raw = TRUE) - 1, data = .x |> filter(!outlier_scan_manual)), otherwise = em())),
                    lm_49 = map(data, purrr::possibly(~ lm(formula = min_49 ~ poly(max_44, 3, raw = TRUE) - 1, data = .x |> filter(!outlier_scan_manual)), otherwise = em())),
                    ## coef_44 = map(lm_44, "coefficients"), #otherwise = NA),
                    ## coef_45 = map(lm_45, "coefficients"), #otherwise = NA),
                    ## coef_46 = map(lm_46, "coefficients"),
                    ## coef_47 = map(lm_47, "coefficients"),
                    ## coef_48 = map(lm_48, "coefficients"),
                    ## coef_49 = map(lm_49, "coefficients"),
                    ## ## intercept_44 = map_dbl(coef_44, 1),
                    ## intercept_45 = map_dbl(coef_45, 1),
                    ## intercept_46 = map_dbl(coef_46, 1),
                    ## intercept_47 = map_dbl(coef_47, 1),
                    ## intercept_48 = map_dbl(coef_48, 1),
                    ## intercept_49 = map_dbl(coef_49, 1),
                    ## ## slope_44 = map_dbl(coef_44, 2),
                    ## slope_45 = map_dbl(coef_45, 2),
                    ## slope_46 = map_dbl(coef_46, 2),
                    ## slope_47 = map_dbl(coef_47, 2),
                    ## slope_48 = map_dbl(coef_48, 2),
                    ## slope_49 = map_dbl(coef_49, 2)
                    ) |>
  ## tidylog::select(-starts_with("lm"), -starts_with("coef")) |>
  arrange(scan_datetime) |>
  tidylog::ungroup(scan_group) |>
  tidylog::mutate(scan_duration = c(lubridate::int_diff(scan_datetime), NA_real_)) |>
  tidylog::mutate(bg_group = cut_scan_groups(scan_datetime, scan_datetime)) |>
  tidylog::filter(!is.na(bg_group))
}

# If the model fails, we return an empty model so we can still call ~coef~ on it
# without problems.
em <- function() {
  out  <- list()
  class(out) <- "lm"
  out$coefficients <- c("(Intercept)" = NA, "max_44" = NA)
  out
}

cut_scan_groups <- function(file, scan) {
  cut(file,
      # we need to make sure oldest and newest scans are also assigned a category
      c(parse_datetime("1990-02-13 12:00:00"), # my birthday!
        unique(scan),
        lubridate::now(tzone = "UTC"))) |>
    as.character()
}

add_scan_group <- function(info, bg) {
  if (nrow(bg) == 0) {
    warning("Could not match background, it's empty")
    return(info)
  }

  info |>
    ## tidylog::select(all_of(c("file_id", "file_datetime"))) |>
    tidylog::mutate(bg_group = cut_scan_groups(file_datetime, bg$scan_datetime)) #|>
    ## tidylog::select(-file_datetime) |>
    ## tidylog::left_join(bg,
    ##                    distinct(scan_datetime, bg_group),
    ##                    # the background scans need to be cut up exactly the same as the files
    ##                    ## mutate(bg_group = cut_scan_groups(scan_datetime, scan_datetime)),
    ##                    by = "bg_group")
}

add_background_info <- function(data, bg) {
  message("Info: adding background models")
  if (nrow(data) == 0L) {
    return(tibble(file_info = character()))
  }

  data |>
    tidylog::left_join(bg |>
                       select(bg_group, #file_id,
                              starts_with("scan_"),
                              ## starts_with("intercept_"),
                              starts_with("lm_")#,
                              ## starts_with("slope_"),
                              ## bg_fac
                              ),
                       by = "bg_group") #"file_id"
}

# Apply the background corrections to the raw measurement intensities at the cycle level.
correct_backgrounds_scn <- function(data, fac) {  #  = 0.91, masses = c(44:49, 54)
  message("Info: correcting backgrounds using scan models")
  if (nrow(data) == 0L) {
    return(tibble(file_info = character()))
  }

  ## if (!all(c(paste0("slope_", 45:49),
  ##            paste0("intercept_", 45:49)) %in% colnames(data))) {
  ##   warning("Columns needed for background scans not found!")
  ##   data <- data |>
  ##     mutate(slope_45 = NA_real_,
  ##            slope_46 = NA_real_,
  ##            slope_47 = NA_real_,
  ##            slope_48 = NA_real_,
  ##            slope_49 = NA_real_,
  ##            intercept_45 = NA_real_,
  ##            intercept_46 = NA_real_,
  ##            intercept_47 = NA_real_,
  ##            intercept_48 = NA_real_,
  ##            intercept_49 = NA_real_)
  ## }

  out <- data |>
    mutate(s44_bg45 = map2_dbl(s44, lm_45, ~ purrr::possibly(predict, NA_real_)(.y, newdata = list(max_44 = .x))),
           s44_bg46 = map2_dbl(s44, lm_46, ~ purrr::possibly(predict, NA_real_)(.y, newdata = list(max_44 = .x))),
           s44_bg47 = map2_dbl(s44, lm_47, ~ purrr::possibly(predict, NA_real_)(.y, newdata = list(max_44 = .x))),
           s44_bg48 = map2_dbl(s44, lm_48, ~ purrr::possibly(predict, NA_real_)(.y, newdata = list(max_44 = .x))),
           s44_bg49 = map2_dbl(s44, lm_49, ~ purrr::possibly(predict, NA_real_)(.y, newdata = list(max_44 = .x)))) |>
    mutate(r44_bg45 = map2_dbl(r44, lm_45, ~ purrr::possibly(predict, NA_real_)(.y, newdata = list(max_44 = .x))),
           r44_bg46 = map2_dbl(r44, lm_46, ~ purrr::possibly(predict, NA_real_)(.y, newdata = list(max_44 = .x))),
           r44_bg47 = map2_dbl(r44, lm_47, ~ purrr::possibly(predict, NA_real_)(.y, newdata = list(max_44 = .x))),
           r44_bg48 = map2_dbl(r44, lm_48, ~ purrr::possibly(predict, NA_real_)(.y, newdata = list(max_44 = .x))),
           r44_bg49 = map2_dbl(r44, lm_49, ~ purrr::possibly(predict, NA_real_)(.y, newdata = list(max_44 = .x)))) |>
    ## mutate_at(.vars = vars(one_of("s44", "r44")),
    ##           .funs = list(#bg44 = ~ {{fac}} * (. * slope_44 + intercept_44),
    ##             bg45 = ~ {{fac}} * predict(lm_45, newdata = list(max_44 = .)),
    ##             bg46 = ~ {{fac}} * predict(lm_46, newdata = list(max_44 = .)),
    ##             bg47 = ~ {{fac}} * predict(lm_47, newdata = list(max_44 = .)),
    ##             bg48 = ~ {{fac}} * predict(lm_48, newdata = list(max_44 = .)),
    ##             bg49 = ~ {{fac}} * predict(lm_49, newdata = list(max_44 = .)))) |>
    mutate(
      ## s44_bg = ifelse(is.na(s44_bg44), s44, s44 - s44_bg44),
      s45_bg = ifelse(is.na(s44_bg45), s45, s45 - {{fac}} * s44_bg45),
      s46_bg = ifelse(is.na(s44_bg46), s46, s46 - {{fac}} * s44_bg46),
      s47_bg = ifelse(is.na(s44_bg47), s47, s47 - {{fac}} * s44_bg47),
      s48_bg = ifelse(is.na(s44_bg48), s48, s48 - {{fac}} * s44_bg48),
      s49_bg = ifelse(is.na(s44_bg49), s49, s49 - {{fac}} * s44_bg49),
      ## r44_bg = ifelse(is.na(r44_bg44), r44, r44 - r44_bg44),
      r45_bg = ifelse(is.na(r44_bg45), r45, r45 - {{fac}} * r44_bg45),
      r46_bg = ifelse(is.na(r44_bg46), r46, r46 - {{fac}} * r44_bg46),
      r47_bg = ifelse(is.na(r44_bg47), r47, r47 - {{fac}} * r44_bg47),
      r48_bg = ifelse(is.na(r44_bg48), r48, r48 - {{fac}} * r44_bg48),
      r49_bg = ifelse(is.na(r44_bg49), r49, r49 - {{fac}} * r44_bg49))

  if (sum(is.na(out$s44_bg47)) > 0) {
    warning(glue::glue("{sum(is.na(out$s44_bg47))} out of {nrow(out)} intensities could not be assigned a background scan! Investigate!"))
  }

  out
}

parse_preparation_number <- function(data, col = sheet) {
  sheet <- NULL
  data |>
    tidylog::mutate(Preparation = str_extract({{col}}, "_\\d+_") |>
             str_replace_all("_", "") |>
             parse_double())
}

# This convers the list to a simple string vector for easier export.
string_scan_files <- function(data) {
  data |>
    tidylog::mutate(scan_files = paste0(scan_files) |>
             stringr::str_replace_all("c?\\(?\\\\?\",?\\)?", ""))
}

# a special version of clumpedr's add_info that does not rely on Analysis
add_scan_info <- function(data, .info, cols, quiet = clumpedr:::default(quiet)) {
  if (nrow(data) == 0) {
    return(tibble(file_id = character()))
  }

  if (!"file_id" %in% cols) {
    cols <- c("file_id", cols)
  }

  if (!quiet) {
    message("Info: appending measurement information.")
  }

  left_join(x = data, y = .info %>% select(tidyselect::all_of(cols)), 
            by = "file_id")
}

# This was the easiest way I could find to create consistent output with the
# desired order of columns.
export_scan_metadata <- function(data, meta, file) {
  
  if ("file_datetime" %in% colnames(data)) {
    data <- data |> 
      tidylog::filter(!file_datetime %in% meta$file_datetime)
  } else {
    warning("Column `file_datetime` not found in data.", call. = FALSE)
  }

  if ("file_id" %in% colnames(data)) {
    data <- data |>
      tidylog::filter(!file_id %in% meta$file_id)
  }  else {
    warning("Column `file_id` not found in data.", call. = FALSE)
  }
  
   data |>
     #tidylog::filter(scan_datetime > max(meta$scan_datetime, na.rm = TRUE)) |>
     tidylog::select(any_of(c("file_id",
                              "file_root",
                              "file_datetime",
                              "voltage",
                              "voltage_max",
                              "min_44",
                              "min_45",
                              "min_46",
                              "min_47",
                              "min_48",
                              "min_49",
                              "min_54",
                              "max_44",
                              "max_45",
                              "max_46",
                              "max_47",
                              "max_48",
                              "max_49",
                              "max_54",
                              "scan_group",
                              "scan_datetime",
                              "bg_group",
                              "scan_files",
                              "scan_n",
                              "scan_duration",
                              ## "intercept_45",
                              ## "intercept_46",
                              ## "intercept_47",
                              ## "intercept_48",
                              ## "intercept_49",
                              ## "slope_45",
                              ## "slope_46",
                              ## "slope_47",
                              ## "slope_48",
                              ## "slope_49",
                              "min",
                              "max",
                              "min_start_44",
                              "min_end_44",
                              "min_start_45_49",
                              "min_end_45_49",
                              "max_start",
                              "max_end",
                              "manual_outlier",
                              "manual_notes",
                              "fix_software",
                              "scan_group_overwrite",
                              "voltage_overwrite",
                              "checked_by",
                              "checked_date",
                              "checked_comment"))) |>
     writexl::write_xlsx(file)
   file
}


# raw deltas --------------------------------------------------------------

# Most functions to calculate raw deltas are already a part of ~clumpedr~
# https://github.com/isoverse/clumpedr/
filter_duplicated_raw_cycles <- function(data) {
  if (nrow(data) == 0L) {
    return(tibble(file_id = character()))
  }
  tidylog::distinct(data, Analysis, file_id, type, cycle, v44.mV, .keep_all = TRUE)
}

add_mineralogy <- function(data, info) {
  if (nrow(data) == 0L) {
    return(tibble(file_id = character()))
  }

  data |>
    tidylog::left_join(select(info, file_id, Mineralogy), by = "file_id")
}

add_R18 <- function(data, min = Mineralogy) {
  if (nrow(data) == 0L) {
    return(tibble(file_id = character()))
  }

  data |>
    tidylog::mutate(R18_PDB = case_when(is.na({{min}}) ~ #{
      ## warning("No mineralogy specified, defaulting to Calcite") ;
      clumpedr:::default(R18_PDB), #},
      {{min}} %in% "Calcite" ~ clumpedr:::default(R18_PDB),
      {{min}} %in% "Aragonite" ~ 1.00909,
      {{min}} %in% "Dolomite" ~ NA_real_, #{ warning("No R18 available for Dolomite"); NA_real_ },
      !is.na({{min}}) ~ NA_real_ #{ warning("Incorrect Mineralogy"); NA_real_ }
      ))
}

# summarize d45 d46 d47 d48 d49 d13C d18O D45 D46 D47 D48 D49 param_49
summarize_d13C_d18O_D47 <- function(data) {
  if (nrow(data) == 0L) {
    return(tibble(file_id = character()))
  }

  if (!"cycle_data" %in% names(data)) {
    stop("'cycle_data' not found in data.")
  }

  data |>
    ## group_by(file_id) |>
    mutate(summaries = map(.data$cycle_data,
                           .f = ~ .x |>
                             filter(!outlier, !outlier_cycle) |>
                             dplyr::select(d45, d46, d47, d48, d49,
                                           d13C_PDB, d18O_PDB,
                                           D45_raw, D46_raw, D47_raw, D48_raw, D49_raw,
                                           param_49) |>
                             dplyr::summarize_all(list(
                                      n = ~ length(.),  # get the number of cycles excluding the outliers
                                      mean = ~ mean(., na.rm = TRUE),
                                      sd = ~ sd(., na.rm = TRUE))) |>
                             # TODO: rewrite using dplyr 1.0.0's across()
                             mutate(n_ok = d45_n, d45_n = NULL, d46_n = NULL, # n is the same for all
                                    d47_n = NULL, d48_n = NULL,  d49_n = NULL,
                                    d13C_PDB_n = NULL, d18O_PDB_n = NULL,
                                    D45_raw_n = NULL, D46_raw_n = NULL,
                                    D47_raw_n = NULL, D48_raw_n = NULL,
                                    D49_raw_n = NULL, param_49_n = NULL,
                                    d13C_PDB_sem = d13C_PDB_sd / sqrt(n_ok - 1),
                                    d18O_PDB_sem = d18O_PDB_sd / sqrt(n_ok - 1),
                                    D47_raw_sem = D47_raw_sd / sqrt(n_ok - 1),
                                    d13C_PDB_cl = qt((1 - 0.05), n_ok - 1) * d13C_PDB_sem,
                                    d18O_PDB_cl = qt((1 - 0.05), n_ok - 1) * d18O_PDB_sem,
                                    D47_raw_cl = qt((1 - 0.05), n_ok - 1) * D47_raw_sem,
                                    d13C_PDB_lwr = d13C_PDB_mean - d13C_PDB_cl,
                                    d18O_PDB_lwr = d18O_PDB_mean - d18O_PDB_cl,
                                    D47_raw_lwr = D47_raw_mean - D47_raw_cl,
                                    d13C_PDB_upr = d13C_PDB_mean + d13C_PDB_cl,
                                    d18O_PDB_upr = d18O_PDB_mean + d18O_PDB_cl,
                                    D47_raw_upr = D47_raw_mean + D47_raw_cl))) |>
    unnest(cols = c(summaries))
}

##' Rolling offset correction
##'
##' Calculates the offset of standards with respect to their accepted values.
##' Then takes a rolling mean of this offset and applies it to the data. This
##' will get rid of inter-preparation drift. Note that error propagation is not
##' implemented at the moment!
##'
##' @param data
##' @param std The standard(s) to perform offset correction with.
##' @param grp A string with the column name to group by
##' @param exp The expected/accepted values to append to the data.
##' @param raw The raw data column to use for calculation.
##' @param off The name of the new offset column.
##' @param off_good The name of the new column of offset values that are not outliers and are in =std=.
##' @param off_avg The name of the new moving average of the off_good column.
##' @param cor The name of the new offset-corrected column.
##' @param width The width of the moving average of the offset.
##' @param out The name of the outlier_offset column.
##' @param min The minimum offset to determine whether it's an outlier_offset.
##' @param max The maximum offset to determine whether it's an outlier_offset.
offset_correction <- function(data, std = "ETH-3", grp = NULL,
                              exp, raw, off, off_good,
                              off_avg, cor,
                              ## off_bin = offset_bin_D47, dur = 1.5 * 3600,
                              width = 7, out, min = 0.5, max = 0.9, quiet = clumpedr:::default(quiet)) {
  if (nrow(data) == 0L) {
    return(tibble(file_id = character()))
  }

  if (!quiet) {
    message("Info: performing rolling offset correction.\n",
            glue::glue("\ni raw = {quo_name(enquo(raw))}\n"),
            #glue::glue("i std = {glue::glue_collapse(unique(std), sep = ' ', last = '\n')}")
            #glue::glue("i std = {std}\n"),
            #glue::glue("i width = {width}\n"),
            # because of the quosures I don't know how to do this.
            glue::glue("\ni grp = {grp}")
            )
  }

  D47_offset_std <- expected_D47 <- D47_raw_mean <- D47_offset_average <- D47_offset_corrected <- NULL

  prm <- purrr::possibly(zoo::rollmean, NA_real_)

  if (is.null(grp) || is.na(grp)) {
    data |>
      mutate({{off}} := {{exp}} - {{raw}},
             {{out}} := {{off}} < {{min}} | {{off}} > {{max}}) |>
      ## summarize_outlier() |>
      mutate({{off_good}} := ifelse(!.data$outlier & 
                                      (.data$broadid %in% std),
                                    {{off}}, 
                                    NA_real_),
             ## {{off_bin}} := seq_along(findInterval(file_datetime - dur, file_datetime)),
             {{off_avg}} := prm({{off_good}}, 
                                {{width}}, 
                                na.rm = TRUE, 
                                fill = "extend"),
             ## {{off_avg}} := zoo::rollapplyr({{off_good}}, {{off_bin}}, mean, na.rm = TRUE, fill = NA_real_),
             {{cor}} := {{raw}} + {{off_avg}})
  } else {
    data |>
      mutate({{off}} := {{exp}} - {{raw}},
             {{out}} := {{off}} < {{min}} | {{off}} > {{max}}) |>
      ## summarize_outlier() |>
      group_by_at(grp) |>
      mutate({{off_good}} := ifelse(!.data$outlier &
                                      (.data$broadid %in% {{std}}),
                                    {{off}},
                                    NA_real_),
             {{off_avg}} := prm({{off_good}}, 
                                unique({{width}}),
                                na.rm = TRUE, 
                                fill = "extend"),
             {{cor}} := {{raw}} + {{off_avg}}) |>
      ungroup()
  }
}

##' Apply offset correction
##'
##' This applies [offset_correction()] to \eqn{\delta^{13}C}{δ13C},
##' \eqn{\delta^{18}O}{δ18O}, and \eqn{\Delta_{47}}{Δ47}
##'
##' @param acc A tibble/dataframe with accepted values.
##' @param par A tibble/dataframe with paramters `grp`, `width`, and `std`.
offset_correction_wrapper <- function(data, acc) {
  if (nrow(data) == 0L) {
    return(tibble(file_id = character()))
  }

  prm <- purrr::possibly(zoo::rollmean, NA_real_)

  data |>
    clumpedr::append_expected_values(std_names = acc$id, 
                                     by = broadid,
                                     std_values = acc$D47, 
                                     exp = expected_D47) |>
    offset_correction(std = str_split(.data$off_D47_stds, 
                                      " ", 
                                      simplify = TRUE),
                      grp = "preparation",
                      exp = .data$expected_D47,
                      raw = .data$D47_raw_mean,
                      off = D47_offset,
                      off_good = D47_offset_good,
                      off_avg = D47_offset_average,
                      cor = D47_offset_corrected,
                      width = .data$off_D47_width,
                      out = outlier_offset_D47,
                      min = .data$off_D47_min,
                      max = .data$off_D47_max) |>
    group_by(.data$preparation, .data$Line) |>
    mutate(D47_offset_average_line = prm(.data$D47_offset_good, 
                                         .data$off_D47_width * 2, 
                                         na.rm = TRUE, 
                                         fill = "extend"),
           D47_offset_corrected_line = .data$D47_raw_mean + 
             .data$D47_offset_average_line) |>
    ungroup() |>
    clumpedr::append_expected_values(std_names = acc$id, 
                                     by = broadid,
                                     std_values = acc$d13C, 
                                     exp = expected_d13C) |>
    offset_correction(std = str_split(.data$off_d13C_stds, 
                                      " ", 
                                      simplify = TRUE),
                      grp = "preparation",
                      exp = .data$expected_d13C,
                      raw = .data$d13C_PDB_mean,
                      off = d13C_offset,
                      off_good = d13C_offset_good,
                      off_avg = d13C_offset_average,
                      cor = d13C_offset_corrected,
                      width = .data$off_d13C_width,
                      out = outlier_offset_d13C,
                      min = .data$off_d13C_min,
                      max = .data$off_d13C_max) |>
    group_by(.data$Line) |>
    mutate(d13C_offset_average_line = prm(.data$d13C_offset_good, 
                                          .data$off_d13C_width * 2, 
                                          na.rm = TRUE, 
                                          fill = "extend"),
           d13C_offset_corrected_line = .data$d13C_PDB_mean + 
             .data$d13C_offset_average_line) |>
    ungroup() |>
    # d18O
    clumpedr::append_expected_values(std_names = acc$id, 
                                     by = broadid,
                                     std_values = acc$d18O, 
                                     exp = expected_d18O) |>
    offset_correction(std = str_split(.data$off_d18O_stds, 
                                      " ",
                                      simplify = TRUE),
                      grp = "preparation",
                      exp = .data$expected_d18O,
                      raw = .data$d18O_PDB_mean,
                      off = d18O_offset,
                      off_good = d18O_offset_good,
                      off_avg = d18O_offset_average,
                      cor = d18O_offset_corrected,
                      width = .data$off_d18O_width,
                      out = outlier_offset_d18O,
                      min = .data$off_d18O_min,
                      max = .data$off_d18O_max) |>
    group_by(.data$Line) |>
    mutate(d18O_offset_average_line = prm(.data$d18O_offset_good, 
                                          .data$off_d18O_width * 2, 
                                          na.rm = TRUE, 
                                          fill = "extend"),
           d18O_offset_corrected_line = .data$d18O_PDB_mean +
             .data$d18O_offset_average_line) |>
    ungroup()
}

# The empirical transfer function relates the raw D47 values of the standards to
# their expected values. Here we apply a rolling version, that is affected by
# the ~width~ measurements that bracket the current one.
rolling_etf <- function(data,
                        x = expected_D47,
                        y = D47_offset_corrected,
                        slope = etf_slope,
                        intercept = etf_intercept,
                        std = paste0("ETH-", 1:3), width = 201,
                        grp = etf_grp,
                        quiet = clumpedr:::default(quiet)) {
  ## if (nrow(data) == 0L) {
  ##   return(tibble(file_id = character()))
  ## }

  if (!quiet) message(glue::glue("Info: calculating rolling empirical transfer function based on non-outlier standards {glue::glue_collapse(distinct(data, {{std}}), sep = ' ')} {quo_name(enquo(y))} values with width = {glue::glue_collapse(distinct(data, {{width}}), sep = ' ')}, grouped by {quo_name(enquo(grp))}"))

  ## lengths <- pull(data, {{width}})
  ## if (unique(lengths) == 1L) {
  ##   message("only one window size, simplifying parameter")
  ##   lengths <- unique(lengths)
  ## }

  data |>
    group_by({{grp}}) |>
    mutate(
      x_good = ifelse(!outlier & broadid %in% str_split({{std}}, " ", simplify = TRUE),
                      {{x}}, NA_real_),
      y_good = ifelse(!outlier, {{y}}, NA_real_),
      starts = row_number() - floor({{width}} / 2),
      stops = row_number() + floor({{width}} / 2),
      fit = slider::hop(cur_data(), # cur_data ensures I'm within a group
                purrr::possibly(~ lm(y_good ~ x_good, data = .),
                                list(coefficients = c("(Intercept)" = NA, "y_good" = NA))),
                .starts = starts,
                .stops = stops),
      # perhaps these two are the culprits that crash my laptop?
      {{intercept}} := map_dbl(fit, ~ coef(.x)[[1]]),
      {{slope}} := map_dbl(fit, ~ coef(.x)[[2]])) |>
    ungroup({{grp}}) |>
    tidylog::select(-one_of("x_good", "y_good", "fit"))
}

summarise_cycle_outliers <- function(data) {
  message("Info: summarizing cycle outliers")
  data |>
    mutate(
      # the number of cycles, including the outlier cycles (compare to n_ok)
      n_cyc = map_dbl(cycle_data,
                      purrr::possibly(~ .x |>
                                        select(cycle) |>
                                        max(na.rm = TRUE),
                                      NA_real_)),
      prop_bad_cycles = map_dbl(cycle_data,
                                purrr::possibly(~ sum(.$outlier_cycle, na.rm = TRUE), NA_real_)) / n_cyc,
      outlier_noscan = is.na(scan_group),
      outlier_nodelta = is.na(d47_mean),
      outlier_cycles = prop_bad_cycles > .data$prop_bad_cyc,
      ## prop_bad_param49s = map_dbl(cycle_data,
      ##                             purrr::possibly(~ sum(.$outlier_param49, na.rm = TRUE), NA_real_)) / n_cyc,
      ## outlier_param49 = param_49_mean > p49_crit | param_49_mean < -p49_crit,
      outlier_internal_sd_D47_raw = D47_raw_sd > .data$sd_D47,
      outlier_internal_sd_d13C_PDB = d13C_PDB_sd > .data$sd_d13C,
      outlier_internal_sd_d18O_PDB = d18O_PDB_sd > .data$sd_d18O) #|>
    ## mutate(manual_outlier = ifelse(is.na(manual_outlier), FALSE, manual_outlier)) |>
    ## rename(outlier_manual = manual_outlier) |>
    ## clumpedr::summarise_outlier(quiet = TRUE)
    ## mutate(outlier = outlier_noscan | outlier_nodelta | (!is.na(outlier_cycles) & outlier_cycles))
}

# This is to simply represent in one column why a particular measurement could
# be an outlier.
create_reason_for_outlier <- function(data) {
  data |>
    tidylog::mutate(reason_for_outlier =
                      paste0(ifelse(outlier_manual, paste("manual", ifelse(!is.na(checked_comment), checked_comment, " no_comment "), "\n"), ""),
                             ifelse(outlier_nodelta, "  noδ\n", ""),
                             ifelse(outlier_noscan, "  noscn\n", ""),
                             ifelse(is.na(outlier_init), "  init_NA\n", ""),
                             ifelse(!is.na(outlier_init) & outlier_init, "  init\n", ""),
                             ifelse(!is.na(outlier_s44_init_low) & outlier_s44_init_low, "    s44_low\n", ""),
                             ifelse(!is.na(outlier_r44_init_low) & outlier_r44_init_low, "    r44_low\n", ""),
                             ifelse(!is.na(outlier_s44_init_high) & outlier_s44_init_high, "    s44_high\n", ""),
                             ifelse(!is.na(outlier_r44_init_high) & outlier_r44_init_high, "    r44_high\n", ""),
                             ifelse(!is.na(outlier_i44_init_diff) & outlier_i44_init_diff, "    i44_diff\n", ""),
                             ## ifelse(is.na(outlier_cycles), "  cyc_NA\n", ""),
                             ifelse(!is.na(outlier_cycles) & outlier_cycles, "  cyc\n", ""),
                             ## ifelse(is.na(outlier_param49), "  p49_NA\n", ""),
                             ifelse(!is.na(outlier_param49) & outlier_param49, "  p49\n", ""),
                             ifelse(!is.na(outlier_internal_sd_D47_raw) & outlier_internal_sd_D47_raw, "  D47_sd\n", ""),
                             ifelse(!is.na(outlier_internal_sd_d13C_PDB) & outlier_internal_sd_d13C_PDB, "  d13C_sd\n", ""),
                             ifelse(!is.na(outlier_internal_sd_d18O_PDB) & outlier_internal_sd_d18O_PDB, "  d18O_sd\n", ""),
                             ifelse(!is.na(outlier_offset_D47) & outlier_offset_D47, "  D47_off\n", ""),
                             ifelse(!is.na(outlier_offset_d13C) & outlier_offset_d13C, "  d13C_off\n", ""),
                             ifelse(!is.na(outlier_offset_d18O) & outlier_offset_d18O, "  d18O_off\n", "")))
}


# export ------------------------------------------------------------------
order_columns <- function(data, extra = NULL) {
  data |>
    tidylog::select(tidyselect::one_of(c(
      # we want these all the way in the beginning for easy access and column blocking
      "Analysis",
      "file_id",
      "broadid",
      "masspec",

      # metadata from file_info
      "file_datetime",
      "time_diff",
      "file_root",
      "file_path",
      "file_subpath",
      "file_size",
      "timeofday",
      "Row",
      "Peak Center",
      "Background",
      "Pressadjust",
      "Reference Refill",
      "Line",
      "Sample",
      "Weight [mg]",
      "weight",
      "Identifier 1",
      "identifier_1",
      "Identifier 2",
      "identifier_2",
      "Comment",
      "comment",

      "Preparation",
      "preparation",
      "time_prep",
      "dir_prep",
      "Method",

      # meas_info and it's parsed components
      "measurement_info",
      "acid_temperature",
      "ref_mbar",
      "ref_pos",
      "bellow_pos_smp",
      "init_int",
      "background",
      "PC",
      "VM1_aftr_trfr",
      "CO2_after_exp",
      "no_exp",
      "total_CO2",
      "p_gases",
      "p_no_acid",
      "extra_drops",
      "leak_rate",
      "MS_integration_time.s",

      # background scan components
      "bg_group",
      "scan_group",
      "scan_datetime",
      "bg_fac",

      ## "intercept_45",
      ## "intercept_46",
      ## "intercept_47",
      ## "intercept_48",
      ## "intercept_49",
      ## "slope_45",
      ## "slope_46",
      ## "slope_47",
      ## "slope_48",
      ## "slope_49",
      "outlier_noscan",

      "cycle_data",

      # anything related to cycle disabling
      "dis_min",
      "dis_max",
      "dis_fac",
      "dis_rel",
      "cycle_has_drop",
      "n_ok",
      "n_cyc",
      "prop_bad_cycles", # proportion of outlier_cycle
      "prop_bad_cyc",
      "outlier_cycles",

      # raw values
      "d45_mean",
      "d46_mean",
      "d47_mean",
      "d48_mean",
      "d49_mean",
      # little delta
      "d45_sd",
      "d46_sd",
      "d47_sd",
      "d48_sd",
      "d49_sd",

      "outlier_nodelta",

      "R18_PDB", # the value used in calculations, based on mineralogy

      "d13C_PDB_mean",
      "d18O_PDB_mean",

      "d13C_PDB_sd",
      "d18O_PDB_sd",
      "d13C_PDB_sem",
      "d18O_PDB_sem",
      "d13C_PDB_cl",
      "d18O_PDB_cl",
      "d13C_PDB_lwr",
      "d18O_PDB_lwr",
      "d13C_PDB_upr",
      "d18O_PDB_upr",

      # ref gas values
      "d13C_PDB_wg",
      "d18O_PDBCO2_wg",

      # internal sd
      "sd_d13C",
      "outlier_internal_sd_d13C_PDB",
      "sd_d18O",
      "outlier_internal_sd_d18O_PDB",

      # offset correction
      "accepted_d13C",
      "d13C_offset",
      "off_d13C_min",
      "off_d13C_max",
      "outlier_offset_d13C",
      "d13C_offset_good",
      "off_d13C_grp",
      "off_d13C_width",
      "off_d13C_stds",
      "d13C_offset_average",
      "d13C_offset_corrected",
      "d13C_offset_average_line",
      "d13C_offset_corrected_line",

      "accepted_d18O",
      "d18O_offset",
      "off_d18O_min",
      "off_d18O_max",
      "outlier_offset_d18O",
      "d18O_offset_good",
      "off_d18O_grp",
      "off_d18O_width",
      "off_d18O_stds",
      "d18O_offset_average",
      "d18O_offset_corrected",
      "d18O_offset_average_line",
      "d18O_offset_corrected_line",

      "D45_raw_mean",
      "D46_raw_mean",
      "D47_raw_mean",
      "D48_raw_mean",
      "D49_raw_mean",

      "D45_raw_sd",
      "D46_raw_sd",
      "D47_raw_sd",
      "D48_raw_sd",
      "D49_raw_sd",
      "D47_raw_sem",
      "D47_raw_cl",
      "D47_raw_lwr",
      "D47_raw_upr",

      # internal sd outliers
      "sd_D47",
      "outlier_internal_sd_D47_raw",

      "expected_D47",
      "D47_offset",
      "off_D47_min",
      "off_D47_max",
      "outlier_offset_D47",
      "off_D47_grp",
      "off_D47_stds",
      "D47_offset_good",
      "off_D47_width",
      "D47_offset_average",
      "D47_offset_corrected",
      "D47_offset_average_line",
      "D47_offset_corrected_line",

      "param_49_mean",
      "param_49_sd",
      # param 49 related stuff
      "p49_crit",
      "prop_bad_param49s",
      "prop_bad_param49",
      "outlier_param49",

      # anything related to initial intensity
      # values
      "s44_init",
      "r44_init",
      # criteria
      "init_low",
      "init_high",
      "init_diff",
      # outlier
      "outlier_s44_init_low",
      "outlier_r44_init_low",
      "outlier_s44_init_high",
      "outlier_r44_init_high",
      "outlier_i44_init_diff",
      "outlier_init",

      # empirical transfer function
      "etf_grp",
      "etf_stds",
      "etf_width",
      "etf_slope_raw", # rolling no offset
      "etf_intercept_raw",
      "etf_slope", # rolling + offset correction
      "etf_intercept",
      "etf_slope_grp", # sessions
      "etf_intercept_grp",
      "etf_slope_grp_off", # sessions + offset correction
      "etf_intercept_grp_off",

      ## "D47_70_deg",
      ## "D47_70_deg_raw",

      # acid fractionation
      "acid_fractionation_factor",
      "D47_final", # session + offset correction
      "D47_final_roll", # rolling + offset correction
      "D47_final_no_offset", # session
      "D47_final_roll_no_offset", # rolling

      "temperature_slope",
      "temperature_intercept",
      "temperature",
      "temperature_no_offset",

      ## extra
      "outlier",
      "reason_for_outlier",

      # metadata fixes that we need to be at the end for easy inspection
      "outlier_manual",
      "Preparation_overwrite",
      "Identifier 1_overwrite",
      "Identifier 2_overwrite",
      "Weight [mg]_overwrite",
      "Comment_overwrite",
      "scan_group_overwrite",
      "Mineralogy",
      "checked_by",
      "checked_date",
      "checked_comment")))
}

add_remaining_meta <- function(data, meta) {
  if (nrow(data) == 0L) {
    return(tibble(file_id = character()))
  }

  prefer_data <- c("bg_group", "bg_fac", "scan_group", "scan_datetime", "scan_group_overwrite",
                   "scan_files", "scan_n", "scan_duration",
                   "d13C_PDB_wg", "d18O_PDBCO2_wg",
                   "Mineralogy")

  data |>
    ## select(-one_of("Analysis")) |> # some are giving us issues!
    # we use a full join so that files that don't have any raw data are still included in the list!
    tidylog::full_join(meta |>
                       # remove the columns that are already there in the file info itself
                       select(-any_of(prefer_data))
                       # ,
                       # we intentionally do NOT specify by what the matching occurs, since the columns
                       # differ slightly between machines and file types
                       ## by = c("file_id",
                       ## ##        ## "Analysis",
                       ##        "bg_group",
                       ## ##        ## "bg_fac",
                       ##        "scan_group_overwrite",
                       ##        "scan_group",
                       ##        "scan_datetime",
                       ##        "scan_files",
                       ##        "scan_n",
                       ##        "scan_duration",
                       ##        ## "lm_45",
                       ##        ## "lm_46",
                       ##        ## "lm_47",
                       ##        ## "lm_48",
                       ##        ## "lm_49",
                       ##        ## "intercept_45",
                       ##        ## "intercept_46",
                       ##        ## "intercept_47",
                       ##        ## "intercept_48",
                       ##        ## "intercept_49",
                       ##        ## "slope_45",
                       ##        ## "slope_46",
                       ##        ## "slope_47",
                       ##        ## "slope_48",
                       ##        ## "slope_49",
                       ## "d13C_PDB_wg",
                       ## "d18O_PDBCO2_wg",
                       ## "Mineralogy")
                       )
}

tar_excel <- function(dat, file) {
  dat |>
    tidylog::filter(!is.na(Analysis)) |>
    rename(manual_outlier = outlier_manual) |>
    writexl::write_xlsx(path = file)
  file
}

tar_csv <- function(dat, file) {
  dat |>
    tidylog::filter(!is.na(Analysis)) |>
    ## rename(manual_outlier = outlier_manual) |> # do not rename for widget
    readr::write_csv(file = file)
  file
}

tar_write  <- function(dat, file) {
  readr::write_rds(dat, file)
  file
}

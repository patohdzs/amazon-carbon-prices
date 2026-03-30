# > PROJECT INFO
# NAME: CARBON PRICES AND FOREST PRESERVATION OVER SPACE AND TIME IN THE BRAZILIAN AMAZON
# LEAD: JULIANO ASSUNÇÃO, LARS PETER HANSEN, TODD MUNSON, JOSÉ A. SCHEINKMAN
#
# > THIS SCRIPT
# AIM: CLEAN RAW DATA SHAPEFILE OF MUNI DIVISION (2015)
# AUTHOR: JOÃO PEDRO VIEIRA
#
# > NOTES
# 1: -

library(tidyverse)
library(tictoc)
library(sf)
library(stringi)

# START TIMER
tic(msg = "muniDivision2015_raw2clean.R script", log = TRUE)

# Read shapefile
raw_muni <- st_read(
  dsn = "data/raw/ibge/muni_division_2015",
  layer = "BRMUE250GC_SIR"
)

# Translate column names
raw_muni <-
  raw_muni %>%
  rename(
    muni_code = CD_GEOCMU,
    muni_name = NM_MUNICIP
  )

# Add state_uf column
raw_muni <-
  raw_muni %>%
  mutate(state_uf = case_when(
    str_sub(muni_code, 1, 2) == 11 ~ "RO",
    str_sub(muni_code, 1, 2) == 12 ~ "AC",
    str_sub(muni_code, 1, 2) == 13 ~ "AM",
    str_sub(muni_code, 1, 2) == 14 ~ "RR",
    str_sub(muni_code, 1, 2) == 15 ~ "PA",
    str_sub(muni_code, 1, 2) == 16 ~ "AP",
    str_sub(muni_code, 1, 2) == 17 ~ "TO",
    str_sub(muni_code, 1, 2) == 21 ~ "MA",
    str_sub(muni_code, 1, 2) == 22 ~ "PI",
    str_sub(muni_code, 1, 2) == 23 ~ "CE",
    str_sub(muni_code, 1, 2) == 24 ~ "RN",
    str_sub(muni_code, 1, 2) == 25 ~ "PB",
    str_sub(muni_code, 1, 2) == 26 ~ "PE",
    str_sub(muni_code, 1, 2) == 27 ~ "AL",
    str_sub(muni_code, 1, 2) == 28 ~ "SE",
    str_sub(muni_code, 1, 2) == 29 ~ "BA",
    str_sub(muni_code, 1, 2) == 31 ~ "MG",
    str_sub(muni_code, 1, 2) == 32 ~ "ES",
    str_sub(muni_code, 1, 2) == 33 ~ "RJ",
    str_sub(muni_code, 1, 2) == 35 ~ "SP",
    str_sub(muni_code, 1, 2) == 41 ~ "PR",
    str_sub(muni_code, 1, 2) == 42 ~ "SC",
    str_sub(muni_code, 1, 2) == 43 ~ "RS",
    str_sub(muni_code, 1, 2) == 50 ~ "MS",
    str_sub(muni_code, 1, 2) == 51 ~ "MT",
    str_sub(muni_code, 1, 2) == 52 ~ "GO",
    str_sub(muni_code, 1, 2) == 53 ~ "DF"
  ))

# Class - muni_code should be numeric
lapply(raw_muni, class)

raw_muni <- raw_muni %>% mutate(muni_code = as.numeric(muni_code))

# LATIN CHARACTER TREATMENT
raw_muni <-
  raw_muni %>%
  mutate(
    across(
      where(is.character),
      \(x) stri_trans_general(
        stri_enc_toutf8(x, is_unknown_8bit = TRUE),
        "Latin-ASCII"
      )
    )
  )

# LETTERS CAPITALIZATION
raw_muni <-
  raw_muni %>%
  mutate(muni_name = toupper(muni_name))

# PROJECTION
# SIRGAS 2000 / Brazil Polyconic (https://epsg.io/5880)
raw_muni <- st_transform(x = raw_muni, crs = 5880)

# Save as Rdata
save(raw_muni, file = "data/clean/raw_muni.Rdata")

# REMOVE POLYGONS IDENTIFIED AS BODY OF WATERS AND NOT MUNICIPALITIES
# see muniDivision2007_raw2clean
raw_muni <- raw_muni %>% filter(!muni_code %in% c(4300001, 4300002))

# GEOMETRY CLEANUP
raw_muni <- st_make_valid(raw_muni)


# Set labels
set_label(raw_muni$muni_code) <- "municipality code (7-digit, IBGE)"
set_label(raw_muni$muni_name) <- "municipality name"
set_label(raw_muni$state_uf) <- "state name (abbreviation)"

# Change object name before saving
muni_division_2015 <- raw_muni

# Save data set
out_path <- "data/clean/muni_division_2015.Rdata"
save(muni_division_2015, file = out_path)

# END TIMER
toc(log = TRUE)

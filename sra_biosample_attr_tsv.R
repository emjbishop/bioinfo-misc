# This script takes a partially completed biosample attributes TSV file and
# edits it to be compatible with SRA requirements

library(tidyverse)

raw <- read_tsv("~/Downloads/Model.organism.animal.1.0.tsv", skip = 11)

out <- raw %>%
  separate(
    `*sample_name`,
    sep = "_",
    remove = FALSE,
    # The sample_name column is ignored when checking that rows are unique,
    # so split that into new columns so they'll be unique
    into = c("ptid", "replicate_id", NA, "condition", "data_type")
  ) %>%
  # Convert to the correct date format (YYYY-MM-DD)
  mutate(`*collection_date` = lubridate::mdy(`*collection_date`))

write_tsv(out, file = "~/Downloads/model_organism.tsv")


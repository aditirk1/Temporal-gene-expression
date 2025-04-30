# Load GEO IDs from metadata file
library(tidyverse)
dataset_meta <- read.csv("dataset_meta.csv")
library(dplyr)
library(tidyr)

# Load the data (assuming it's already read into a variable called 'data')
dataset_meta <- read.csv("dataset_meta.csv")

data$Timepoint..after.stressor. <- as.character(data$Timepoint..after.stressor.)
dataset_meta <- dataset_meta %>%
  mutate(
    Sample_IDs = strsplit(`Relevant.Sample.IDs`, ","),
    Control_Treatment = strsplit(`Control.Treatment`, ","),
  ) %>%
  unnest(cols = c(Sample_IDs, Control_Treatment)) %>%
  select(-`Relevant.Sample.IDs`, -`Control.Treatment`)

# Save the new expanded data to a CSV file
write.csv(expanded_data, "path/to/your/expanded_dataset_meta.csv", row.names = FALSE)


  geo_ids <- read.csv("dataset_meta.csv") %>%
  pull(GEO.Dataset.ID) %>%
  unique()

# Main function
driver_script <- function(geo_ids) {
  source("worker_script.R")
  source("filtering.R")
  
  # Process each GEO dataset
  for (geo_id in geo_ids) {
    tryCatch({
      download_and_extract_counts(geo_id)  # Worker script
      process_gse(geo_id)                 # Filtering script
    }, error = function(e) {
      cat("Error processing GEO ID:", geo_id, "-", e$message, "\n")
    })
  }
}

driver_script(geo_ids)
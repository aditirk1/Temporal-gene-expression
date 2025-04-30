# Load required libraries
library(httr)
library(R.utils)

# Function to download both raw counts and series matrix files for given GSE IDs
download_files <- function(gse_ids, output_dir = "raw_counts") {
  # Read GSE IDs from the input file
  geo_ids <- readLines(gse_ids)
  
  # Remove any empty lines or whitespace
  geo_ids <- trimws(geo_ids)
  geo_ids <- geo_ids[geo_ids != ""]
  
  # Verify we have GSE IDs to process
  if (length(geo_ids) == 0) {
    stop("No GSE IDs found in the input file")
  }
  
  cat("Found", length(geo_ids), "GSE IDs to process\n")
  
  # Create main output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Loop through each GEO ID to download and extract the files
  for (geo_id in geo_ids) {
    # Create a directory for this specific GSE ID
    gse_dir <- file.path(output_dir, geo_id)
    if (!dir.exists(gse_dir)) {
      dir.create(gse_dir, recursive = TRUE)
    }
    
    # 1. Download raw counts file
    raw_counts_url <- paste0(
      "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=", 
      geo_id, 
      "&format=file&file=", 
      geo_id, 
      "_raw_counts_GRCh38.p13_NCBI.tsv.gz"
    )
    
    raw_counts_gz <- file.path(gse_dir, paste0(geo_id, "_raw_counts.tsv.gz"))
    
    # Attempt download with error handling for raw counts
    tryCatch({
      # Download the raw counts file
      GET(raw_counts_url, write_disk(raw_counts_gz, overwrite = TRUE))
      cat("Successfully downloaded raw counts for", geo_id, "\n")
      
      # Extract the .gz file to .tsv format
      raw_counts_tsv <- file.path(gse_dir, paste0(geo_id, "_raw_counts.tsv"))
      gunzip(raw_counts_gz, destname = raw_counts_tsv, remove = TRUE)
      cat("Successfully extracted raw counts file:", raw_counts_tsv, "\n")
    }, error = function(e) {
      warning("Error downloading or extracting raw counts file for ", geo_id, ": ", e$message)
    })
    
    # 2. Download series matrix file
    # Extract GSE number from the ID
    gse_number <- as.numeric(gsub("GSE", "", geo_id))
    
    # Determine the appropriate subdirectory based on GSE number
    gse_prefix <- substr(gse_number, 1, nchar(gse_number) - 3)  # Extract first few digits
    subdir <- paste0("GSE", gse_prefix, "nnn/", geo_id)
    
    
    series_matrix_url <- paste0(
      "https://ftp.ncbi.nlm.nih.gov/geo/series/", 
      subdir, 
      "/matrix/", 
      geo_id, 
      "_series_matrix.txt.gz"
    )
    
    series_matrix_gz <- file.path(gse_dir, paste0(geo_id, "_series_matrix.txt.gz"))
    
    # Attempt download with error handling for series matrix
    tryCatch({
      # Download the series matrix file
      GET(series_matrix_url, write_disk(series_matrix_gz, overwrite = TRUE))
      cat("Successfully downloaded series matrix for", geo_id, "\n")
      
      # Extract the .gz file to .txt format
      series_matrix_txt <- file.path(gse_dir, paste0(geo_id, "_series_matrix.txt"))
      gunzip(series_matrix_gz, destname = series_matrix_txt, remove = TRUE)
      cat("Successfully extracted series matrix file:", series_matrix_txt, "\n")
    }, error = function(e) {
      warning("Error downloading or extracting series matrix file for ", geo_id, ": ", e$message)
    })
  }
}

# Main execution
gse_ids <- file.path(getwd(), "gse_ids.txt")
download_files(gse_ids)

filtering_data <- function(raw_counts_dir = "raw_counts", output_dir = "filtered_counts", min_samples = 0.5, cpm_threshold = NULL) {
  # Load required libraries
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Please install the 'edgeR' package first: install.packages('BiocManager'); BiocManager::install('edgeR')")
  }
  library(edgeR)
  
  # Create output directory for filtered data if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  # Create a log directory within output dir
  log_dir <- file.path(output_dir, "logs")
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE)
    cat("Created log directory:", log_dir, "\n")
  }
  
  # Create a directory for plots
  plot_dir <- file.path(output_dir, "plots")
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
    cat("Created plot directory:", plot_dir, "\n")
  }
  
  # List GSE folders in the raw counts directory
  gse_folders <- list.dirs(raw_counts_dir, full.names = TRUE, recursive = FALSE)
  
  # If no GSE folders found, look for TSV files directly in the main directory
  if (length(gse_folders) == 0) {
    cat("No GSE folders found, looking for TSV files in the main directory.\n")
    gse_folders <- raw_counts_dir
  }
  
  # Process each GSE folder
  for (gse_folder in gse_folders) {
    # Get GSE ID from folder name
    gse_id <- basename(gse_folder)
    cat("\n===== Processing", gse_id, "=====\n")
    
    # Create output subfolder for this GSE
    gse_output_dir <- file.path(output_dir, gse_id)
    if (!dir.exists(gse_output_dir)) {
      dir.create(gse_output_dir, recursive = TRUE)
    }
    
    # List all .tsv files in the GSE folder
    tsv_files <- list.files(gse_folder, pattern = "\\.tsv$", full.names = TRUE)
    
    # If no files found in the GSE folder
    if (length(tsv_files) == 0) {
      cat("No TSV files found in", gse_folder, "\n")
      next
    }
    
    # Log file for this GSE
    log_file <- file.path(log_dir, paste0(gse_id, "_filtering_log.txt"))
    sink(log_file, split = TRUE)
    
    # Process each TSV file
    for (file in tsv_files) {
      file_basename <- basename(file)
      cat("\nProcessing file:", file_basename, "\n")
      cat("Full path:", file, "\n")
      
      # Read the raw counts file
      tryCatch({
        count_data <- read.table(file, header = TRUE, sep = "\t", check.names = FALSE)
        
        # Check if the first column is GeneID and handle accordingly
        if (colnames(count_data)[1] == "GeneID") {
          cat("GeneID column detected as first column.\n")
          # Save gene IDs and use them as row names
          gene_ids <- count_data$GeneID
          # Remove GeneID column and set row names
          count_data <- count_data[, -1, drop = FALSE]
          rownames(count_data) <- gene_ids
        }
        
        # Ensure all data is numeric
        count_data_numeric <- apply(count_data, 2, function(x) as.numeric(as.character(x)))
        rownames(count_data_numeric) <- rownames(count_data)
        
        # Check for NAs and report
        na_count <- sum(is.na(count_data_numeric))
        if (na_count > 0) {
          cat("Warning: Found", na_count, "NA values in the data. These will be treated as zeros.\n")
          count_data_numeric[is.na(count_data_numeric)] <- 0
        }
        
        # Basic data stats
        cat("Dimensions of raw count matrix:", nrow(count_data_numeric), "genes x", ncol(count_data_numeric), "samples\n")
        
        # Calculate summary statistics for raw counts
        gene_totals <- rowSums(count_data_numeric)
        sample_totals <- colSums(count_data_numeric)
        
        cat("Sample total counts summary:\n")
        print(summary(sample_totals))
        
        cat("Gene total counts summary:\n")
        print(summary(gene_totals))
        
        # Calculate genes with zero counts
        zero_count_genes <- sum(gene_totals == 0)
        cat("Number of genes with zero counts:", zero_count_genes, "\n")
        
        # Calculate CPM using edgeR
        y <- DGEList(counts = count_data_numeric)
        cpm_data <- cpm(y)
        
        # If no threshold provided, calculate one
        if (is.null(cpm_threshold)) {
          # Calculate CPM percentiles
          cpm_vector <- as.vector(cpm_data)
          cpm_vector <- cpm_vector[cpm_vector > 0]  # exclude zeros
          
          # Calculate percentile-based thresholds
          percentiles <- c(0.1, 1, 5, 10, 25, 50)
          threshold_candidates <- quantile(cpm_vector, probs = percentiles/100, na.rm = TRUE)
          names(threshold_candidates) <- paste0(percentiles, "%")
          
          cat("\nCPM Percentile Thresholds:\n")
          print(threshold_candidates)
          
          # Choose 1st percentile as default threshold
          chosen_threshold <- threshold_candidates["1%"]
          cat("\nAutomatically selected CPM threshold:", chosen_threshold, "(1st percentile)\n")
        } else {
          chosen_threshold <- cpm_threshold
          cat("\nUsing user-provided CPM threshold:", chosen_threshold, "\n")
        }
        
        # Generate CPM distribution plot
        png(file.path(plot_dir, paste0(gse_id, "_", gsub("\\.tsv$", "", file_basename), "_cpm_dist.png")), 
            width = 1200, height = 800, res = 100)
        par(mfrow = c(2, 1))
        
        # Histogram of log2 CPM
        hist(log2(cpm_data + 1), 
             breaks = 100,
             main = paste0("Distribution of log2(CPM + 1) for ", gse_id, " - ", file_basename),
             xlab = "log2(CPM + 1)",
             col = "steelblue")
        abline(v = log2(chosen_threshold + 1), col = "red", lwd = 2)
        text(log2(chosen_threshold + 1), 0, paste("Threshold:", round(chosen_threshold, 3)), 
             pos = 4, col = "red", cex = 1.2)
        
        # Density plot of log2 CPM with threshold line
        plot(density(log2(cpm_data[cpm_data > 0] + 1)), 
             main = paste0("Density of log2(CPM + 1) for ", gse_id, " - ", file_basename),
             xlab = "log2(CPM + 1)",
             col = "blue", lwd = 2)
        abline(v = log2(chosen_threshold + 1), col = "red", lwd = 2)
        dev.off()
        cat("Created CPM distribution plot in", plot_dir, "\n")
        
        # Calculate minimum number of samples required
        min_samples_required <- ceiling(min_samples * ncol(count_data_numeric))
        cat("Minimum samples required to pass filter:", min_samples_required, 
            "(", min_samples * 100, "% of", ncol(count_data_numeric), "samples)\n")
        
        # Filter genes based on CPM threshold
        pass_filter <- rowSums(cpm_data >= chosen_threshold) >= min_samples_required
        filtered_data <- count_data_numeric[pass_filter, , drop = FALSE]
        
        # Reporting
        cat("Number of genes before filtering:", nrow(count_data_numeric), "\n")
        cat("Number of genes after filtering:", nrow(filtered_data), "\n")
        cat("Percentage of genes retained:", round(nrow(filtered_data) / nrow(count_data_numeric) * 100, 2), "%\n")
        
        # Create a table showing filtering stats for each gene
        filtering_stats <- data.frame(
          Gene = rownames(count_data_numeric),
          TotalCounts = rowSums(count_data_numeric),
          MaxCPM = apply(cpm_data, 1, max),
          SamplesAboveThreshold = rowSums(cpm_data >= chosen_threshold),
          PassFilter = pass_filter
        )
        
        # Save filtering stats
        stats_file <- file.path(gse_output_dir, paste0("stats_", file_basename))
        write.table(filtering_stats, stats_file, sep = "\t", quote = FALSE, row.names = FALSE)
        cat("Saved gene filtering statistics to", stats_file, "\n")
        
        # Save filtered data 
        output_file <- file.path(gse_output_dir, paste0("filtered_", file_basename))
        write.table(filtered_data, output_file, sep = "\t", quote = FALSE)
        cat("Saved filtered data to", output_file, "\n")
        
        # Create a summary plot of filtering
        png(file.path(plot_dir, paste0(gse_id, "_", gsub("\\.tsv$", "", file_basename), "_filtering_summary.png")), 
            width = 1000, height = 800, res = 100)
        par(mfrow = c(2, 2))
        
        # Plot 1: Filter pass/fail
        barplot(c(Retained = sum(pass_filter), Filtered = sum(!pass_filter)), 
                main = "Genes Retained vs Filtered",
                col = c("darkgreen", "darkred"),
                ylab = "Number of genes")
        
        # Plot 2: Samples above threshold per gene (histogram)
        hist(filtering_stats$SamplesAboveThreshold, 
             breaks = seq(0, ncol(count_data_numeric), by = 1),
             main = "Number of Samples Above CPM Threshold Per Gene",
             xlab = "Number of samples",
             col = "steelblue")
        abline(v = min_samples_required, col = "red", lwd = 2)
        
        # Plot 3: Total counts before/after
        boxplot(list(Before = log2(gene_totals + 1), 
                     After = log2(rowSums(filtered_data) + 1)),
                main = "Gene Total Counts (log2)",
                col = c("lightblue", "darkgreen"),
                ylab = "log2(Total Counts + 1)")
        
        # Plot 4: CPM values before/after
        boxplot(list(Before = log2(rowMeans(cpm_data) + 1), 
                     After = log2(rowMeans(cpm(filtered_data)) + 1)),
                main = "Mean CPM Values (log2)",
                col = c("lightblue", "darkgreen"),
                ylab = "log2(Mean CPM + 1)")
        
        dev.off()
        cat("Created filtering summary plots in", plot_dir, "\n")
        
      }, error = function(e) {
        cat("ERROR processing file", file_basename, ":", conditionMessage(e), "\n")
      })
    }
    
    # Close the log file
    sink()
  }
  
  cat("\nFiltering complete. Results saved to", output_dir, "\n")
  cat("Log files saved to", log_dir, "\n")
  cat("Plot files saved to", plot_dir, "\n")
}
filtering_data()
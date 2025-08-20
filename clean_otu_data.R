# OTU Data Cleaning Script
# This script cleans sparse microbiome OTU data for faster analysis

# Set working directory to the correct location
setwd("C:/MicrobiomeBot/Microbot-")

# Load required libraries (only dplyr, using base R for file I/O)
library(dplyr)

# Function to clean OTU data
clean_otu_data <- function(input_file = "otu_table_with_species.csv", 
                          output_file = "otu_table_cleaned.csv",
                          min_prevalence = 0.05,  # Feature must be present in at least 5% of samples
                          min_abundance = 10,     # Minimum total abundance per feature
                          max_features = 200,     # Maximum number of features to keep
                          remove_low_samples = TRUE) {
  
  cat("Loading OTU data from:", input_file, "\n")
  
  # Read the OTU table using base R
  otu_data <- read.csv(input_file, row.names = 1, check.names = FALSE)
  
  cat("Original data dimensions:", nrow(otu_data), "features x", ncol(otu_data), "samples\n")
  
  # Calculate sparsity
  total_cells <- nrow(otu_data) * ncol(otu_data)
  zero_cells <- sum(otu_data == 0)
  sparsity <- (zero_cells / total_cells) * 100
  cat("Data sparsity:", round(sparsity, 1), "% zeros\n")
  
  # Step 1: Remove features with very low prevalence
  min_samples <- max(2, round(ncol(otu_data) * min_prevalence))
  feature_prevalence <- rowSums(otu_data > 0)
  prevalent_features <- feature_prevalence >= min_samples
  
  cat("Step 1: Removing features present in <", min_samples, "samples\n")
  otu_cleaned <- otu_data[prevalent_features, ]
  cat("Kept", nrow(otu_cleaned), "features after prevalence filtering\n")
  
  # Step 2: Remove features with very low total abundance
  feature_abundance <- rowSums(otu_cleaned)
  abundant_features <- feature_abundance >= min_abundance
  
  cat("Step 2: Removing features with total abundance <", min_abundance, "\n")
  otu_cleaned <- otu_cleaned[abundant_features, ]
  cat("Kept", nrow(otu_cleaned), "features after abundance filtering\n")
  
  # Step 3: Keep only top abundant features if still too many
  if (nrow(otu_cleaned) > max_features) {
    cat("Step 3: Keeping top", max_features, "most abundant features\n")
    feature_totals <- rowSums(otu_cleaned)
    top_features <- names(sort(feature_totals, decreasing = TRUE)[1:max_features])
    otu_cleaned <- otu_cleaned[top_features, ]
    cat("Kept", nrow(otu_cleaned), "top abundant features\n")
  }
  
  # Step 4: Remove low-abundance samples (optional)
  if (remove_low_samples) {
    sample_totals <- colSums(otu_cleaned)
    abundance_threshold <- quantile(sample_totals, 0.1)  # Bottom 10%
    good_samples <- sample_totals > abundance_threshold
    
    cat("Step 4: Removing samples with abundance <=", round(abundance_threshold), "\n")
    otu_cleaned <- otu_cleaned[, good_samples]
    cat("Kept", ncol(otu_cleaned), "samples after sample filtering\n")
  }
  
  # Step 5: Final cleanup - remove any remaining all-zero rows/columns
  otu_cleaned <- otu_cleaned[rowSums(otu_cleaned) > 0, ]
  otu_cleaned <- otu_cleaned[, colSums(otu_cleaned) > 0]
  
  # Calculate final sparsity
  final_total_cells <- nrow(otu_cleaned) * ncol(otu_cleaned)
  final_zero_cells <- sum(otu_cleaned == 0)
  final_sparsity <- (final_zero_cells / final_total_cells) * 100
  
  cat("\nFinal cleaned data dimensions:", nrow(otu_cleaned), "features x", ncol(otu_cleaned), "samples\n")
  cat("Final sparsity:", round(final_sparsity, 1), "% zeros\n")
  cat("Data reduction:", round((1 - (final_total_cells / total_cells)) * 100, 1), "% size reduction\n")
  
  # Save cleaned data using base R
  write.csv(otu_cleaned, output_file, row.names = TRUE)
  cat("Cleaned data saved to:", output_file, "\n")
  
  # Return summary statistics
  return(list(
    original_dims = c(nrow(otu_data), ncol(otu_data)),
    cleaned_dims = c(nrow(otu_cleaned), ncol(otu_cleaned)),
    original_sparsity = sparsity,
    cleaned_sparsity = final_sparsity,
    size_reduction = (1 - (final_total_cells / total_cells)) * 100
  ))
}

# Run the cleaning function
cat("=== OTU Data Cleaning Script ===\n")
results <- clean_otu_data()

cat("\n=== Cleaning Summary ===\n")
cat("Original:", results$original_dims[1], "x", results$original_dims[2], "\n")
cat("Cleaned:", results$cleaned_dims[1], "x", results$cleaned_dims[2], "\n")
cat("Sparsity reduced from", round(results$original_sparsity, 1), "% to", round(results$cleaned_sparsity, 1), "%\n")
cat("Dataset size reduced by", round(results$size_reduction, 1), "%\n")
cat("Cleaning completed successfully!\n")

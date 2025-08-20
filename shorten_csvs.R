# R Script to Shorten CSV Files for Microbiome Analysis
# Creates meaningful subsets of large datasets while preserving biological diversity

# Using base R functions - no additional packages required

# Set working directory
setwd("C:/MicrobiomeBot/Microbot-")

# Function to calculate coefficient of variation for OTU selection
calculate_cv <- function(x) {
  if(mean(x, na.rm = TRUE) == 0) return(0)
  sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE)
}

# 1. SHORTEN INTERACTIONS.CSV (2671 -> 100 rows)
# Strategy: Select top interactions by compression value + diverse interaction types
cat("Processing interactions.csv...\n")
interactions <- read.csv("interactions.csv", stringsAsFactors = FALSE)

# Get top interactions by compression value (most significant)
interactions_ordered <- interactions[order(interactions$comp, decreasing = TRUE), ]
top_by_compression <- interactions_ordered[1:60, ]

# Get diverse interaction types
interaction_types <- unique(interactions$lnk)
diverse_interactions <- data.frame()

for(type in interaction_types) {
  type_subset <- interactions[interactions$lnk == type, ]
  n_sample <- min(10, nrow(type_subset))
  if(n_sample > 0) {
    sampled <- type_subset[sample(nrow(type_subset), n_sample), ]
    diverse_interactions <- rbind(diverse_interactions, sampled)
  }
}

# Take first 40 from diverse interactions
diverse_interactions <- diverse_interactions[1:min(40, nrow(diverse_interactions)), ]

# Combine and ensure exactly 100 rows
interactions_combined <- rbind(top_by_compression, diverse_interactions)
interactions_short <- interactions_combined[!duplicated(interactions_combined), ]
interactions_short <- interactions_short[1:min(100, nrow(interactions_short)), ]

write.csv(interactions_short, "interactions_short.csv", row.names = FALSE)
cat("Created interactions_short.csv with", nrow(interactions_short), "rows\n")

# 2. SHORTEN TAXONOMY_TABLE.CSV (1505 -> 100 rows)
# Strategy: Maximize taxonomic diversity across all levels
cat("Processing taxonomy_table.csv...\n")
taxonomy <- read.csv("taxonomy_table.csv", stringsAsFactors = FALSE)

# Calculate taxonomic diversity score for each OTU
taxonomy$diversity_score <- 
  (!is.na(taxonomy$Kingdom)) * 1 +
  (!is.na(taxonomy$Phylum)) * 2 +
  (!is.na(taxonomy$Class)) * 3 +
  (!is.na(taxonomy$Order)) * 4 +
  (!is.na(taxonomy$Family)) * 5 +
  (!is.na(taxonomy$Genus)) * 6

taxonomy_scored <- taxonomy[order(taxonomy$diversity_score, decreasing = TRUE), ]

# Select representatives from each phylum
phylum_list <- unique(taxonomy_scored$Phylum[!is.na(taxonomy_scored$Phylum)])
phylum_representatives <- data.frame()

for(phylum in phylum_list) {
  phylum_subset <- taxonomy_scored[!is.na(taxonomy_scored$Phylum) & taxonomy_scored$Phylum == phylum, ]
  n_select <- min(5, nrow(phylum_subset))
  if(n_select > 0) {
    selected <- phylum_subset[1:n_select, ]
    phylum_representatives <- rbind(phylum_representatives, selected)
  }
}

# Fill remaining slots with high diversity score OTUs
remaining_slots <- 100 - nrow(phylum_representatives)
if(remaining_slots > 0) {
  # Get OTUs not already selected
  selected_ids <- phylum_representatives[, 1]
  remaining_otus <- taxonomy_scored[!(taxonomy_scored[, 1] %in% selected_ids), ]
  additional_otus <- remaining_otus[1:min(remaining_slots, nrow(remaining_otus)), ]
  taxonomy_short <- rbind(phylum_representatives, additional_otus)
} else {
  taxonomy_short <- phylum_representatives[1:100, ]
}

# Remove diversity score column
taxonomy_short$diversity_score <- NULL

write.csv(taxonomy_short, "taxonomy_table_short.csv", row.names = FALSE)
cat("Created taxonomy_table_short.csv with", nrow(taxonomy_short), "rows\n")

# 3. SHORTEN SPECIES_MAPPING.CSV (1505 -> 100 rows)
# Strategy: Match with selected taxonomy OTUs
cat("Processing species_mapping.csv...\n")
species_mapping <- read.csv("species_mapping.csv", stringsAsFactors = FALSE)

# Get OTU IDs from shortened taxonomy
selected_otu_ids <- taxonomy_short[, 1]  # First column contains OTU IDs

# Filter species mapping to match selected OTUs
species_mapping_short <- species_mapping[species_mapping$otu_id %in% selected_otu_ids, ]
species_mapping_short <- species_mapping_short[1:min(100, nrow(species_mapping_short)), ]

write.csv(species_mapping_short, "species_mapping_short.csv", row.names = FALSE)
cat("Created species_mapping_short.csv with", nrow(species_mapping_short), "rows\n")

# 4. SHORTEN OTU_TABLE_NUMERIC.CSV (1505 rows x ~400 cols -> 100x100)
# Strategy: High variance OTUs + representative samples
cat("Processing otu_table_numeric.csv...\n")
otu_numeric <- read.csv("otu_table_numeric.csv", stringsAsFactors = FALSE, check.names = FALSE)

# Calculate variance for each OTU (row)
otu_data <- otu_numeric[, -1]  # Remove first column (OTU names)
otu_variances <- apply(otu_data, 1, function(x) var(as.numeric(x), na.rm = TRUE))

# Select top 100 OTUs by variance (most variable = most informative)
high_var_indices <- order(otu_variances, decreasing = TRUE)[1:min(100, length(otu_variances))]
selected_otus <- otu_numeric[high_var_indices, ]

# Select 100 most diverse samples (columns)
selected_otu_data <- selected_otus[, -1]
sample_sums <- colSums(selected_otu_data, na.rm = TRUE)

# Select samples with good coverage and diversity
sample_indices <- order(sample_sums, decreasing = TRUE)[1:min(100, length(sample_sums))]
sample_names <- names(sample_sums)[sample_indices]

# Create final shortened table
otu_numeric_short <- cbind(selected_otus[, 1, drop = FALSE], selected_otus[, sample_names])

write.csv(otu_numeric_short, "otu_table_numeric_short.csv", row.names = FALSE)
cat("Created otu_table_numeric_short.csv with", nrow(otu_numeric_short), "rows and", ncol(otu_numeric_short)-1, "sample columns\n")

# 5. SHORTEN OTU_TABLE_WITH_SPECIES.CSV (1505 rows x ~400 cols -> 100x100)
# Strategy: Match with selected high-variance OTUs
cat("Processing otu_table_with_species.csv...\n")
otu_species <- read.csv("otu_table_with_species.csv", stringsAsFactors = FALSE, check.names = FALSE)

# Get the row indices that correspond to high variance OTUs
# Match by position since both files should have same OTU order
otu_species_short <- cbind(otu_species[high_var_indices, 1, drop = FALSE], 
                          otu_species[high_var_indices, sample_names])

write.csv(otu_species_short, "otu_table_with_species_short.csv", row.names = FALSE)
cat("Created otu_table_with_species_short.csv with", nrow(otu_species_short), "rows and", ncol(otu_species_short)-1, "sample columns\n")

# 6. CREATE ENHANCED INTERACTIONS (if interactions_enhanced.csv exists)
if(file.exists("interactions_enhanced.csv")) {
  cat("Processing interactions_enhanced.csv...\n")
  interactions_enhanced <- read.csv("interactions_enhanced.csv", stringsAsFactors = FALSE)
  
  # Match with shortened interactions
  if("sp1" %in% names(interactions_enhanced) && "sp2" %in% names(interactions_enhanced)) {
    # Create matching logic
    match_indices <- c()
    for(i in 1:nrow(interactions_enhanced)) {
      for(j in 1:nrow(interactions_short)) {
        if(interactions_enhanced$sp1[i] == interactions_short$sp1[j] && 
           interactions_enhanced$sp2[i] == interactions_short$sp2[j]) {
          match_indices <- c(match_indices, i)
          break
        }
      }
    }
    
    if(length(match_indices) > 0) {
      interactions_enhanced_short <- interactions_enhanced[match_indices, ]
      interactions_enhanced_short <- interactions_enhanced_short[1:min(100, nrow(interactions_enhanced_short)), ]
      
      write.csv(interactions_enhanced_short, "interactions_enhanced_short.csv", row.names = FALSE)
      cat("Created interactions_enhanced_short.csv with", nrow(interactions_enhanced_short), "rows\n")
    }
  }
}

# Summary Report
cat("\n=== SHORTENING SUMMARY ===\n")
cat("All files reduced to ~100 rows/columns while preserving:\n")
cat("✓ Interactions: Top compression values + interaction type diversity\n")
cat("✓ Taxonomy: Maximum taxonomic diversity across all levels\n") 
cat("✓ Species Mapping: Matched with selected taxonomic OTUs\n")
cat("✓ OTU Tables: High variance OTUs + high coverage samples\n")
cat("✓ Biological meaningfulness maintained through strategic selection\n")
cat("\nFiles created:\n")
cat("- interactions_short.csv\n")
cat("- taxonomy_table_short.csv\n") 
cat("- species_mapping_short.csv\n")
cat("- otu_table_numeric_short.csv\n")
cat("- otu_table_with_species_short.csv\n")
if(file.exists("interactions_enhanced.csv")) {
  cat("- interactions_enhanced_short.csv\n")
}
cat("\nReady for fast microbiome analysis!\n")

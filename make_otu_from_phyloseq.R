suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(phyloseq)
})

# -------------------- CLI --------------------
option_list <- list(
  make_option(c("--rda"),     type="character",            help="Path to phyloseq .rda file"),
  make_option(c("--outdir"),  type="character", default=".", help="Output directory")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$rda)) stop("Please provide --rda")

rda_path <- normalizePath(opt$rda, mustWork = TRUE)
outdir   <- normalizePath(opt$outdir, mustWork = FALSE)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# -------------------- Load phyloseq --------------------
env <- new.env()
load(rda_path, envir = env)

ps <- NULL
for (obj in ls(env)) {
  if (methods::is(env[[obj]], "phyloseq")) { ps <- env[[obj]]; break }
}
if (is.null(ps)) stop("No phyloseq object found in RDA")

# Get OTU/ASV matrix with taxa as rows
otu <- as(otu_table(ps), "matrix")
if (taxa_are_rows(ps) == FALSE) {
  otu <- t(otu)
}
otu <- as.data.frame(otu, check.names = FALSE)

# Ensure row/col names are clean & unique
if (is.null(rownames(otu))) rownames(otu) <- paste0("taxon_", seq_len(nrow(otu)))
rownames(otu) <- make.names(rownames(otu), unique = TRUE)
colnames(otu) <- make.names(colnames(otu), unique = TRUE)

# Coerce to numeric safely
otu[] <- lapply(otu, function(x) as.numeric(as.character(x)))
otu[is.na(otu)] <- 0

# -------------------- Write FULL (unchanged) --------------------
full_path <- file.path(outdir, "otu_table_full.csv")
fwrite(cbind(Taxon = rownames(otu), otu), full_path)
cat(sprintf("✅ Wrote full OTU table: %s  (%d taxa × %d samples)\n",
            full_path, nrow(otu), ncol(otu)))

# -------------------- Build 50×10 IMPORTANT subset --------------------
message("📌 Creating 50×10 important subset for otu_table.csv ...")

mat <- as.matrix(otu)
# Rank taxa by abundance + variability
row_sums <- rowSums(mat)
row_vars <- apply(mat, 1, var)
# Avoid division by zero
rs_max <- ifelse(max(row_sums) > 0, max(row_sums), 1)
rv_max <- ifelse(max(row_vars) > 0, max(row_vars), 1)
taxa_score <- (row_sums / rs_max) + (row_vars / rv_max)

n_taxa_target   <- min(50, nrow(mat))
top_taxa_names  <- names(sort(taxa_score, decreasing = TRUE))[seq_len(n_taxa_target)]

# Rank samples by abundance + diversity (non-zero count)
col_sums <- colSums(mat)
col_div  <- apply(mat > 0, 2, sum)
cs_max <- ifelse(max(col_sums) > 0, max(col_sums), 1)
cd_max <- ifelse(max(col_div)  > 0, max(col_div),  1)
sample_score <- (col_sums / cs_max) + (col_div / cd_max)

n_samples_target <- min(10, ncol(mat))
top_sample_names <- names(sort(sample_score, decreasing = TRUE))[seq_len(n_samples_target)]

sub <- mat[top_taxa_names, top_sample_names, drop = FALSE]

# Drop any-empty taxa/samples (should be rare with the scoring, but safe)
if (nrow(sub) > 0 && ncol(sub) > 0) {
  sub <- sub[rowSums(sub) > 0, , drop = FALSE]
  sub <- sub[, colSums(sub) > 0, drop = FALSE]
}

# If after pruning we lost too many, fall back to a looser selection:
if (nrow(sub) == 0 || ncol(sub) == 0) {
  message("⚠️ Strict filter removed everything. Falling back to abundance-only ranking.")
  # Top taxa by abundance
  top_taxa_names  <- names(sort(row_sums, decreasing = TRUE))[seq_len(n_taxa_target)]
  # Top samples by abundance
  top_sample_names <- names(sort(col_sums, decreasing = TRUE))[seq_len(n_samples_target)]
  sub <- mat[top_taxa_names, top_sample_names, drop = FALSE]
  # Keep non-zero rows/cols if possible
  if (nrow(sub) > 0 && ncol(sub) > 0) {
    sub <- sub[rowSums(sub) > 0, , drop = FALSE]
    sub <- sub[, colSums(sub) > 0, drop = FALSE]
  }
}

# -------------------- Write SHORT (important 50×10) --------------------
short_path <- file.path(outdir, "otu_table.csv")
short_df <- data.frame(Taxon = rownames(sub), check.names = FALSE)
short_df <- cbind(short_df, as.data.frame(sub, check.names = FALSE))
fwrite(short_df, short_path)

cat(sprintf("✅ Wrote SHORT important OTU table: %s  (%d taxa × %d samples, non-zero prioritized)\n",
            short_path, nrow(sub), ncol(sub)))

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(InfIntE)
})

# ---- Options ----
option_list <- list(
  make_option(c("--otu"), type="character", help="Path to OTU table CSV"),
  make_option(c("--outdir"), type="character", default=".", help="Output directory")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$otu)) stop("Please provide --otu")

otu_path <- normalizePath(opt$otu, mustWork = TRUE)
outdir   <- normalizePath(opt$outdir, mustWork = FALSE)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ---- Load OTU table ----
otu <- fread(otu_path, check.names = FALSE, data.table = FALSE)

# If first column is taxa IDs
if (!is.numeric(otu[[1]])) {
  rownames(otu) <- make.names(otu[[1]], unique = TRUE)
  otu[[1]] <- NULL
} else {
  rownames(otu) <- paste0("taxon_", seq_len(nrow(otu)))
}
otu[] <- lapply(otu, function(x) as.numeric(as.character(x)))
otu <- as.matrix(otu)
otu[is.na(otu)] <- 0

# ---- Run InfIntE ----
ncores <- if (.Platform$OS.type == "windows") 1 else min(6, parallel::detectCores())
nperms <- 10   # fixed at 10

res <- tryCatch({
  infinte(otu_tb = otu, exclusion = TRUE, ncores = ncores, nperms = nperms)
}, error = function(e) {
  warning("InfIntE error, retrying with safe settings...")
  infinte(otu_tb = otu, exclusion = TRUE, ncores = 1, nperms = 5)
})

# ---- Save selected interactions ----
sel <- res$selected_interactions
sel <- sel[, c("sp1", "sp2", "lnk", "comp")]
sel$comp <- as.numeric(sel$comp)

out_path <- file.path(outdir, "interact_short.csv")
fwrite(sel, out_path)

cat(sprintf("✅ Wrote %d interactions to %s (from %d taxa × %d samples, nperms=%d)\n",
            nrow(sel), out_path, nrow(otu), ncol(otu), nperms))

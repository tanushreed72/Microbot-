suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(phyloseq)
  library(InfIntE)
})

option_list <- list(
  make_option(c("--rda"),    type="character", help="Path to .rda containing a phyloseq object"),
  make_option(c("--outdir"), type="character", default=".", help="Output directory"),
  make_option(c("--nperms"), type="integer",   default=10, help="InfIntE permutations [default: 10]")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$rda)) stop("Please provide --rda")

rda_path <- normalizePath(opt$rda, mustWork = TRUE)
outdir   <- normalizePath(opt$outdir, mustWork = FALSE)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ---- Load phyloseq object ----
env <- new.env()
loaded <- load(rda_path, envir = env)
ps_names <- loaded[ vapply(loaded, function(nm) inherits(env[[nm]], "phyloseq"), logical(1)) ]
if (length(ps_names) == 0) stop("No phyloseq object found in the .rda")
ps <- env[[ ps_names[1] ]]

# ---- Extract OTU table (taxa as rows) ----
otu <- as(otu_table(ps), "matrix")
if (!taxa_are_rows(ps)) otu <- t(otu)
rownames(otu) <- make.names(rownames(otu), unique = TRUE)
colnames(otu) <- make.names(colnames(otu), unique = TRUE)
otu <- apply(otu, 2, as.numeric)
otu_df <- as.data.frame(otu)

# ---- Run InfIntE with Windows-safe cores ----
ncores <- if (.Platform$OS.type == "windows") 1L else min(6L, parallel::detectCores())
nperms <- as.integer(opt$nperms)

res <- tryCatch({
  infinte(otu_tb = otu_df, exclusion = TRUE, ncores = ncores, nperms = nperms)
}, error = function(e) {
  warning(sprintf("InfIntE error: %s; retrying with ncores=1, nperms=5", conditionMessage(e)))
  infinte(otu_tb = otu_df, exclusion = TRUE, ncores = 1, nperms = 5)
})

sel <- res$selected_interactions
if (is.null(sel) || nrow(sel) == 0) {
  sel <- data.frame(sp1=character(), sp2=character(), lnk=character(), comp=numeric())
} else {
  sel <- sel[, c("sp1", "sp2", "lnk", "comp")]
  sel$comp <- as.numeric(sel$comp)
}

out_path <- file.path(outdir, "interact.csv")
fwrite(sel, out_path)
cat(sprintf("✅ Wrote %d interactions to %s\n", nrow(sel), out_path))

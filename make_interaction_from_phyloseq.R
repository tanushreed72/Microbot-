suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(phyloseq)
  library(InfIntE)
})

option_list <- list(
  make_option(c("--rda"), type="character", help="Path to phyloseq .rda file"),
  make_option(c("--outdir"), type="character", default=".", help="Output directory"),
  make_option(c("--cells"), type="integer", default=500, help="Target number of nonzero cells (subset OTU before running InfIntE)")
)
opt <- parse_args(OptionParser(option_list=option_list))
if (is.null(opt$rda)) stop("Please provide --rda")

rda_path <- normalizePath(opt$rda, mustWork = TRUE)
outdir   <- normalizePath(opt$outdir, mustWork = FALSE)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ---- Load phyloseq ----
env <- new.env()
load(rda_path, envir = env)
ps <- NULL
for (obj in ls(env)) {
  if (methods::is(env[[obj]], "phyloseq")) {
    ps <- env[[obj]]
    break
  }
}
if (is.null(ps)) stop("No phyloseq object found in RDA")

otu <- as(otu_table(ps), "matrix")
otu <- as.data.frame(otu)

# ---- Subset to ~cells nonzero values ----
nz <- which(otu > 0, arr.ind = TRUE)
if (nrow(nz) > opt$cells) {
  set.seed(123)
  keep_idx <- sample(seq_len(nrow(nz)), opt$cells)
  nz_keep <- nz[keep_idx, , drop = FALSE]
  taxa_keep <- unique(rownames(otu)[nz_keep[,1]])
  samp_keep <- unique(colnames(otu)[nz_keep[,2]])
  otu <- otu[taxa_keep, samp_keep, drop = FALSE]
}
cat(sprintf("⚡ Using reduced OTU: %d taxa x %d samples (nonzero cells = %d)\n",
            nrow(otu), ncol(otu), sum(otu > 0)))

# ---- Run InfIntE (safe + fast settings) ----
ncores <- if (.Platform$OS.type == "windows") 1 else min(4, parallel::detectCores())
nperms <- 5  # keep small for speed

res <- tryCatch({
  infinte(otu_tb = otu, exclusion = TRUE, ncores = ncores, nperms = nperms)
}, error = function(e) {
  warning("InfIntE error, retrying with minimal settings...")
  infinte(otu_tb = otu, exclusion = TRUE, ncores = 1, nperms = 3)
})

sel <- res$selected_interactions
sel <- sel[, c("sp1", "sp2", "lnk", "comp")]
sel$comp <- as.numeric(sel$comp)

out_path <- file.path(outdir, "interact.csv")
fwrite(sel, out_path)
cat(sprintf("✅ Wrote %d interactions to %s\n", nrow(sel), out_path))

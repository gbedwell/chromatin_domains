library(Rsamtools)
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFiles)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hs1)

source("funs.R")

n.core <- 8

# Change path to wherever the data resides
bam.files <- list.files(
  path = "~/tsa_seq_hmm/son_data/",
  full.names = TRUE,
  pattern = ".bam"
)

pulldown.files <- bam.files[grepl(pattern = "pulldown", x = bam.files)]
input.files <- bam.files[grepl(pattern = "input", x = bam.files)]

BPPARAM <- MulticoreParam(workers = n.core, exportglobals = FALSE)

pulldown <- bin_counts(
  bam.files = pulldown.files,
  paired = FALSE,
  mid = FALSE,
  yield.size = 1E6,
  mapq.cutoff = 0,
  omit = "chrM",
  width = 20000,
  step = NULL,
  genome = Hsapiens,
  BPPARAM = BPPARAM
)

input <- bin_counts(
  bam.files = input.files,
  paired = FALSE,
  mid = FALSE,
  yield.size = 1E6,
  mapq.cutoff = 0,
  omit = "chrM",
  width = 20000,
  step = NULL,
  genome = Hsapiens,
  BPPARAM = BPPARAM
)

samp.names <- gsub(
  pattern = "_rmdup",
  replacement = "",
  x = gsub(
    pattern = c("_input_replicate|_pulldown_replicate"),
    replacement = "",
    x = colnames(pulldown)
  )
)

cell.types <- unique(gsub("_[0-9]+$", "", samp.names))

# 'rel.pulldown' refers to the normalization method described
# in DOI: 10.1101/gr.266239.120
# 'rel.input' performs the more standard log2(cpm_pulldown / cpm_input)
son.tsa <- calc_scores(
  pulldown = pulldown,
  input = input,
  samp.names = samp.names,
  norm.type = "rel.pulldown"
)

# By default, smoothing is accomplished with a Hamming window of length = 21.
# funs.R contains a custom smoothing function that has other options, as well.
son.tsa <- smooth_scores(norm.se = son.tsa, BPPARAM = BPPARAM)

# Resolution defines the granularity of the quantiles.
son.tsa <- calc_quantiles(se = son.tsa, resolution = 0.01, data.type = "smooth")

# Combine replicates by averaging.
combo <- combine_replicates(se = son.tsa, cell.types = cell.types)

combo <- calc_quantiles(se = combo, resolution = 0.01, data.type = "smooth")

# SPADs were defined by Belmont's lab as bins with the top 5% of SON TSA-seq scores
spads <- filter_quantiles(combo, threshold = 0.95, data.type = "smooth", collapse = FALSE)

# Does the same thing as above, but merges adjacent SPADs into a single continuous annotation.
spads.merged <- filter_quantiles(combo, threshold = 0.95, 
                                 data.type = "smooth", collapse = TRUE)

save(list = ls(all.names = TRUE), file = "son_tsa_analysis.RData.gz", compression_level = 6)


library(Rsamtools)
library(Biostrings)
library(GenomicRanges)
library(GenomicAlignments)
library(GenomicFiles)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hs1)
library(depmixS4)

source("funs.R")

# DamID-seq fragments should overlap GATC positions in the genome.
# Build GATC coordinate annotations for filtering mapped
chr.seqs <- get_chromosome_seqs(Hsapiens)

gatc <- lapply(
  X = seq_along(chr.seqs),
  FUN = function(x){
    target <- DNAStringSet("GATC", use.names = TRUE)
    for.pos <- matchPDict(target, chr.seqs[[x]])
    for.pos <- GRanges(
      seqnames = names(chr.seqs[x]),
      ranges = unlist(for.pos),
      strand = "*" 
    )
  }
)

# suppressWarnings() is used here to suppress the
# "No seqlevels in common" warning.
suppressWarnings(
  gatc <- do.call(c, gatc)
)

n.core <- 8

# Change path to wherever the data resides
bam.files <- list.files(
  path = "~/tsa_seq_hmm/laminB1_DamID_data/",
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
  filter.gr = gatc,
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
  filter.gr = gatc,
  BPPARAM = BPPARAM
)

samp.names <- gsub(
  pattern = "_trimmed_rmdup",
  replacement = "",
  x = gsub(
    pattern = c("_input_replicate|_pulldown_replicate"),
    replacement = "",
    x = colnames(pulldown)
  )
)

cell.types <- unique(gsub("_[0-9]+$", "", samp.names))

# 'rel.input' refers to the common normalization method
# log2(cpm_pulldown / cpm_input)
lamin.dam <- calc_scores(
  pulldown = pulldown,
  input = input,
  samp.names = samp.names,
  norm.type = "rel.input"
)

num.col <- ncol(assays(lamin.dam)$raw.scores)

# Fits a 2-state Gaussian HMM to data from each cell-type
# Uses depmixS4, which can utilize both replicates simultaneously
# global_hmm() uses a slightly modified random start procedure than the
# depmixS4 default. 
# See funs.R for more information.
lamin.dam.fits <- bplapply(
  X = seq(1, num.col, by = 2),
  FUN = function(x){
    mat <- assays(lamin.dam)$raw.scores[,x:(x+1)]
    mod <- global_hmm(
      df = data.frame(mat),
      gr <- rowRanges(lamin.dam),
      n.cond = 1,
      n.samp = 2,
      n.states = 2,
      multi.start = TRUE,
      n.start = 25,
      init.iters = 25
    )
    p.mod <- parse_model(mod.obj = mod)
    return(p.mod)
  },
  BPPARAM = MulticoreParam(workers = num.col, exportglobals = FALSE)
)

save(list = ls(all.names = TRUE), file = "lamin_dam_analysis.RData.gz", compression_level = 6)

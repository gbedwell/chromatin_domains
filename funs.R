# Get positions of aligned fragments
aligned_positions <- function(bam.file, paired = FALSE, mid = FALSE,
                              yield.size = 1E6, mapq.cutoff = 0, omit = NULL,
                              filter.gr = NULL, ...){
  if (isTRUE(paired)){
    as.mates = TRUE
  } else{
    as.mates = FALSE
    if (isTRUE(mid)){
      warning("mid cannot be TRUE for single-end reads. Setting to FALSE.")
      mid <- FALSE
    }
  }

  bam <- BamFile(bam.file, asMates = as.mates, yieldSize = yield.size)

  if (isTRUE(paired)){
    params <- ScanBamParam(
      flag = scanBamFlag(
        isProperPair = TRUE,
        isDuplicate = FALSE,
        isSecondaryAlignment = FALSE,
        isUnmappedQuery = FALSE
      ),
      mapqFilter = mapq.cutoff,
      ...
    )

    yield <- function(x){
      readGAlignmentPairs(x, param = params)
    }
  } else{
    params <- ScanBamParam(
      flag = scanBamFlag(
        isDuplicate = FALSE,
        isSecondaryAlignment = FALSE,
        isUnmappedQuery = FALSE
      ),
      mapqFilter = mapq.cutoff,
      ...
    )

    yield <- function(x){
      readGAlignments(x, param = params)
    }
  }

  map <- identity
  reduce <- c

  out <- reduceByYield(X = bam, YIELD = yield, MAP = map, REDUCE = reduce)

  df <- data.frame(seqnames(out), ranges(out), strand(out))
  colnames(df) <- c("seqnames", "start", "end", "width", "strand")

  gr <- GRanges(df)

  if (!isTRUE(mid)){
    plus <- gr[strand(gr) == "+"]
    minus <- gr[strand(gr) == "-"]
    end(plus) <- start(plus)
    start(minus) <- end(minus)
    gr <- c(plus, minus)
  } else{
    start(gr) <- floor((start(gr) + end(gr)) / 2)
    end(gr) <- start(gr)
  }

  if (!is.null(omit)){
    gr <- gr[!seqnames(gr) %in% omit]
  }

  gr <- sort(gr, ignore.strand = TRUE)
  
  if (!is.null(filter.gr)){
    gr <- gr[gr %over% filter.gr]
  }

  return(gr)
}

# Bin the genome
# The index metadata column serves as a bin counter for easier comparisons
make_bins <- function(genome, width, step = NULL, by.chrom = FALSE, omit = NULL){
  sl <- seqlengths(genome)
  sl <- sl[!names(sl) %in% omit]
  
  genome.gr <- GRanges(
    seqnames = names(sl),
    ranges = IRanges(start = 1, end = sl)
  )

  if (is.null(step) || step == 0){
    step = width
  }

  unindexed.bins <- slidingWindows(
    x = genome.gr,
    width = width,
    step = step
  )
  
  bins = GRangesList()
  start.pos <- 1
  for (i in seq_along(unindexed.bins)){
    tmp.gr <- unindexed.bins[[i]]
    tmp.gr$index <- seq(start.pos, (start.pos + (length(tmp.gr) - 1)), by = 1)
    bins[[i]] <- tmp.gr
    start.pos <- start.pos + length(tmp.gr)
  }

  if (!isTRUE(by.chrom)){
    bins <- unlist(bins)
  } else{
    names(bins) <- names(sl)
  }

  return(bins)
}

# Count the number of fragments in each bin
bin_overlaps <- function(aln.pos, bins){
  if (is.list(bins) || is(bins, "GRangesList")){
    counts <- lapply(
      X = bins,
      FUN = function(x){
        countOverlaps(x, aln.pos)
      }
    )
  } else{
    counts <- countOverlaps(bins, aln.pos)
  }

  return(counts)
}

# Wrapper function to parallelize bin counting
# Hard-code by.chrom = FALSE in binning
counts_wrapper <- function(bam.file, paired, mid, yield.size, mapq.cutoff, filter.gr,
                           omit, width, genome, step, ...){

  aln.gr <- aligned_positions(
    bam.file = bam.file,
    paired = paired,
    mid = mid,
    yield.size = yield.size,
    mapq.cutoff = mapq.cutoff,
    omit = omit,
    filter.gr = filter.gr,
    ...
  )

  genome.bins <- make_bins(
    genome = genome,
    width = width,
    by.chrom = FALSE,
    omit = omit,
    step = step
  )

  counts <- bin_overlaps(aln.pos = aln.gr, bins = genome.bins)

  return(list(total = length(aln.gr), counts = counts))
}

# Create a RangedSummarizedExperiment holding bins and bin counts from BAM files.
bin_counts <- function(bam.files, paired = FALSE, mid = FALSE, yield.size = 1E6,
                       mapq.cutoff = 0, filter.gr = NULL, omit = NULL, width, step = NULL,
                       genome, BPPARAM = SerialParam(), ...){

  file.names <- gsub(paste( ".*/(.*?)\\", ".bam", "(\\.\\w+)*(\\.\\w+)?$", sep = "" ),
                       "\\1", bam.files )

  if (length(width) > 1){
    stop("length(width) cannot be > 1..",
         call. = FALSE)
  }

  if (!bpisup(BPPARAM)){
    bpstart(BPPARAM)
    on.exit(bpstop(BPPARAM))
  }

  out <- bpmapply(
    FUN = counts_wrapper,
    bam.file = bam.files,
    MoreArgs = list(
      paired = paired,
      mid = mid,
      yield.size = yield.size,
      mapq.cutoff = mapq.cutoff,
      filter.gr = filter.gr,
      omit = omit,
      width = width,
      step = step,
      genome = genome,
      ...
    ),
    BPPARAM = BPPARAM,
    SIMPLIFY = FALSE
  )

  count.mat <- do.call(cbind, lapply(out, function(x) x$counts))
  colnames(count.mat) <- file.names

  col.dat <- data.frame(
    total = do.call(c, lapply(out, function(x) x$total))
  )
  rownames(col.dat) <- file.names

  genome.bins <- make_bins(
    genome = genome,
    width = width,
    by.chrom = FALSE,
    omit = omit,
    step = step
  )

  rse <- SummarizedExperiment(
    assays = SimpleList(counts = count.mat),
    colData = col.dat,
    rowRanges = genome.bins
  )
  return(rse)
}

# Calculate TSA-seq scores
# Supports two normalization approaches:
# 1) rel.input: log2((Ntsa/Ttsa) / (Ninput/Tinput))
# 2) rel.pulldown: log2(N' / mean(N')),
# where N' = (Ntsa * mean(Ninput)) / Ninput
calc_scores <- function(pulldown, input, norm.type = "rel.pulldown", samp.names){
  if (!norm.type %in% c("rel.pulldown", "rel.input")){
    stop("norm.type must be one of 'rel.pulldown' or 'rel.input'.")
  }

  pd.mat <- assays(pulldown)$counts
  in.mat <- assays(input)$counts
  colnames(pd.mat) <- samp.names
  colnames(in.mat) <- samp.names

  pd.tot <- colData(pulldown)$total
  in.tot <- colData(input)$total

  col.dat <- data.frame(
    pulldown.total = pd.tot,
    input.total = in.tot
  )
  rownames(col.dat) <- samp.names

  if (norm.type == "rel.pulldown"){
    raw.scores <- sapply(
      X = seq_len(ncol(pd.mat)),
      FUN = function(x){
        mean.input <- mean(in.mat[,x][in.mat[,x] != 0])
        n.prime <- ifelse(((pd.mat[,x] > 0) & (in.mat[,x] > 0)),
                          ((pd.mat[,x]) * mean.input) / (in.mat[,x]),
                          NA)
        mean.prime <- mean(n.prime, na.rm = TRUE)
        norm.sc <- ifelse(!is.na(n.prime), log2(n.prime / mean.prime), NA)
        return(norm.sc)
      }
    )
  } else if (norm.type == "rel.input"){
    raw.scores <- sapply(
      X = seq_len(ncol(pd.mat)),
      FUN = function(x){
        limit <- log2(((1 / pd.tot[x]) * 1E6) / ((1 / in.tot[x]) * 1E6))
        norm.sc <- log2((((pd.mat[,x] + 1) / pd.tot[x]) * 1E6) /
                        (((in.mat[,x] + 1) / in.tot[x]) * 1E6))
        norm.sc[norm.sc == limit] <- NA
        return(norm.sc)
      }
    )
  }

  colnames(raw.scores) <- samp.names

  rse <- SummarizedExperiment(
    assays = list(pulldown.counts = pd.mat,
                  input.counts = in.mat,
                  raw.scores = raw.scores),
    rowRanges = rowRanges(pulldown),
    colData = col.dat
  )

  return(rse)
}

# Smooth normalized TSA-seq scores
# Window-based smoothing by convolution
# By default, uses a Hanning window of length = 21
smooth_scores <- function(norm.se, smooth.type = "hanning",
                          win.len = 21, BPPARAM){
  if (!bpisup(BPPARAM)){
    bpstart(BPPARAM)
    on.exit(bpstop(BPPARAM))
  }

  split.se <- split(norm.se, seqnames(norm.se))

  sm.list <- bplapply(
    X = split.se,
    FUN = function(x){
      tmp.se <- x
      mat <- assays(tmp.se)$raw.scores
      sm.scores <- apply(
        X = mat,
        MARGIN = 2,
        FUN = function(y){
          if (length(unique(y)) == 1){
            sm <- rep(NA, length(y))
          } else{
            new.y = y
            new.y[is.na(new.y)] <- 0
            sm <- win_smooth(
              new.y,
              window.len = win.len,
              window = smooth.type
            )
          }
          return(sm)
        }
      )
      colnames(sm.scores) <- colnames(tmp.se)
      assays(tmp.se)$smooth.scores <- sm.scores
      return(tmp.se)
    },
    BPPARAM = BPPARAM
  )

  rse <- do.call(rbind, sm.list)

  assays(rse)$raw.scores[which(is.na(assays(rse)$smooth.scores))] <- NA

  return(rse)
}


# Window smoothing function
win_smooth <- function(x, window.len = 21, window = "hanning") {
  if (!is.numeric(x) || length(dim(x)) >= 2) {
    stop("smooth only accepts 1-dimensional numeric vectors.")
  }
  if (length(x) < window.len) {
    stop("Input vector needs to be bigger than window size.")
  }
  if (window.len < 3) {
    return(x)
  }
  if (!window %in% c("flat", "hanning", "hamming", "bartlett", "blackman",
                     "Hanning", "Hamming", "Bartlett", "Blackman")) {
    stop("Window must be one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'.")
  }

  s <- c(2 * x[1] - rev(x[1:window.len]),
         x,
         2 * x[length(x)] - rev(x[(length(x) - window.len + 1):(length(x) - 1)]))

  if (window == "flat") {
    w <- rep(1, window.len)
  } else if (window == "hanning" | window == "Hanning") {
    w <- 0.5 - 0.5 * cos(2 * pi * seq(0, window.len - 1) / (window.len - 1))
  } else if (window == "hamming" | window == "Hamming") {
    w <- 0.54 - 0.46 * cos(2 * pi * seq(0, window.len - 1) / (window.len - 1))
  } else if (window == "bartlett" | window == "Bartlett") {
    w <- 2 / (window.len - 1) *
      ((window.len - 1) / 2 - abs(seq(0, window.len - 1) - (window.len - 1) / 2))
  } else if (window == "blackman" | window == "Blackman") {
    w <- 0.42 - 0.5 *
      cos(2 * pi * seq(0, window.len - 1) / (window.len - 1)) +
      0.08 * cos(4 * pi * seq(0, window.len - 1) / (window.len - 1))
  }

  y <- as.vector(filter(s, (w / sum(w))))
  y <- y[(window.len + 1):(length(y) - (window.len - 1))]
  
  return(y)
}


# Calculate the respective quantile of each bin
calc_quantiles <- function(se, resolution = 0.01, data.type = "smooth"){
  if (!data.type %in% c("smooth", "raw")){
    stop("data.type must be one of 'smooth' or 'raw'.")
  }
  
  if (data.type == "smooth"){
    if (any(grepl("mean", names(assays(se))))){
      mat <- assays(se)$mean.smooth.scores
      new.name <- "mean.smooth.quantiles"
    } else{
      mat <- assays(se)$smooth.scores
      new.name <- "smooth.quantiles"
    }
  } else{
    if (any(grepl("mean", names(assays(se))))){
      mat <- assays(se)$mean.raw.scores
      new.name <- "mean.raw.quantiles"
    } else{
      mat <- assays(se)$raw.scores
      new.name <- "raw.quantiles"
    }
  }
  
  bins <- seq(0, 1, by = resolution)[-1]

  quantiles <- apply(
    X = mat,
    MARGIN = 2,
    FUN = function(x){
      thresholds <- quantile(x, probs = bins, na.rm = TRUE)
      bucket <- findInterval(x, thresholds, all.inside = TRUE)
      qq <- bins[bucket]
      return(qq)
    }
  )
  
  rse <- se
  assays(rse)[[new.name]] <- quantiles
  
  return(rse)  
}


# Combine replicates
# Averages biological replicates together. NA's are ignored.
combine_replicates <- function(se, cell.types){
  if (ncol(se) %% 2 != 0){
    stop("The number of columns in the SummarizedExperiment must be even.")
  }
  start.cols <- seq(1, (ncol(se) - 1), 2)
  if (any(grepl("smooth", names(assays(se))))){
    mat.list <- list(mean.raw.scores = assays(se)$raw.scores,
                     mean.smooth.scores = assays(se)$smooth.scores)
  } else{
    mat.list <- list(mean.raw.scores = assays(se)$raw.scores)
  }

  mean.scores <- lapply(
    X = mat.list,
    FUN = function(x){
      ms <- sapply(
        X = start.cols,
        FUN = function(y){
          tmp.mat <- x[,(y:(y + 1))]
          rowMeans(tmp.mat, na.rm = TRUE)
        }
      )
      colnames(ms) <- cell.types
      return(ms)
    }
  )

  rse <- SummarizedExperiment(
    assays = mean.scores,
    rowRanges = rowRanges(se)
  )
  return(rse)
}

# Filter quantiles
# Assumes a quantile assay is present in the SE object (calc_quantiles()).
# Done this way to facilitate easy quantile comparisons between datasets.
filter_quantiles <- function(se, threshold = 0.95, data.type = "smooth",
                             collapse = FALSE){
  if (!data.type %in% c("smooth", "raw")){
    stop("data.type must be one of 'smooth' or 'raw'.")
  }
  
  if (data.type == "smooth"){
    if (any(grepl("mean", names(assays(se))))){
      mat <- assays(se)$mean.smooth.quantiles
    } else{
      mat <- assays(se)$smooth.quantiles
    }
  } else{
    if (any(grepl("mean", names(assays(se))))){
      mat <- assays(se)$mean.raw.quantiles
    } else{
      mat <- assays(se)$raw.quantiles
    }
  }
  
  gr <- rowRanges(se)
  
  idx <- apply(
    X = mat,
    MARGIN = 2,
    FUN = function(x){
      which(x >= threshold)
    }
  )
  
  regions <- lapply(
    X = idx,
    FUN = function(x){
      gr2 <- gr[x]
      if (isTRUE(collapse)){
        gr2 <- GenomicRanges::reduce(gr2)
      }
      return(gr2)
    }
  )
  
  return(regions)
}


global_hmm <- function(df,
                       gr,
                       n.cond = 1,
                       n.samp = 2,
                       n.states = 2,
                       init.params = NULL,
                       rand.start = FALSE,
                       multi.start = FALSE,
                       n.start = 100,
                       init.iters = 25,
                       fit.only = FALSE){
  
  if(length(n.samp) != n.cond){
    stop("The stated number of conditions (n.cond) does not match the vector length of
         the number of samples in each condition (n.samp).",
         call. = FALSE)
  }
  
  if(ncol(df) != sum(n.samp)){
    stop("The sum of the number of samples in each condition (n.samp) does not equal the
         number of columns in the given dataframe.",
         call. = FALSE)
  }
  
  if(nrow(df) != length(gr)){
    stop("The number of rows in the dataframe must be equal to the number of
         ranges given.",
         call. = FALSE)
  }
  
  ranges <- lapply(
    X = seq_len(ncol(df)),
    FUN = function(x){
      tmp.ranges <- gr
      tmp.ranges$score <- df[,x]
      return(tmp.ranges)
    }
  )
  
  n.cols <- as.list(n.samp)
  max.len <- nrow(df)
  
  formulas <- list()
  start <- 1
  for (i in seq_along(n.cols)) {
    end <- start + n.cols[[i]] - 1
    formula.string <- paste(names(df)[start:end], collapse = " + ")
    formulas[[i]] <- as.formula(paste(formula.string, "~ 1"))
    start <- end + 1
  }
  
  family <- lapply(X = 1:n.cond,
                   FUN = function(x){
                     gaussian()
                   }
  )
  
  if(is.null(init.params)){
    mod <- depmix(response = formulas,
                  family = family,
                  data = df,
                  nstates = n.states,
                  ntimes = max.len)
  } else{
    mod <- depmix(response = formulas,
                  family = family,
                  data = df,
                  nstates = n.states,
                  ntimes = max.len,
                  respstart = init.params[[2]],
                  trstart = init.params[[1]])
  }
  
  if(!multi.start){
    fm <- fit(mod,
              emcontrol = em.control(maxit = 1000,
                                     tol = 1e-12,
                                     crit = "relative",
                                     random.start = rand.start,
                                     classification = "soft"),
              verbose = FALSE)
  } else{
    fm <- multistart2(mod,
                      nstart = n.start,
                      initIters = init.iters,
                      emcontrol = em.control(maxit = 1000,
                                             tol = 1e-12,
                                             crit = "relative",
                                             random.start = FALSE,
                                             classification = "soft"))
  }
  
  if(fit.only){
    return(fm)
  } else{
    states <- viterbi(fm, na.allow = TRUE)
    df$state <- states$state
    
    obj <- list(regions = ranges,
                fit = fm,
                states = states,
                assigned.df = df)
    return(obj)
  }
}

multistart2 <- function(object,
                        nstart = 10,
                        initIters = 10,
                        emcontrol = NULL,
                        verbose = FALSE,
                        ...) {
  llbest <- as.numeric(logLik(object))
  bestmodel <- object
  nfailed <- 0
  for (i in 1:nstart) {
    fmod <- try(fit(object, emcontrol = em.control(maxit = initIters)), silent = TRUE)
    if (inherits(fmod, "try-error")) {
      nfailed <- nfailed + 1
    } else {
      if (logLik(fmod) > llbest) {
        llbest <- logLik(fmod)
        bestmodel <- fmod
      }
    }
  }
  
  if (is.null(emcontrol)) {
    emcontrol <- em.control(random.start = FALSE)
  }
  bestmodel <- fit(bestmodel, emcontrol = emcontrol, verbose = verbose)
  if (nfailed > 0) {
    warning(
      nfailed,
      "out of",
      nstart,
      "attempts failed; result is based on ",
      nstart - nfailed,
      "starting values.\n"
    )
  }
  return(bestmodel)
}

# mod.obj is the output from global_hmm().
# base.response is the response that corresponds to the "reference"
# target of interest (e.g. SON).
# custom.order allows the specification of the reassigned state order.
parse_model <- function(mod.obj){
  
  ranges <- mod.obj[[1]]
  df <- mod.obj[[4]]
  colnames(df)["state"]
  fit <- mod.obj[[2]]
  ns <- nstates(fit)
  
  pars <- getpars(fit)
  tm.ind <- (ns + 1):(ns + (ns^2))
  rm.ind <- (ns + (ns^2) + 1):(length(pars))
  
  init.ps <- pars[1:ns]
  
  tm <- matrix(
    pars[tm.ind],
    byrow = TRUE,
    nrow = ns,
    ncol = ns
  )
  
  rm <- matrix(
    pars[rm.ind],
    byrow = TRUE,
    nrow = ns,
    ncol = length(rm.ind) / ns
  )
  
  rm.means <- rm[,seq(from = 1,
                      to = ncol(rm),
                      by = 2)]
  
  rm.sds <- rm[,seq(from = 2,
                    to = ncol(rm),
                    by = 2)]
  
  ord <- order(rm.means, decreasing = TRUE)
  
  new.states <- numeric(nrow(df))
  
  for ( i in seq_along(ord) ) {
    new.states[df$state == ord[i]] <- i
  }
  
  df$ordered.state <- new.states
  
  assigned.coords <- ranges[[1]]
  mcols(assigned.coords) <- NULL
  mcols(assigned.coords)$original.state <- df$state
  mcols(assigned.coords)$ordered.state <- new.states
  
  split.coords <- lapply(
    X = 1:length(ord),
    FUN = function(x){
      assigned.coords[mcols(assigned.coords)$ordered.state == x]
    }
  )
  
  names(df)[names(df) == "state"] <- "original.state"
  
  out <- list(coordinates = assigned.coords,
              split.coordinates = split.coords,
              df = df,
              fit = list(init.ps = init.ps,
                         tm = tm,
                         rm.means = rm.means,
                         rm.sds = rm.sds))
  
  return(out)
}

get_chromosome_seqs <- function(genome.obj){
  seqs <- lapply(X = seqlevels(genome.obj),
                 FUN = function(x){
                   string <- DNAString(genome.obj[[ x ]])
                   metadata(string) $name <- x
                   return(string)
                 }
  )
  names(seqs) <- seqnames(genome.obj)
  return(seqs)
}



# paths should be a character vector containing the paths of the respective wig files. Can be arbitrary in length.
# pattern should be string common to the files-of-interest. Expected to be the same across paths.
# filter.chroms allows the removal of certain chromosomes from the wig files before analysis.
# keep.chroms is an alternative to filter.chroms. It tells the function which chromosomes to retain.
# Only one of filter.chroms or keep.chroms can be given in a single function call. 
# filter.names allows the removal of particular files based on the given pattern.
# n.states defines the number of states in the model.

global_hmm <- function( paths, 
                        pattern, 
                        filter.chroms = NULL, 
                        keep.chroms = NULL, 
                        filter.names = NULL, 
                        n.states,
                        init.params = NULL,
                        rand.start = FALSE,
                        multi.start = FALSE,
                        n.start = 100,
                        init.iters = 25,
                        fit.only = FALSE,
                        data.only = FALSE ){
  
  pkgs <- c( "GenomeInfoDb", "rtracklayer", "GenomicRanges" )
  
  invisible(
    lapply( X = pkgs, 
            FUN = function(x){
              suppressPackageStartupMessages( library( x, character.only = TRUE ) )
              }
            )
    )
  
  if( !data.only ){
    suppressPackageStartupMessages( library( "depmixS4" ) ) 
    }
  
  if( !is.null( filter.chroms ) && !is.null( keep.chroms ) ){
    stop( "Only one of filter.chroms or keep.chroms can be given.",
          call. = FALSE )
    } 
  
  n <- length( paths )
  
  ins <- lapply( X = 1:n,
                 FUN = function(x){
                   lf <- list.files( path = paths[x],
                                     pattern = pattern,
                                     full.names = TRUE )
                   
                   if( !is.null( filter.names ) ){
                     if( length( filter.names ) > 1 ){
                       filter.names <- paste( filter.names, collapse = "|" )
                     }
                     lf <- lf[ !grepl( pattern= filter.names, x = lf ) ]
                   }
                   return( lf )
                  }
                 )
  
  n.cols <- lapply( X = 1:n,
                    FUN = function(x){
                      length( ins[[x]] )
                      }
                    )
  
  ins <- do.call( c, ins )
  
  wigs <- lapply( X = ins,
                  FUN = function(x){
                    w <- import(x)
                    
                    if( !is.null( filter.chroms ) ){
                      w <- w[ !seqnames(w) %in% filter.chroms ]
                    }
                    
                    if( !is.null( keep.chroms ) ){
                      w <- w[ seqnames(w) %in% keep.chroms ]
                    }
                    
                    seqlevelsStyle(w) <- "NCBI"
                    return(w)
                    }
                  )
  
  if( any( lapply( wigs, length ) == 0 ) ){
    stop( paste0( "Empty GRanges objects found at position(s) ",
                  paste( which( lapply( wigs, length ) == 0 ), collapse = ", " ), ". " ),
          "Enumerate unnecessary files with 'filter.names'.",
          call. = FALSE )
    }
  
  scores <- lapply( X = wigs,
                    FUN = function(x){
                      mcols(x)$score
                      }
                    )
  
  lens <- do.call( c, lapply( X=scores, FUN=function(x){ length(x) } ) )
  max.len <- max( lens )
  
  scores <- lapply(X=scores,
                   FUN=function(x){
                     if( length(x) < max.len ){
                       length(x) <- max.len
                       return(x)
                     } else{
                       x
                     }
                    }
                   )
  
  mat <- do.call( cbind, scores )
  df <- data.frame( mat )
  
  if( isTRUE( data.only ) ){
    return( df )
  }
  
  # set.seed(1)
  
  formulas <- list()
  start <- 1
  for (i in seq_along(n.cols)) {
    end <- start + n.cols[[i]] - 1
    formula.string <- paste( names( df )[ start:end ], collapse = " + " )
    formulas[[i]] <- as.formula( paste( formula.string, "~ 1" ) )
    start <- end + 1
  }
  
  family <- lapply( X = 1:n,
                    FUN = function(x){
                      gaussian()
                      }
                    )
  
  if( is.null( init.params ) ){
    mod <- depmix( response = formulas, 
                   family = family, 
                   data = df, 
                   nstates = n.states,
                   ntimes = max.len )
    } else{
      mod <- depmix( response = formulas, 
                     family = family, 
                     data = df, 
                     nstates = n.states,
                     ntimes = max.len,
                     respstart = init.params[[2]],
                     trstart = init.params[[1]] )
    }
  
  if( !multi.start ){
    fm <- fit( mod,
               emcontrol = em.control( maxit = 1000,
                                       tol = 1e-12,
                                       crit = "relative",
                                       random.start = rand.start,
                                       classification = "soft" ),
               verbose = FALSE )
  } else{
    fm <- multistart2( mod, 
                       nstart = n.start, 
                       initIters = init.iters, 
                       emcontrol = em.control( maxit = 1000,
                                               tol = 1e-12,
                                               crit = "relative",
                                               random.start = FALSE,
                                               classification = "soft" ) )
  }
  
  if( fit.only ){
    return( fm )
  } else{
    states <- viterbi( fm, na.allow = TRUE )
    df$state <- states$state
    
    obj <- list( regions = wigs,
                 fit = fm,
                 states = states,
                 assigned.df = df )
    
    return( obj )
  }
}

# Elbow function for calculating k-means over multiple k's.
elbow <- function( data, k.max, nstart, iter.max, algorithm = "Hartigan-Wong" ){
  wcss <- sapply( X = 1:k.max, 
                  FUN = function(k){
                    kmeans( x = data, 
                            centers = k, 
                            nstart = nstart, 
                            iter.max = iter.max,
                            algorithm = algorithm )$tot.withinss
                  }
  )
  
  return( wcss )
}


# mod.obj is the output from global_hmm().
# base.response is the response that corresponds to the "reference" target of interest (e.g. SON).
# custom.order allows the specification of the reassigned state order.
parse_model <- function( mod.obj, base.response = NULL, custom.order = NULL ){
  
  wigs <- mod.obj[[1]]
  df <- mod.obj[[4]]
  colnames(df)["state"]
  fit <- mod.obj[[2]]
  ns <- nstates(fit)
  
  pars <- getpars(fit)
  tm.ind <- ( ns + 1 ):( ns + ( ns^2 ) )
  rm.ind <- ( ns + ( ns^2 ) + 1 ):( length(pars) )
  
  init.ps <- pars[ 1:ns ]
  
  tm <- matrix( pars[ tm.ind ], 
                byrow = TRUE,
                nrow = ns, 
                ncol = ns )
  
  rm <- matrix( pars[ rm.ind ], 
                byrow = TRUE,
                nrow = ns, 
                ncol = length( rm.ind ) / ns )
  
  rm.means <- rm[, seq( from = 1, 
                        to = ncol(rm),
                        by = 2 ) ]
  
  rm.sds <- rm[, seq( from = 2, 
                      to = ncol(rm),
                      by = 2 ) ]
  
  if( is.null( custom.order ) ){
    ord <- order( rm.means[, base.response ], decreasing = TRUE )
  } else{
    ord <- custom.order
  }
  
  new.states <- numeric( nrow(df) )
  
  for ( i in seq_along(ord) ) {
    new.states[ df$state == ord[i] ] <- i
  }
  
  df$ordered.state <- new.states
  
  assigned.coords <- wigs[[1]]
  mcols( assigned.coords ) <- NULL
  mcols( assigned.coords )$original.state <- df$state
  mcols( assigned.coords )$ordered.state <- new.states
  
  split.coords <- lapply( X = 1:length(ord),
                          FUN = function(x){
                            assigned.coords[ mcols( assigned.coords )$ordered.state == x ]
                          }
  )
  
  names(df)[names(df) == "state"] <- "original.state"
  
  out <- list( coordinates = assigned.coords,
               split.coordinates = split.coords,
               df = df,
               fit = list( init.ps = init.ps,
                           tm = tm, 
                           rm.means = rm.means,
                           rm.sds = rm.sds ) )
  
  return( out )
  
}

# Lifted from factoextra
pca_lift <- function( prcomp.obj ){
  .get_pca_var_results <- function(var.coord){
    
    var.cor <- var.coord # correlation
    var.cos2 <- var.cor^2 # variable qualities 
    
    # variable contributions (in percent)
    # var.cos2*100/total Cos2 of the component
    comp.cos2 <- apply(var.cos2, 2, sum)
    contrib <- function(var.cos2, comp.cos2){var.cos2*100/comp.cos2}
    var.contrib <- t(apply(var.cos2,1, contrib, comp.cos2))
    
    colnames(var.coord) <- colnames(var.cor) <- colnames(var.cos2) <-
      colnames(var.contrib) <- paste0("Dim.", 1:ncol(var.coord)) 
    
    # Variable coord, cor, cos2 and contrib
    list(coord = var.coord, cor = var.cor, cos2 = var.cos2, contrib = var.contrib)
  }
  
  var_cor_func <- function(var.loadings, comp.sdev){var.loadings*comp.sdev}
  var.cor <- t(apply(prcomp.obj$rotation, 1, var_cor_func, prcomp.obj$sdev))
  var <- .get_pca_var_results(var.cor)
  
  return(var)
}


# Estimates transition probabilities and creates a transition matrix for HMMs
# Generates reasonable starting values for the Baum-Welch algorithm.
estimate_trans_mat <- function( x, n.states ) {

  # Enumerate the transitions present in the vector of states
  transpairs <- cbind( x[ -length( x ) ], x[ -1 ] )
  unique.states <- sort( unique( c( x, x + 1 ) ) )

  transprobs <- matrix( 0, nrow = n.states, ncol = n.states )

  # Count the number of times that state i transitions to state j
  # and create the transition matrix
  for( i in 1:n.states ) {
    state.i <- unique.states[ i ]
    state.i.count <- length( which( x == state.i ) )

    for(j in 1:n.states) {
      state.j <- unique.states[ j ]

      state.j.count <- suppressWarnings(
        length( which(
          rowSums( transpairs ) == ( state.i + state.j ) & x == state.i ) )
        )

      if( state.i.count == 0 ) {
        transprobs[ i,j ] <- 0
        } else {
          transprobs[ i,j ] <- state.j.count / state.i.count
        }
      }
    }

  colnames( transprobs ) <- paste( "State", 1:n.states )
  rownames( transprobs ) <- paste( "State", 1:n.states )

  return( transprobs )

}


# Lifted from https://atrebas.github.io/post/2019-06-08-lightweight-dendrograms/
dendro_data_k <- function(hc, k) {
  
  hcdata <- ggdendro::dendro_data(hc, type = "rectangle")
  seg <- hcdata$segments
  labclust <- cutree(hc, k)[hc$order]
  segclust <- rep(0L, nrow(seg))
  heights <- sort(hc$height, decreasing = TRUE)
  height <- mean(c(heights[k], heights[k - 1L]), na.rm = TRUE)
  
  for (i in 1:k) {
    xi <- hcdata$labels$x[labclust == i]
    idx1 <- seg$x    >= min(xi) & seg$x    <= max(xi)
    idx2 <- seg$xend >= min(xi) & seg$xend <= max(xi)
    idx3 <- seg$yend < height
    idx <- idx1 & idx2 & idx3
    segclust[idx] <- i
  }
  
  idx <- which(segclust == 0L)
  segclust[idx] <- segclust[idx + 1L]
  hcdata$segments$clust <- segclust
  hcdata$segments$line <- as.integer(segclust < 1L)
  hcdata$labels$clust <- labclust
  
  hcdata
}


# This is a very slight modification of the multistart() function/method in depmixS4.
# This function gives more control over the EM parameters in the final optimization.
multistart2 <- function(object, nstart = 10, initIters = 10, emcontrol = NULL, verbose = FALSE, ...) {
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

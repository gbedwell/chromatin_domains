# An adaptation of the first part of global_hmm.R
# Useful for preparing data for e.g. cluster analysis.
#
# paths should be a character vector containing the paths of the respective wig files. Can be arbitrary in length.
# pattern should be string common to the files-of-interest. Expected to be the same across paths.
# filter.chroms allows the removal of certain chromosomes from the wig files before analysis.
# keep.chroms is an alternative to filter.chroms. It tells the function which chromosomes to retain.
# Only one of filter.chroms or keep.chroms can be given in a single function call. 
# filter.names allows the removal of particular files based on the given pattern.

prepare_data <- function( paths, 
                          pattern, 
                          filter.chroms = NULL, 
                          keep.chroms = NULL, 
                          filter.names = NULL
                          ){
  
  pkgs <- c( "GenomeInfoDb", "rtracklayer", "GenomicRanges", "depmixS4" )
  
  invisible(
    lapply( X = pkgs, 
            FUN = function(x){
              suppressPackageStartupMessages( library( x, character.only = TRUE ) )
              }
            )
    )
  
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
  
  df <- data.frame( do.call( cbind, scores ) )
}





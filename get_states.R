get_states <- function( path, nstates, pattern = "_sm_20kb.wig", filter = NULL ){
  
  require( GenomicRanges )
  require( rtracklayer )
  
  files <- list.files( path = path, 
                       full.names = TRUE,
                       pattern = pattern )
  
  if( !is.null( filter ) ){
    filter <- 
    files <- files[ !grepl(pattern = filter, files ) ]
  }
  
  
  
  nfiles <- length( files )
  
  wigs <- lapply( X = files,
                  FUN = function(x){
                    rtracklayer::import( x )
                    }
                  )
  
  names( wigs ) <- gsub( paste0( pattern, "*" ), "", gsub( ".*/", "", files ) )
  
  fits <- lapply( X = wigs,
                  FUN = function(x){
                    w <- GenomicRanges::mcols(x)$score
                    states <- fit_hmm( data = w, nstates = nstates, return.states = TRUE )
                    GenomicRanges::mcols(x)$state <- states
                    return(x)
                    }
                  )
  
  return( fits )
}

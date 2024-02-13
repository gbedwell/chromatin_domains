prune_states <- function( gr.list, state ){
  
  require( GenomicRanges )
  
  nfiles <- length( gr.list )
  
  states <- lapply( X = gr.list,
                    FUN = function(x){
                      GenomicRanges::mcols(x)$state
                      }
                    )
  
  state.mat <- t( do.call( rbind, states ) )
  
  conserved <- which( rowSums( state.mat ) == nfiles * state )
  
  foi <- lapply( X = gr.list,
                 FUN = function(x){
                   nonreg <- x[ -conserved ]
                   mcols( nonreg )$region <- FALSE
                   
                   reg <- x[ conserved ]
                   mcols( reg )$region <- TRUE
                   c( reg, nonreg )
                   }
                 )
  return( foi )
}

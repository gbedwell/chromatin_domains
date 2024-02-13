# Called within get_states.R

fit_hmm <- function( data, nstates, return.states = TRUE ){
  
  require( HiddenMarkov )
  
  quants <- quantile( data, seq( 0.05, 0.95, length.out = nstates ) )
  
  Pi <- matrix( rep( 1/nstates, nstates^2 ),
                byrow = TRUE, nrow = nstates )
  
  delta <- rep( 1/nstates, nstates )
                                        
  obj <- dthmm(x = data,
               Pi = Pi,
               delta = delta,
               distn = "norm",
               pm <- list( mean = sort( quants[ 1:nstates ] ),
                           sd = rep( 1/nstates, nstates ) ) )
  
  obj <- BaumWelch( obj, control = bwcontrol(prt=FALSE) )
  
  if( return.states ){
    states <- Viterbi( obj )
    return( states )
  } else{
    return( obj )
  }
}

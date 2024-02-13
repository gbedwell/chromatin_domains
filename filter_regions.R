filter_regions <- function( gr.list, ntargets, common.target.name, cutoff = 0.28 ){
  
  targs <- gr.list[ grepl( pattern = target.name, x = names( gr.list ) ) ]
  
  if( length( targs ) != ntargets ){
    stop( "The number of identified target datasets does not match the stated number of targets.",
          "\n",
          "Target datasets should either have a common element in their name or different names separated by a pipe (|).",
          call. = FALSE )
  }
  
  rmv <- lapply( X = targs,
                 FUN = function(x){
                   gr <- x[ mcols( x )$region == TRUE ]
                   gr[ mcols( gr )$score <= cutoff ]
                   }
                 )
  
  rmv <- unlist( as( rmv, "GRangesList" ) )
  
  filt <- lapply( X = gr.list,
                  FUN = function(x){
                    gr <- x[ mcols( x )$region == TRUE ]
                    sub <- subtract( x=gr, y = rmv, minoverlap = 1L )
                    unlist( sub )
                    }
                  )
  
  return( filt )
}

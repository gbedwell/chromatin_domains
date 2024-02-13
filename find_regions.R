library( rtracklayer )
library( BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0 )
library( GenomicRanges )
library( HiddenMarkov )

source( "get_states.R" )
source( "fit_hmm.R" )
source( "prune_states.R" )
source( "filter_hmm.R" )

gs <- sum( seqlengths( BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0 ) )

laminB1.states <- get_states( path = "laminB1_data",
                              nstates = 2,
                              pattern = "_sm_20kb.wig" )

lads.pruned <- prune_states( gr.list = laminB1.states, state = 2 )

lads <- filter_regions( gr.list = lads.pruned,
                        ntargets = 2,
                        target.name = "K562",
                        cutoff = 0.28 )

lads <- reduce( unlist( as( lads, "GRangesList" ) ) )

sum( width( lads ) ) / gs


son.states <- get_states( path = "son_data",
                          nstates = 3,
                          pattern = "_sm_20kb.wig" )

spads.pruned <- prune_states( gr.list = son.states, state = 3 )

spads <- filter_regions( gr.list = spads.pruned,
                         ntargets = 2,
                         target.name = "K562",
                         cutoff = 0.28 )

lads <- reduce( unlist( as( spads, "GRangesList" ) ) )

sum( width( spads ) ) / gs

rtracklayer::export( object = lads, con = "../chm13v2.0_cLADs.bed", format =  "BED" )
rtracklayer::export( object = spads, con = "../chm13v2.0_cSPADs.bed", format =  "BED" )



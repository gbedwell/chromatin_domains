# source("funs.R")
load("slld_cluster.RData.gz")

# Model fitted using initial values for emission and transition matrices estimated from HC.
mod1 <- global_hmm( paths = c( "son_data", "laminB1_data", "laminB1_DamID_data" ), 
                    pattern = "combined.wig", 
                    filter.chroms = c( "chrM", "chrY", "MT", "Y" ),
                    n.states = 5,
                    init.params = list( tr = t(hc.trans.mat),
                                        em = do.call(rbind, hc.stats) ),
                    multi.start = FALSE,
                    rand.start = FALSE )

# Model fitted using randomized initial parameters.
mod2 <- global_hmm( paths = c( "son_data", "laminB1_data", "laminB1_DamID_data" ), 
                    pattern = "combined.wig", 
                    filter.chroms = c( "chrM", "chrY", "MT", "Y" ),
                    n.states = 5,
                    init.params = NULL,
                    multi.start = TRUE,
                    n.start = 100,
                    init.iters = 25 )

mod1.crit <- c( AIC = AIC(mod1[[2]]), BIC = BIC(mod1[[2]]),  LL = logLik(mod1[[2]]) )

mod2.crit <- c( AIC = AIC(mod2[[2]]), BIC = BIC(mod2[[2]]), LL = logLik(mod2[[2]]) )

save( list = ls( all.names = TRUE ), 
      file = "slld_model.RData.gz", 
      compression_level = 6 )
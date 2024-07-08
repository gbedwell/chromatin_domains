source("funs.R")

# Import the smoothed, normalized, and averaged data.
dat <- global_hmm( paths = c( "son_data", "laminB1_data", "laminB1_DamID_data" ),
                   pattern = "combined.wig",
                   filter.chroms = c( "chrM", "chrY", "MT", "Y" ),
                   n.states = NULL,
                   fit.only = FALSE,
                   data.only = TRUE )

# Do principle component analysis on the data.
pca <- prcomp( dat, scale. = TRUE )

# Get the contributions of each PC to the data variance.
prop.var <- summary(pca)$importance

# Define the threshold for which PCs to take for downstream analyses.
# This is the minimum cumulative percentage of the variance explained desired.
# This is not used currently.
# cutoff <- 0.90

# Get the index of the PC that matches the cutoff.
# Not currently used.
# dims <- sort( which( prop.var[3,] <= cutoff ), decreasing = TRUE )[1]
# pc.dat <- pca$x[,1:dims]

samp.rows <- sort( sample( 1:nrow(pc.dat), size = round( nrow(pc.dat) / 2 ) ) )

hc <- hclust( dist( dat[samp.rows,] ), method = "complete" )

nc <- 5

hc.states <- cutree( hc, k = nc )

hc.stats <- lapply( X = list(1:4, 5:8, 9:12),
                    FUN = function(x){
                      sapply( X = 1:nc,
                              FUN = function(y){
                                ind <- which( hc.states == y )
                                rbind(
                                  mean( rowMeans( dat[samp.rows,][,x][ind,] ) ),
                                  sqrt( var( rowMeans( dat[samp.rows,][,x][ind,] ) ) )
                                )
                              }
                      )
                    }
)

hc.state.order <- order(hc.stats[[1]][1,], decreasing = TRUE)

hc.stats <- lapply( X = hc.stats,
                    FUN = function(x){
                      x[,hc.state.order]
                    }
)

hc.states.new <- numeric( length(hc.states) )
for ( i in seq_along(hc.state.order) ) {
  hc.states.new[ hc.states == hc.state.order[i] ] <- i
}

hc.trans.mat <- estimate_trans_mat( x = hc.states.new, n.states = nc )

dend.dat <- dendro_data_k( hc = hc, k = nc )

# Calculate within-cluster sum of squares across k's.
# k.max <- 15
#
# wcss <- elbow( data = dat,
#                k.max = k.max,
#                nstart = 100,
#                iter.max = 1E6 )

save( list = ls( all.names = TRUE ),
      file = "slld_cluster.RData.gz",
      compression_level = 6 )

# save( dend.dat,
#       file = "dend_dat_7.RData.gz",
#       compression_level = 6 )

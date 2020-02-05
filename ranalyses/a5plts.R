load( 'a5.rda' )


library( treedater ) 
library( lubridate )
library( treestructure )
library( phydynR )

library( WVPlots )
library( ggplot2 )


# doubling 
pbdblsci = quantile( pbdbls, c(.025 , .975) )
dbldf <- data.frame( dbls = pbdbls )
pl = ShadedDensityCenter(
	frame = dbldf 
	,xvar = 'dbls'
	,boundaries = pbdblsci
	,title = 'Doubling time' 
	,linecolor = "black"
	,shading = "darkblue"
	,annotate_area = FALSE
)
pl <- pl + theme_minimal() + xlab('Days') + ylab('') + geom_vline( aes( xintercept = mean( tddbls )) ,col='red' ) 
ggsave( plot=pl, file = 'a5dbl.svg', width = 3.85, height = 3.4 )

svg('a5dbl2.svg',width = 3.85, height = 3.4 )
plot( density( pbdbls ), xlab = 'Days', ylab = '', main = 'Doubling time', bty='n', axes=FALSE)
axis(1)
abline( v= pbdblsci, col = 'red', lty = 3 )
abline( v = mean(tddbls), col = 'red', lty = 1 )
dev.off() 



#~ -----------------
#tmrca 
tmrca
pbtmrca
ci = quantile( pbtmrca, c(.025 , .975) )
tmrcadf <- data.frame( tmrca = pbtmrca )
pl = ShadedDensityCenter(
	frame = tmrcadf 
	,xvar = 'tmrca'
	,boundaries = ci
	,title = 'TMRCA' 
	,linecolor = "black"
	,shading = "darkblue"
	,annotate_area = FALSE
)
pl <- pl + theme_minimal() + xlab('') + ylab('') + geom_vline( aes( xintercept = tmrca) ,col='red' ) 
ggsave( plot=pl, file = 'a5tmrca.svg', width = 3.85, height = 3.4 )

svg('a5tmrca2.svg',width = 3.85, height = 3.4 )
plot( density( pbtmrca ), xlab = '', ylab = '', main = 'TMRCA', bty='n', axes=FALSE)
axis(1)
abline( v= ci, col = 'red', lty = 3 )
abline( v = tmrca, col = 'red', lty = 1 )
dev.off() 

invisible('
feb 3
- rate & likelihoods 
')

library( treedater ) 
library( lubridate )

tr <- read.tree( '../results/tree_raw.nwk' )
md <- read.table( '../data/metadata.tsv', header=TRUE, sep = '\t', stringsAsFactors=FALSE) 
md$tx <- decimal_date( ymd( md$date ))

# collect sample time info 
sts <- setNames( rep(NA, Ntip(tr)), tr$tip.label )
sts[ tr$tip.label %in% md$strain ] <- md$tx[ match( tr$tip.label, md$strain)]

RATES <- seq( .0002, .001, length = 50)

rate2td <- function( omega ){
	dater( 
		unroot(tr) # ML tree 
		, sts = sts # sample times
		, s = 29e3 # sequence length
		, omega0 = omega
		, numStartConditions= 0
		, meanRateLimits = c(omega , omega + 1e-6) # fixing the rate to omega 
		, searchRoot = 15
	)
}

tds  <- lapply ( RATES, rate2td )

tmrcas <- sapply( tds, '[[', 'timeOfMRCA' )
logliks <- sapply( tds, '[[', 'loglik' )

pldf <- data.frame( rate = RATES, tmrca = tmrcas, loglik = logliks )
library( ggplot2) 

p0 = qplot( x = rate , y = loglik, data = pldf , geom = 'line' ) + theme_minimal()
p1 = qplot( x = rate , y = tmrcas, data = pldf , geom = 'point' ) + theme_minimal()

ggsave( file='rateAndLikelihood-tmrca.png', plot = p1, width = 3.5, height = 3, unit='in', dpi = 300 )
ggsave( file='rateAndLikelihood-loglik.png', plot = p0, width = 3.5, height = 3, unit='in', dpi = 300 )


# also compare to lsd: 
# make lsd input 
fileConn<-file("sts.tab")
writeLines( as.character(length(sts)), fileConn)
close(fileConn)
stsdf = data.frame( taxon = names( sts) , time = unname( sts ) )
write.table( stsdf, file = 'sts.tab', col.names=FALSE , sep = ' ', quote=FALSE, row.names=FALSE, append=TRUE)

invisible("
lsd -i '../results/tree_raw.nwk' -s 29000 -d 'sts.tab' -r a
lsd 0.2 

rate 0.000523, tMRCA 2019.852

> date_decimal( 2019.852)
[1] 2019-11-07 23:31:12 UTC

"

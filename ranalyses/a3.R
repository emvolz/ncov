invisible('
feb 2
- note no missing sample times in this set 
- estimating root position
- estimating growth rates with exponential growth coalescent
- parametric bootstrap for CIs 
- parm exp growth coalescent
')

library( treedater ) 
library( lubridate )

tr <- read.tree( '../results/tree_raw.nwk' )
md <- read.table( '../data/metadata.tsv', header=TRUE, sep = '\t', stringsAsFactors=FALSE) 
md$tx <- decimal_date( ymd( md$date ))

# collect sample time info 
sts <- setNames( rep(NA, Ntip(tr)), tr$tip.label )
sts[ tr$tip.label %in% md$strain ] <- md$tx[ match( tr$tip.label, md$strain)]


# search root: 
td <- dater( 
	unroot(tr) # ML tree 
	, sts = sts # sample times
	, s = 29e3 # sequence length
	, omega0 = c( 0.0005, 0.00075, .001) # initial guess of clock rate 
	, numStartConditions= 0
	, meanRateLimits = c(.0004 , .0012)
	, searchRoot = 15
)

#~ pb <- parboot( td , overrideTemp=FALSE, ncpu = 8, nreps = 500) 
pb <- parboot( td , overrideTemp=FALSE, overrideSearchRoot=FALSE, ncpu = 8, nreps = 500) 

print( td )
print( pb )

print( date_decimal( td$timeOf ))
print( date_decimal( pb$timeOf ))

invisible('
> print( td )

Phylogenetic tree with 36 tips and 35 internal nodes.

Tip labels:
	Wuhan/IVDC-HB-01/2019, Wuhan/IVDC-HB-05/2019, Wuhan/IPBCAMS-WH-01/2019, Wuhan/IPBCAMS-WH-03/2019, Wuhan-Hu-1/2019, Wuhan/HBCDC-HB-01/2019, ...

Rooted; includes branch lengths.

 Time of common ancestor 
2019.92648059518 

 Time to common ancestor (before most recent sample) 
0.147289896622169 

 Weighted mean substitution rate (adjusted by branch lengths) 
0.000894427190999916 

 Unadjusted mean substitution rate 
0.000894427190999916 

 Clock model  
strict 

 Coefficient of variation of rates 
0 
> print( pb )
                           pseudo ML        2.5 %      97.5 %
Time of common ancestor 2.019926e+03 2.019804e+03 2.01995e+03
Mean substitution rate  8.944272e-04 5.562857e-04 1.43811e-03

> print( date_decimal( td$timeOf ))
[1] "2019-12-05 03:58:12 UTC"
> print( date_decimal( pb$timeOf ))
                     2.5%                     97.5% 
"2019-10-21 15:19:30 UTC" "2019-12-13 19:44:22 UTC" 


')


# sum of branch lengths( approx number of snps)
print( sum( tr$edge.length ) * 29e3 )
invisible( '
40.29
')

# plots: 
if (TRUE)
{
	png( 'treedaterTimeTree-nextstrain03feb2020.png', width = 480*3, height = 480*2) 
	par( mfrow = c( 1,2 ))
	plot( tr, main = 'Maximum likelihood tree' );  add.scale.bar(x = max(node.depth.edgelength(tr))-1/29e3,  y= 0)
	plot( td , main = 'Time tree') ; axisPhylo( root.time = td$timeOf , backward = FALSE)
	dev.off() 

	png('treedaterTMRCA-nextstrain03feb2020.png' , width = 480, height = .75*480 )
	plot( density( pb$tmrcas ), xlab = '', ylab = '', main = 'Estimated TMRCA (bootstrap distribution)', bty='n', axes=FALSE)
	axis(1)
	abline( v= pb$timeOf, col = 'red', lty = 3 )
	abline( v = td$timeOf, col = 'red', lty = 1 )
	dev.off() 
}


# root to tip 

png( 'treedater-rtt-nextstrain03feb2020.png', width = 2* 480, height = 1.7 * 480, pointsize=24) 
rootToTipRegressionPlot( td , show.tip.labels=F, pch = 20, cex = 1,bty='n' ) 
dev.off()
invisible('
Root-to-tip mean rate: 0.000698860851846057 
Root-to-tip p value: 0.0147443295365193 
Root-to-tip R squared (variance explained): 0.162602535784956 
') 



#~ -------------------------------
#~ Exponential growth coalescent 
library( phydynR) 

#~ phylodynamically redundant sets:
rs = list( 
c( 'Guangdong/20SF040/2020', 'Guangdong/20SF028/2020' , 'Guangdong/20SF174/2020') # family , identical 
#, c( 'HKU-SZ-002a_2020', 'HKU-SZ-005b_2020' ) # family but _not identical_
, c('France/IDF0373/2020', 'France/IDF0372/2020') #identical
, c('Foshan/20SF210/2020', 'Foshan/20SF211/2020' ) # guessing based on location, identical
#~ , c( 'Nonthaburi/74/2020', 'Nonthaburi/61/2020' ) #thailand  different events (rambaut)
, c( 'WHU01', 'WHU02' ) #identical, same time , guessing
, c('Wuhan/WIV04/2019', 'Wuhan/WIV06/2019' ) # identical,  guessing based on time and id
, c( 'Guangdong/20SF012/2020', 'Guangdong/20SF025/2020', 'Guangdong/20SF013/2020' )
)

remove_phylodynamically_redundant = rpr <- function(td){
	tr2 <- td
	class(tr2) <- 'phylo'
	for ( k in 1:length( rs ) ){
		m <- rs[[k]][ rs[[k]] %in% tr2$tip.label ]
		if ( length( m ) > 1){
			tr2 <- drop.tip( tr2, m[-1] )
		}
	}
	tr2
}
	
	
# simple sir model under exponential growth
gamma <- 1/8.4 

# model equations 
infections = c(i = ' parms$beta * s * i/(s+i ) ' )
removals = c( i = ' parms$gamma * i' )
other = c( s = '-parms$beta * s * i / (s + i )' 
 , r = 'parms$gamma * i' ) 

parms = list(
	beta =  .25* 365 # arbitrary initial condition
	, gamma = gamma * 1/8.4 
)
x0 = c( i = 1 , s = 1e9, r = 0 )

# make the simulator 
model <- build.demographic.process( births = infections 
  , nonDeme = other
  , deaths = removals 
  , parameterNames = names(parms)
  , sde = FALSE
  , rcpp = FALSE 
)

g0 <- rpr( td )
ssts = cbind ( i = rep(1, Ntip(g0)), rep(0, Ntip(g0)))
bdt0 = DatedTree( g0, sampleTimes = td$sts[g0$tip.label] , sampleStates = ssts )

of <- function(theta, t0 = td$timeOf, bdt = bdt0 ){
	logr <- theta[1]
	logi0 <- theta[2]
	r <- exp( logr )
	i0 <- exp( logi0 )
	.theta <- parms 
	beta <- parms$gamma + r 
	.theta$beta <- beta
	.x0 <- x0  
	.x0['i'] <- unname(i0)
	rv = colik( bdt, .theta, model, x0 = .x0, t0 = t0, res = 100 , forgiveAgtY = TRUE , AgtY_penalty = 0)
	print( c( r, i0, rv ))
	rv
}
#~ of( c( log(60), log(1) ) )

f0 <- optim( par = c(log(60), log(1)), fn = of , control = list( fnscale = -1) )
mle_r <- exp( f0$par[1] )
bootreplicate_to_growthrate <- function( pbtree ){
	g1 <- rpr( pbtree )
	bdt1 = DatedTree( g1, sampleTimes = td$sts[g1$tip.label] , sampleStates = ssts )
	f1 = optim( par = c(log(60), log(1)), fn = of, t0 = pbtree$timeOf, bdt = bdt1 , control = list( fnscale = -1) )
	exp( f1$par[1] )
}

pbgrs <- parallel::mclapply( pb$trees, bootreplicate_to_growthrate , mc.cores = 8)
pbgrs <- unlist( pbgrs ) 
pbr <- quantile( pbgrs, c(.025, .975 ) )
dbl <- ( 365 * log(2) / c( mle_r , rev(pbr) ) )
print( dbl)

png('phydynRDoublingTimes-a3-nextstrain03feb2020.png' , width = 480, height = .75*480 )
plot( density( ( 365 * log(2) / pbgrs ) ),ylab = '', main = 'Estimated doubling time (bootstrap distribution)', bty='n', axes=FALSE, xlab = 'Days')
axis(1)
abline( v= ( 365 * log(2) / pbr ) , col = 'red', lty = 3 )
abline( v =  ( 365 * log(2) / mle_r), col = 'red', lty = 1 )
dev.off() 

invisible('
 > print( dbl)
              97.5%      2.5% 
 7.190648  4.318249 16.994222 
7.2 days( 95% CI: 4.3-17.0)
')

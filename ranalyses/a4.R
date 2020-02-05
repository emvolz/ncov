invisible('
feb 4
- note no missing sample times in this set 
- estimating root position
- estimating growth rates with exponential growth coalescent
- parametric bootstrap for CIs 
- parm exp growth coalescent
')

NCPU <- 40

library( treedater ) 
library( lubridate )
library( treestructure )
library( phydynR )

tr <- read.tree( 'nextstrain_ncov_tree_04feb.nwk' )
md <- read.csv( '../data/metadata.csv', header=TRUE, stringsAsFactors=FALSE) 
md$tx <- decimal_date( ymd( md$date ))
tr <- unroot( multi2di( tr )  ) 
tr$edge.length <- pmax( 1e-6, tr$edge.length / 29e3 ) # note branch length is subst per genome ; we translate it to subst per site 

# collect sample time info 
sts <- setNames( rep(NA, Ntip(tr)), tr$tip.label )
sts[ tr$tip.label %in% md$strain ] <- md$tx[ match( tr$tip.label, md$strain)]

# search root: 
td <- dater( 
	tr # ML tree 
	, sts = sts # sample times
	, s = 29e3 # sequence length
	, omega0 = c( 0.00075, .001) # initial guess of clock rate 
	, numStartConditions= 0
	, meanRateLimits = c(.00075 , .00125)
	, searchRoot = 10
	, temp=FALSE
)



#~ pb <- parboot( td , overrideTemp=FALSE, ncpu = 8, nreps = 500) 
pb <- parboot( td , overrideTemp=FALSE, overrideSearchRoot=FALSE, ncpu = NCPU, nreps = 500) 

print( td )
print( pb )

print( date_decimal( td$timeOf ))
print( date_decimal( pb$timeOf ))

invisible('
 Time of common ancestor 
2019.9174394413 

 Time to common ancestor (before most recent sample) 
0.159063290944459 

 Weighted mean substitution rate (adjusted by branch lengths) 
0.000822949016875158 

 Unadjusted mean substitution rate 
0.000822949016875158 

 Clock model  
strict 

 Coefficient of variation of rates 
0 
> print( pb )
                           pseudo ML        2.5 %       97.5 %
Time of common ancestor 2.019917e+03 2.019860e+03 2.019944e+03
Mean substitution rate  8.229490e-04 7.012536e-04 9.657634e-04

 For more detailed output, $trees provides a list of each fit to each simulation 
> 
> print( date_decimal( td$timeOf ))
[1] "2019-12-01 20:46:10 UTC"
> print( date_decimal( pb$timeOf ))
                     2.5%                     97.5% 
"2019-11-10 19:33:17 UTC" "2019-12-11 10:53:21 UTC" 

')


# sum of branch lengths( approx number of snps)
print( sum( tr$edge.length ) * 29e3 )
invisible( '
55.1
')

# plots: 
if (TRUE)
{
	png( 'treedaterTimeTree-nextstrain04feb2020.png', width = 480*3, height = 480*2) 
	par( mfrow = c( 1,2 ))
	plot( tr, main = 'Maximum likelihood tree' );  add.scale.bar(x = max(node.depth.edgelength(tr))-1/29e3,  y= 0)
	plot( td , main = 'Time tree') ; axisPhylo( root.time = td$timeOf , backward = FALSE)
	dev.off() 

	png('treedaterTMRCA-nextstrain04feb2020.png' , width = 480, height = .75*480 )
	plot( density( pb$tmrcas ), xlab = '', ylab = '', main = 'Estimated TMRCA (bootstrap distribution)', bty='n', axes=FALSE)
	axis(1)
	abline( v= pb$timeOf, col = 'red', lty = 3 )
	abline( v = td$timeOf, col = 'red', lty = 1 )
	dev.off() 
}


# root to tip 

svg( 'treedater-rtt-nextstrain04feb2020.png', width = 3.9, height = 3.5) 
rootToTipRegressionPlot( td , show.tip.labels=F, pch = 20, cex = 1,bty='n' ) 
dev.off()
invisible('
Root-to-tip mean rate: 0.00066530112774178 
Root-to-tip p value: 0.00223352421981411
') 



#~ -------------------------------
#~ Exponential growth coalescent 
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
, c('England/01/2020', 'England/02/2020' )
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

g0 <- rpr( td )

library( treestructure )
ts0 <- trestruct( g0, minClade = 5 )
coargs <- lapply( levels( ts0$clustering ), function(k) {
	list(
		g  = keep.tip( g0, ts0$clusterSets[[k]] )
		, pbtrees = lapply( pb$trees, function(x){
			class(x) <- 'phylo'
			keep.tip( x,  ts0$clusterSets[[k]] )  
		})
	)
})

expCoalescent <- function( g0 , pbtrees )
{

	library( phydynR) 
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


	ssts = cbind ( i = rep(1, Ntip(g0)), rep(0, Ntip(g0)))
	bdt0 = DatedTree( g0, sampleTimes = td$sts[g0$tip.label] , sampleStates = ssts )

	of <- function(theta, t0 = NULL, bdt = bdt0 ){
		if ( is.null( t0 )){
			stop()
		}
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
	
	.t0 <-  max(sts[g0$tip.label])  - max( node.depth.edgelength( g0 )  )
	f0 <- optim( par = c(log(60), log(1)), fn = of , control = list( fnscale = -1) , t0 = .t0 )
	mle_r <- exp( f0$par[1] )
	bootreplicate_to_growthrate <- function( pbtree ){
		bdt1 = DatedTree( pbtree, sampleTimes = td$sts[g0$tip.label] , sampleStates = ssts )
		.t0 <-  max(sts[pbtree$tip.label])  - max( node.depth.edgelength( pbtree)  )
		f1 = optim( par = c(log(60), log(1)), fn = of, t0 = .t0, bdt = bdt1 , control = list( fnscale = -1) )
		exp( f1$par[1] )
	}

	pbgrs <- parallel::mclapply(pbtrees, bootreplicate_to_growthrate , mc.cores = NCPU)
	pbgrs <- unlist( pbgrs ) 
	pbr <- quantile( pbgrs, c(.025, .975 ) )
	dbl <- ( 365 * log(2) / c( mle_r , rev(pbr) ) )
	print( dbl)
	
	list( 
	pbgrs = pbgrs
	, dbl = dbl
	, dbls = ( 365 * log(2) / pbgrs )
	, mle_r = mle_r 
	)

}

res1 = expCoalescent( coargs[[1]][[1]], coargs[[1]][[2]] ) 
res2 = expCoalescent( coargs[[2]][[1]], coargs[[2]][[2]] ) 
res3 = expCoalescent( coargs[[3]][[1]], coargs[[3]][[2]] ) 

ress <- list( res1, res2, res3 ) 
dbls <- do.call( cbind, lapply( ress, '[[', 'dbls' ))
cdbls <- rowMeans(dbls )

print( quantile( cdbls, c(.025 , .5, .975) )  )
print(  apply( dbls, MAR=2, FUN=function(x) quantile(x, c(.025, .5, .975)) ) )


save.image( 'a4.rda' ) 

#~ png('phydynRDoublingTimes-a4-nextstrain04feb2020.png' , width = 480, height = .75*480 )
#~ plot( density( ( 365 * log(2) / pbgrs ) ),ylab = '', main = 'Estimated doubling time (bootstrap distribution)', bty='n', axes=FALSE, xlab = 'Days')
#~ axis(1)
#~ abline( v= ( 365 * log(2) / pbr ) , col = 'red', lty = 3 )
#~ abline( v =  ( 365 * log(2) / mle_r), col = 'red', lty = 1 )
#~ dev.off() 
#~ svg('phydynRDoublingTimes-a4-nextstrain04feb2020.png' , width = 3.9, height = 3.6)
#~ plot( density( ( 365 * log(2) / pbgrs ) ),ylab = '', main = 'Estimated doubling time (bootstrap distribution)', bty='n', axes=FALSE, xlab = 'Days')
#~ axis(1)
#~ abline( v= ( 365 * log(2) / pbr ) , col = 'red', lty = 3 )
#~ abline( v =  ( 365 * log(2) / mle_r), col = 'red', lty = 1 )
#~ dev.off() 

invisible('
              97.5%      2.5% 
10.143656  5.875914 17.150667 

')




if (FALSE)
{
	library( skygrowth ) 
	tr2 <- bdt0
	class( tr2 ) <- 'phylo' 
	sg = skygrowth.map( tr2, res = 2 )
	sg2 <- skygrowth.mcmc( tr2, res = 10, tau0 = 1e-4 )
}

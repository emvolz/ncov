invisible('
feb 4
- a5 average over different resolutions of polytomies 
- apply trestructure to each bootstrap tree and estimate growth rate for each 
- note no missing sample times in this set 
- estimating root position
- estimating growth rates with exponential growth coalescent
- parametric bootstrap for CIs 
- parm exp growth coalescent
')

NCPU <- 40
set.seed(20200204)
ntres <- 10 
nparboot <- 50 

library( treedater ) 
library( lubridate )
library( treestructure )
library( phydynR )


md <- read.csv( '../data/metadata.csv', header=TRUE, stringsAsFactors=FALSE) 
md$tx <- decimal_date( ymd( md$date ))
tr <- read.tree( 'nextstrain_ncov_tree_04feb.nwk' )
tr <- di2multi( tr, tol = 1e-5 )
trs <- lapply( 1:ntres, function(i) {
	tr = unroot( multi2di( tr )  )  
	tr$edge.length <- pmax( 1e-6, tr$edge.length / 29e3 ) # note branch length is subst per genome ; we translate it to subst per site 
	tr
})


# collect sample time info 
sts <- setNames( rep(NA, Ntip(tr)), tr$tip.label )
sts[ tr$tip.label %in% md$strain ] <- md$tx[ match( tr$tip.label, md$strain)]

tds <- lapply( trs, function(tr){
dater( 
	tr # ML tree 
	, sts = sts # sample times
	, s = 29e3 # sequence length
	, omega0 = c( 0.00075, .001) # initial guess of clock rate 
	, numStartConditions= 0
	, meanRateLimits = c(.0005 , .00125)
	, searchRoot = 10
	, temp=FALSE
)	
})

pbs <- lapply( tds, function(td){
 parboot( td , overrideTemp=FALSE, overrideSearchRoot=FALSE, ncpu = NCPU, nreps = nparboot ) 
})

pbtrees <- do.call( c, lapply( pbs, '[[', 'trees' ))


tmrca <- mean( sapply( tds, '[[', 'timeOfMRCA' ))
pbtmrca <- sapply( pbtrees, '[[', 'timeOfMRCA' )
pbtmrcaci <- quantile( pbtmrca , c(.025 , .975) )

print( date_decimal( tmrca ) )
print( date_decimal( pbtmrcaci ) )

invisible('
> print( date_decimal( tmrca ) )
[1] "2019-11-28 20:58:33 UTC"
> print( date_decimal( pbtmrcaci ) )
                     2.5%                     97.5% 
"2019-10-21 14:43:25 UTC" "2019-12-11 15:22:02 UTC" 
')

tdrate <- mean( sapply( tds, '[[', 'mean.rate'))
pbrates <- sapply( pbtrees, '[[', 'mean.rate')
print( tdrate )
print( quantile( pbrates, c(.025, .975 )))

invisible('
> print( tdrate )
[1] 0.0005898439
> print( quantile( pbrates, c(.025, .975 )))
        2.5%        97.5% 
0.0005676275 0.0009160608 
')

# sum of branch lengths( approx number of snps)
print( sum( tr$edge.length )  )
invisible( '
55.1
')


# root to tip 
if (FALSE){
svg( 'treedater-a5-rtt-nextstrain04feb2020.png', width = 3.9, height = 3.5) 
rootToTipRegressionPlot( tds[[1]] , show.tip.labels=F, pch = 20, cex = 1,bty='n' ) 
dev.off()
invisible('
Root-to-tip mean rate: 0.000728853432953303 
Root-to-tip p value: 0.00223649284876427 
Root-to-tip R squared (variance explained): 0.168898502544352 
') 
}


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
, c("Zhejiang/WZ-01/2020"  ,  "Zhejiang/WZ-02/2020"  ) 
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

alltrees <- lapply( c(tds, pbtrees ), rpr )


# model 
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



expCoalescent <- function( g0 )
{
	ssts = cbind ( i = rep(1, Ntip(g0)), rep(0, Ntip(g0)))
	bdt0 = DatedTree( g0, sampleTimes = sts[g0$tip.label] , sampleStates = ssts )

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
	mle_r 
}

trestruct_co <- function( tree ){
	ts0 <- trestruct( tree, minClade = 10 )
	nclust <- length( levels( ts0$clustering )  )
	if ( nclust  > 1 ){
		tres <- lapply( levels( ts0$clustering ), function(k) {
			 g0 = keep.tip( tree, ts0$clusterSets[[k]] )
			 class( g0 ) <- 'phylo'
			 g0 
		})
	} else{
		tres <- list( tree ) 
	}
	grs <- sapply( tres, expCoalescent )
	grs
}
#~ trestruct_co( alltrees[[23]] )

alltrees_grs <- parallel::mclapply( alltrees, trestruct_co, mc.cores = NCPU )

grs <- sapply( alltrees_grs, mean )
tdgrs <- grs[1:ntres]
pbgrs <- grs[-(1:ntres)]
tddbls <- ( 365 * log(2) /  tdgrs )
pbdbls <- ( 365 * log(2) /  pbgrs )

print( mean( tddbls ) )
print( quantile( pbdbls, c(.025 , .5, .975) )  )


save.image( 'a5.rda' ) 


invisible('
> print( mean( tddbls ) )
[1] 7.064405
> print( quantile( pbdbls, c(.025 , .5, .975) )  )
     2.5%       50%     97.5% 
 3.002793  8.717966 20.521134 
> 
')



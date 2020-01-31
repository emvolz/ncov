# adding samples in first half of February 

library( phydynR ) 
library( treedater ) 
library( skygrowth ) 
library ( lubridate )

set.seed(1111)

n_extra_samples <- 50
datefin = as.Date('2020-02-14')

# sample dates in nextstrain metadata
x = '2020-Jan-22
2020-Jan-22
2020-Jan-22
2020-Jan-23
2020-Jan-23
2020-Jan-14
2020-Jan-15
2020-Jan-15
2020-Jan-15
2020-Jan-17
2020-Jan-18
2020-Jan-22
2020-Jan-23
2020-Jan-22
2020-Jan-08
2020-Jan-13
2020-Jan-13
2020-Jan-13
2020-Jan-16
2020-Jan-16
2020-Jan-23
2020-Jan-19
2020-Jan-22
2020-Jan-23
2020-Jan-22
2020-Jan-21
2020-Jan-02
2020-Jan-02
2019-Dec-26
2019-Dec-30
2019-Dec-24
2019-Dec-30
2019-Dec-30
2019-Dec-30
2020-Jan-01
2019-Dec-30
2020-Jan-01
2019-Dec-30
2019-Dec-30
2019-Dec-30
2019-Dec-30
2019-Dec-30
2019-Dec-30
2020-Jan-16
2020-Jan-17'
x = strsplit( x , split = '\n' )
sts = lapply( x, function( xx ) as.Date(xx, format = '%Y-%b-%d'  ))[[1]]
sts <- decimal_date( sts )

# approx effective sample size because many current sequences were not randomly sampled (eg households)
n <- length( sts ) - 10 
sts <- sample( sts, size = n, replace=F) 

# additional samples 

sts <- c( sts, 
 seq( decimal_date(as.Date('2020-01-31')), decimal_date(datefin), length = n_extra_samples )
)

# tmrca from treedater
t0 = decimal_date( as.Date( '2019-11-28') )
# end simulation approx end of january: 
Tfin <-  decimal_date(datefin)

# simple sir model. values taken from 
#~ [ https://www.imperial.ac.uk/media/imperial-college/medicine/sph/ide/gida-fellowships/Imperial-2019-nCoV-transmissibility.pdf ]
R0  <- 2.2
# removal rate gamma is chosen to make the generation time of this SIR model match SARS 
gamma <- 1/8.4 
beta <- R0 * gamma 
r <- beta - gamma 
doubling_time <- log(2) / r 

# model equations 
infections = c(i = ' parms$beta * s * i/(s+i ) ' )
removals = c( i = ' parms$gamma * i' )
other = c( s = '-parms$beta * s * i / (s + i )' 
 , r = 'parms$gamma * i' ) 

parms = list(
	beta = beta * 365 
	, gamma = gamma * 365 
)
x0 = c( i = 1 , s = 1e9, r = 0 )

# make the simulator 
model <- build.demographic.process( births = infections 
  , nonDeme = other
  , deaths = removals 
  , parameterNames = names(parms)
  , sde = TRUE #stochastic differential equation 
  , rcpp = FALSE 
)

# compute conf interval for estimated growth rates 
sg_to_growthrate_ci <- function(sg) 
{
	deltat <- diff( range( sg$time ) )
	ci <- sg$ne_ci
	k <- nrow(ci)
	c( 
		log( ci[k,1] / ci[1,1] ) / deltat
		, log( ci[k,3] / ci[1,3] ) / deltat
	)
}

# function to run simulation replicate
#~ 1. simulate epidmeic trajectory
#~ 2. simulate coalescent tree conditional on avail sample times 
#~ 3. simulate genetic diversity tree 
#~ 4. estimate time tree
#~ 5. estimate growth rate 
process_replicate_baseline <- function(tfin = Tfin, clockrate = 4.418156e-04, shouldplot = FALSE, ncpu = 1, parboot = FALSE )
{
#~ tfin <- 1/12 + 2020
#~ clockrate = .0005
	
	# simulate epidemic trajectory
	fini <- -Inf 
	while( fini < 1e4 )
	{
		s <- model( parms, x0
		  , t0 = t0
		  , t1 = tfin
		  , res = 1000
		)
		fini = tail(s[[5]][, 'i'], 1)
	}
	
	
	# sample times
	traj <- as.data.frame( s[[5]] )
	finrecovered <- tail( traj$r, 1)
	n <- length( sts )

	# simulate tree
	ssts <-  cbind( rep( 1, n ), rep(0,n)) 
	cotree <- sim.co.tree.fgy(s, sts, ssts,finiteSizeCorrections =FALSE, substitutionRates = clockrate, sequenceLength = 29e3)
	dtre <- cotree[[2]]
	ctre <- cotree[[1]]
	
	# estimate time tree 
	dtre2 <- nj( cophenetic.phylo(dtre ))
	dtr = dater( dtre2, sts = cotree[[1]]$sampleTimes, s = 29e3 , omega0 = c(.0005,.001), numStartConditions = 0)
	
	if ( shouldplot )
	{
		png('sim1compareTrees.png', width = 480*4, height = 480*2) 
		par(mfrow = c(1,3))
		plot( ctre, main = 'Simulated coalescent tree', show.tip=F) ; axisPhylo(root.time = t0, backward=FALSE) 
		plot( dtre, main = 'Genetic distance tree' , show.tip=F); add.scale.bar()  
		plot( dtr, main = 'Estimated time tree', show.tip=F ) ; axisPhylo(root.time = t0, backward=FALSE) 
	}
	
	print( dtr$mean.rate )
	print( date_decimal( dtr$timeOf) ) 
	
	dtr$cotree = ctre
	dtr$dtree = dtre2 
	dtr$snps <- sum( 29e3*dtr$dtree$edge.length )
	
	# estimate growth rate using skygrowth 
	tr <- dtr
	class(tr) <- 'phylo' 
	# use res 2 so trajectory will be approx exponential 
	sg = skygrowth.map( tr, res = 2)
	growthrate <- median( na.omit( sg$growthrate ))
	dtr$sg <- sg
	dtr$growthrate <- growthrate 
	dtr$growthrate_CI <- sg_to_growthrate_ci( sg )
	
	# par bootstrap for rates and dates 
	if ( parboot ) {
		pb = parboot( dtr , ncpu = ncpu, overrideTemp=FALSE, nrep=50)
		dtr$rateci = pb$meanRate_CI
		dtr$dateci = pb$timeOfMRCA_CI
	}
	
	dtr
}

dtr = process_replicate_baseline() 
nreps <- 100
ncpu <- 8
result = list() 
for ( i in 1:nreps)
	result[[i]] <- process_replicate_baseline(ncpu = ncpu , parboot = TRUE) 



# dist of snps 
png(paste(sep='_', 'sim2', as.character(datefin), 'snpDistribution.png') )
snps <- sapply( result, '[[', 'snps' )
plot( density(snps ), xlab = '', ylab = '', main = 'Simulation distribution of SNPS', bty='n', axes=FALSE)
axis(1)
abline( v= 39.35938, col = 'red', lty = 1 )
dev.off()

# dist of tmrcas 
png(paste(sep='_', 'sim2', as.character(datefin), 'tmrcaEstimate.png'), width = 480*2.5)
tmrcas <- sapply( result, '[[', 'timeOfMRCA' )
tmrcalbs <- sapply( result, function(x) x$dateci[1])
tmrcaubs <- sapply( result, function(x) x$dateci[2])
i <- order( tmrcas )
library( Hmisc )
errbar( 1:nreps , tmrcas[i],  tmrcaubs[i] ,tmrcalbs[ i] , col = 'red', axes=FALSE, ylab = 'Simulation estimated TMRCA', ylim = c(2019,Tfin))
axis(2)
points( 1:nreps, tmrcas[i] )
abline( h= t0) 
dev.off() 

# dist of growth rates 
png( paste(sep='_', 'sim2', as.character(datefin),'grateDistribution.png') )
grates <-  sapply( result, '[[', 'growthrate' )
grates <- pmin( grates, 1000 ) #note a couple of outliers 
plot( density( grates ), xlab = '', ylab = '', main = 'Simulation estimated growth rates', bty='n', axes=FALSE, xlim = c(-10, 200) )
axis(1)
abline( v= r*365, col = 'red', lty = 1 )
dev.off()





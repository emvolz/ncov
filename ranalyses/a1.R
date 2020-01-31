library( treedater ) 
library( lubridate )

tr <- read.tree( 'tree_raw.nwk' )
md <- read.table( '../data/metadata.tsv', header=TRUE, sep = '\t', stringsAsFactors=FALSE) 
md$tx <- decimal_date( ymd( md$date ))

# collect sample time info 
sts <- setNames( rep(NA, Ntip(tr)), tr$tip.label )
sts[ tr$tip.label %in% md$strain ] <- md$tx[ match( tr$tip.label, md$strain)]

# some sample times are missing; here we compile a range of dates for these
est <- data.frame(  lower = rep(2020,Ntip(tr)), upper = rep(max(na.omit(sts)), Ntip(tr)) )
rownames(est) <- names(sts)
est <- est[ names(sts)[is.na(sts)] ,]

td <- dater( 
	tr # ML tree 
	, sts = sts # sample times
	, s = 29e3 # sequence length
	, estimateSampleTimes = est # ranges for tips with missing sample time
	, omega0 = c( 0.0005, 0.00075, .001) # initial guess of clock rate 
	, numStartConditions= 0
)

pb <- parboot( td , overrideTemp=FALSE, ncpu = 8, nreps = 500) 

print( td )
print( pb )
invisible('
> pb
                           pseudo ML        2.5 %       97.5 %
Time of common ancestor 2.019908e+03 2.019123e+03 2.019943e+03
Mean substitution rate  4.418156e-04 1.305774e-04 1.494906e-03
')

print( date_decimal( td$timeOf ))
print( date_decimal( pb$timeOf ))

invisible('
 "2019-11-28 08:42:49 UTC"
                     2.5%                     97.5% 
"2019-02-15 01:31:16 UTC" "2019-12-11 08:51:06 UTC" 
')

# sum of branch lengths( approx number of snps)
print( sum( tr$edge.length ) * 29e3 )
invisible( '
[1] 39.35938
')

png( 'treedaterTimeTree-nextstrain30jan2020.png', width = 480*3, height = 480*2) 
par( mfrow = c( 1,2 ))
plot( tr, main = 'Maximum likelihood tree' );  add.scale.bar(x = max(node.depth.edgelength(tr))-1/29e3,  y= 0)
plot( td , main = 'Time tree') ; axisPhylo( root.time = td$timeOf , backward = FALSE)
dev.off() 

png('treedaterTMRCA-nextstrain30jan2020.png' , width = 480, height = .75*480 )
plot( density( pb$tmrcas ), xlab = '', ylab = '', main = 'Estimated TMRCA (bootstrap distribution)', bty='n', axes=FALSE)
axis(1)
abline( v= pb$timeOf, col = 'red', lty = 3 )
abline( v = td$timeOf, col = 'red', lty = 1 )
dev.off() 

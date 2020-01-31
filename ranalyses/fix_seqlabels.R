library( seqinr )

d <- read.fasta( 'sequences.fasta') 
md <- read.table( 'metadata.tsv', sep = '\t', stringsAsFactors=FALSE, header=TRUE)
md[ md == '?' ] <- NA 
write.fasta( d, names(d), file = 'sequences.fasta.bak')


sids <- names(d)
.sids <- sapply(  strsplit( sids, '\\|'), function(x) head(x,1))
.sids <- sapply(  strsplit( .sids, '\\.'), function(x) head(x,1))
.sids <- gsub( 'BetaCoV/', '', .sids)

write.fasta( d, .sids, file = 'sequences.fasta' )

library( seqinr )

d <- read.fasta( '../data/sequences.fasta') 
md <- read.table( '../data/metadata.tsv', sep = '\t', stringsAsFactors=FALSE, header=TRUE)
md[ md == '?' ] <- NA 
write.fasta( d, names(d), file = '../data/sequences.fasta.2')

sids <- names(d)
sid2strain <- setNames(sids, sids )

# gisaid sequences have |
gisaid_id <- sapply(  strsplit( sids, '\\|'), function(x) tail(x,1))
i <- which ( gisaid_id %in% md$gisaid_epi_isl )
sid2strain[ i] <- md$strain[ match( gisaid_id[i], md$gisaid_epi_isl ) ] 

# genbank 
l <- which ( sids %in%  md$genbank_accession)
sid2strain[ l] <- md$strain[ match( sids[l], md$genbank_accession ) ] 

# stragglers ?
.sids <- sapply(  strsplit( sids, '\\|'), function(x) head(x,1))
.sids <- gsub( 'BetaCoV/', '', .sids)
.sids <- gsub( 'BetaCov/', '', .sids)
k <- setdiff( which( (.sids %in% md$strain) ) , c(i, l ))
sid2strain[ k] <- md$strain[k ]


# remove missing 
j <- which( !(sid2strain %in% md$strain ) )
sid2strain <- sid2strain[-j]
dd <- d[-j] 
names(dd) <- sid2strain[names(dd)] 
write.fasta( dd, names(dd),  file = '../data/sequences.fasta' )

#~ setdiff( md$strain, sid2strain )

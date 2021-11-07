FixNames <- function(reduced, cid) {
	txids <- get_txids(phylota = reduced, cid = cid, rnk = 'species')
	scientific_names <- get_tx_slot(phylota = reduced, txid = txids, slt_nm = 'scnm')
# clean the names
	scientific_names <- gsub('\\.', '', scientific_names)
	scientific_names <- gsub('\\s+', '_', scientific_names)
	sids <- reduced@clstrs[[cid]]@sids
	
}

RunPhylotaR <- function(wd, txid, ncbi_dr) {
	phylotaR::setup(wd = wd, txid = txid, ncbi_dr = ncbi_dr)
	phylotaR::run(wd = wd)
	return(phylotaR::read_phylota(wd))
}
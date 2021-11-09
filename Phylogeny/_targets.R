library(targets)
source("functions.R")
tar_option_set(packages=c("phylotaR", "ape"))

list(
 tar_target(wd, file.path(getwd(), 'sponges')),
 tar_target(ncbi_dr, '/usr/local/bin/'),
 tar_target(txid, 6040), # Porifera
 tar_target(all_clusters, RunPhylotaR(wd=wd, txid=txid, ncbi_dr=ncbi_dr)),
 tar_target(cids, all_clusters@cids),
 tar_target(n_taxa, phylotaR::get_ntaxa(phylota = all_clusters, cid = cids)),
 tar_target(keep, cids[n_taxa >= 10]),
 tar_target(selected, phylotaR::drop_clstrs(phylota = all_clusters, cid = keep)),
 tar_target(reduced, phylotaR::drop_by_rank(phylota = selected, rnk = 'species', n = 1)),
 tar_target(save_genes, SaveGenes(reduced)),
 tar_target(process_genes, ProcessSequencesByGeneSingle(save_genes)),
 tar_target(gaps_removed, RemoveGappy(process_genes)),
 tar_target(concatenate_all, ConcatenateAll(gaps_removed))
)
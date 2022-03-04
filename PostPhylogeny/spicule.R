library(ape)
library(geiger)
library(phytools)
library(datelife) #github

phy <- ape::read.tree("~/Downloads/raxml.tre")[[1]]
chars <- read.csv("/Users/bomeara/Documents/MyDocuments/GitClones/SpiculeComplexity/Matrix/CharacterMatrix.csv")
charSpecies <- gsub(" ", "_", chars$Species)
charSpecies[which(charSpecies %in% phy$tip.label)]

getGenus <- function(x) {
  genera <- rep(NA, length(x)) 
    for(i in seq_along(x)) {
      genera[i] <- strsplit(x[i], "_")[[1]][1]
    }
  
  return(genera)
}

charGenera <- getGenus(charSpecies)
charGenera[which(charGenera %in% getGenus(phy$tip.label))]

char_cleaned <- chars[,-1]
rownames(char_cleaned) <- gsub(" ", "_", chars$Species)
for (col_index in sequence(ncol(char_cleaned))) {
  for (row_index in sequence(nrow(char_cleaned))) {
    char_cleaned[row_index, col_index] <- gsub("\\?", "", char_cleaned[row_index, col_index])
  }
}

for (col_index in sequence(ncol(char_cleaned))) {
  char_cleaned[,col_index] <- as.numeric(char_cleaned[,col_index])
}
present <- t(t(rowSums(char_cleaned)))
pruned <- geiger::treedata(phy, present, sort=TRUE)
phytools::contMap(ape::compute.brlen(pruned$phy), pruned$data[,1])


phyGenus <- phy
phyGenus$tip.label <- getGenus(phy$tip.label)
phyGenus <- drop.tip(phyGenus, which(!duplicated(phyGenus$tip.label)))
char_cleaned_genus <- char_cleaned
genera <- getGenus(rownames(char_cleaned_genus))
char_cleaned_genus <- char_cleaned_genus[!duplicated(genera),]

rownames(char_cleaned_genus) <- genera[!duplicated(genera)]
present_genus <- t(t(rowSums(char_cleaned_genus)))

pruned_genus <- geiger::treedata(phyGenus, present_genus, sort=TRUE)

phytools::contMap(ape::compute.brlen(pruned_genus$phy), pruned_genus$data[,1], type="fan")

#sponge_datelife <- datelife_search(gsub("_", " ", charSpecies), summary_format="newick_all")

ott_taxa <- tnrs_match_names(charSpecies, context = "Animals")
ott_tree <- tol_induced_subtree(ott_id(ott_taxa)[is_in_tree(ott_id(ott_taxa))], label_format="name")
missing_in_ott_tree <- charSpecies[!charSpecies %in% ott_tree$tip.label]

ott_taxa2 <- tnrs_match_names(c(charSpecies, getGenus(missing_in_ott_tree)), context = "Animals")

ott_tree2 <- tol_induced_subtree(ott_id(ott_taxa2)[is_in_tree(ott_id(ott_taxa2))], label_format="name")

ott_taxonomy <- taxonomy_subtree(ott_id=ott_id(tnrs_match_names("Porifera")), label_format="name", output_format="newick")
ott_taxonomy_tree <- ape::read.tree(text=ott_taxonomy)

missing_in_raxml_tree <- charSpecies[!charSpecies %in% phy$tip.label]
missing_in_both_trees <- intersect(missing_in_raxml_tree, missing_in_ott_tree)

potential_missing_genera <- getGenus(missing_in_both_trees)
genus_missing_in_ott_tree <- potential_missing_genera[!potential_missing_genera %in% getGenus(ott_tree2$tip.label)]
genus_missing_in_raxml_tree <- potential_missing_genera[!potential_missing_genera %in% getGenus(phy$tip.label)]
genera_missing_in_both_trees <- intersect(genus_missing_in_ott_tree, genus_missing_in_raxml_tree)
genera_missing_in_taxonomy_tree_too <- genera_missing_in_both_trees[! genera_missing_in_both_trees %in% getGenus(ott_taxonomy_tree$tip.label)]

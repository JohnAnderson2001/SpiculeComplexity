library(phytools)
library(geiger)
library(rotl)
library(phangorn)
library(paleotree)

getGenus <- function(x) {
  genera <- rep(NA, length(x)) 
  for(i in seq_along(x)) {
    genera[i] <- strsplit(x[i], "_")[[1]][1]
  }
  
  return(genera)
}

badLabels <- function(x) {
  x[is.na(x)] <- ""
  badones <- which(nchar(x)==0)
  return(badones)
}


complexity <- read.csv("../Matrix/Complexity.csv")
complexity$Species <- gsub(" ", "_", complexity$Species)
charSpecies <- complexity$Species
charGenus <- getGenus(charSpecies)

phy <- ape::read.tree("~/Downloads/raxml.tre")[[1]]

paleodb_phy <- paleotree::makePBDBtaxonTree(paleotree::getCladeTaxaPBDB("Porifera"), rankTaxon = "species")

ott_taxa <- tnrs_match_names(charSpecies, context = "Animals")
ott_tree <- tol_induced_subtree(ott_id(ott_taxa)[is_in_tree(ott_id(ott_taxa))], label_format="name")
missing_in_ott_tree <- charSpecies[!charSpecies %in% ott_tree$tip.label]
genera_missing_in_ott_tree <- charGenus[!charGenus %in% getGenus(ott_tree$tip.label)]

ott_taxa2 <- tnrs_match_names(c(charSpecies, getGenus(missing_in_ott_tree)), context = "Animals")

ott_tree2 <- tol_induced_subtree(ott_id(ott_taxa2)[is_in_tree(ott_id(ott_taxa2))], label_format="name")

ott_taxonomy <- taxonomy_subtree(ott_id=ott_id(tnrs_match_names("Porifera")), label_format="name", output_format="newick")
ott_taxonomy_tree <- ape::read.tree(text=ott_taxonomy)

missing_in_raxml_tree <- charSpecies[!charSpecies %in% phy$tip.label]
genera_missing_in_raxml_tree <- charGenus[!charGenus %in% getGenus(phy$tip.label)]
missing_in_both_trees <- intersect(missing_in_raxml_tree, missing_in_ott_tree)
genera_missing_in_both_trees <- intersect(genera_missing_in_raxml_tree, genera_missing_in_ott_tree)



matching_species <- complexity$Species[complexity$Species %in% phy$tip.label]
#unmatching_species <- complexity$Species[!complexity$Species %in% phy$tip.label]
#genera_of_unmatching_species <- getGenus(unmatching_species)

for (i in seq_along(complexity$Species)) {
  if(!(complexity$Species[i] %in% matching_species)) {
    complexity$Species[i] <- getGenus(complexity$Species[i])
  }
}

complexity_vector <- complexity$SpiculeTypes
names(complexity_vector) <- complexity$Species

for(i in seq_along(phy$tip.label)) {
  if(!(phy$tip.label[i] %in% matching_species)) {
    if(!any(phy$tip.label==getGenus(phy$tip.label[i]))) {
      phy$tip.label[i] <- getGenus(phy$tip.label[i])
    }
  }
}

cleaned <- geiger::treedata(phy, complexity_vector)
cleaned$data <- cleaned$data[!duplicated(rownames(cleaned$data)),]

phytools::contMap(ape::compute.brlen(cleaned$phy), cleaned$data)

ott_taxonomy_tree_no_duplicates <- ott_taxonomy_tree
ott_taxonomy_tree_no_duplicates <- ape::drop.tip(ott_taxonomy_tree_no_duplicates, tip=sequence(length(ott_taxonomy_tree_no_duplicates$tip.label))[duplicated(ott_taxonomy_tree_no_duplicates$tip.label)])
badones <- badLabels(ott_taxonomy_tree_no_duplicates$tip.label)
if(length(badones)>0) {
  ott_taxonomy_tree_no_duplicates <- ape::drop.tip(ott_taxonomy_tree_no_duplicates, tip=badones)
}

phy_no_duplicates <- phy
phy_no_duplicates <- ape::drop.tip(phy_no_duplicates, tip=sequence(length(phy_no_duplicates$tip.label))[duplicated(ott_taxonomy_tree_no_duplicates$tip.label)])
badones <- badLabels(phy_no_duplicates$tip.label)
if(length(badones)>0) {
  phy_no_duplicates <- ape::drop.tip(phy_no_duplicates, tip=badones)
}

super <- phangorn::superTree(c(phy_no_duplicates, ott_taxonomy_tree_no_duplicates), trace=4)


cleaned_super <- geiger::treedata(super, complexity_vector)

phytools::contMap(ape::compute.brlen(cleaned_super$phy), cleaned_super$data)


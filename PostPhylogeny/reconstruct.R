library(phytools)
library(geiger)

getGenus <- function(x) {
  genera <- rep(NA, length(x)) 
  for(i in seq_along(x)) {
    genera[i] <- strsplit(x[i], "_")[[1]][1]
  }
  
  return(genera)
}

phy <- ape::read.tree("~/Downloads/raxml.tre")[[1]]

complexity <- read.csv("../Matrix/Complexity.csv")
complexity$Species <- gsub(" ", "_", complexity$Species)

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

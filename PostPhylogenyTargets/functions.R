

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

getComplexity <- function() {
	complexity <- read.csv("../Matrix/Complexity.csv")
	complexity$Species <- gsub(" ", "_", complexity$Species)
	return(complexity)
}

getRaxmlTree <- function() {
	phy <- ape::read.tree("raxml.tre")[[1]]
	return(phy)
}


getPaleoDBTree <- function() {
	paleodb_phy <- paleotree::makePBDBtaxonTree(paleotree::getCladeTaxaPBDB("Porifera"), rankTaxon = "species")
	return(paleodb_phy)
}


getCharSpecies <- function(complexity) {
	charSpecies <- complexity$Species
	return(charSpecies)
}

getCharGenus <- function(complexity) {
	charSpecies <- getCharSpecies(complexity)
	charGenus <- getGenus(charSpecies)

	return(charGenus)
}

getOTTTrees <- function(charSpecies, charGenus) {
	ott_taxa <- tnrs_match_names(charSpecies, context = "Animals")
	ott_tree <- tol_induced_subtree(ott_id(ott_taxa)[is_in_tree(ott_id(ott_taxa))], label_format="name")
	missing_in_ott_tree <- charSpecies[!charSpecies %in% ott_tree$tip.label]
	genera_missing_in_ott_tree <- charGenus[!charGenus %in% getGenus(ott_tree$tip.label)]

	ott_taxa2 <- tnrs_match_names(c(charSpecies, getGenus(missing_in_ott_tree)), context = "Animals")

	ott_tree2 <- tol_induced_subtree(ott_id(ott_taxa2)[is_in_tree(ott_id(ott_taxa2))], label_format="name")

	ott_taxonomy <- taxonomy_subtree(ott_id=ott_id(tnrs_match_names("Porifera")), label_format="name", output_format="newick")
	ott_taxonomy_tree <- ape::read.tree(text=ott_taxonomy)
	
	
	ott_taxonomy_tree_no_duplicates <- ott_taxonomy_tree
	ott_taxonomy_tree_no_duplicates <- ape::drop.tip(ott_taxonomy_tree_no_duplicates, tip=sequence(length(ott_taxonomy_tree_no_duplicates$tip.label))[duplicated(ott_taxonomy_tree_no_duplicates$tip.label)])
	badones <- badLabels(ott_taxonomy_tree_no_duplicates$tip.label)
	if(length(badones)>0) {
	ott_taxonomy_tree_no_duplicates <- ape::drop.tip(ott_taxonomy_tree_no_duplicates, tip=badones)
	}
	
	return(list(missing_in_ott_tree=missing_in_ott_tree, genera_missing_in_ott_tree=genera_missing_in_ott_tree, ott_tree=ott_tree, ott_tree2=ott_tree2, ott_taxonomy_tree=ott_taxonomy_tree, ott_taxonomy_tree_no_duplicates=ott_taxonomy_tree_no_duplicates))
}

getCleanedRaxmlTreeData <- function(charSpecies, raxml_phy, charGenus, missing_in_ott_tree, complexity, genera_missing_in_ott_tree) {
	phy <- raxml_phy
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
	return(cleaned)
}

plotComplexity <- function(cleaned_raxml) {
	cleaned <- cleaned_raxml
	phytools::contMap(ape::compute.brlen(cleaned$phy), cleaned$data)

}


# phy_no_duplicates <- phy
# phy_no_duplicates <- ape::drop.tip(phy_no_duplicates, tip=sequence(length(phy_no_duplicates$tip.label))[duplicated(ott_taxonomy_tree_no_duplicates$tip.label)])
# badones <- badLabels(phy_no_duplicates$tip.label)
# if(length(badones)>0) {
#   phy_no_duplicates <- ape::drop.tip(phy_no_duplicates, tip=badones)
# }

#super <- phangorn::superTree(c(phy_no_duplicates, ott_taxonomy_tree_no_duplicates), trace=4)


#cleaned_super <- geiger::treedata(super, complexity_vector)

#phytools::contMap(ape::compute.brlen(cleaned_super$phy), cleaned_super$data)

MakeConstraints <- function(fossil_phy) {
	constraints <- data.frame()
	unique_nodes <- unique(fossil_phy$edge[,1])
	if(is.null(fossil_phy$node.label)) {
		fossil_phy$node.label <- fossil_phy$edge[,1]
	}
	for (i in sequence(length(unique_nodes))) {
		descendant_numbers <- phangorn::Descendants(fossil_phy, unique_nodes[i], type="tips")[[1]]
		descendant_names <- fossil_phy$tip.label[descendant_numbers]
		non_descendant_names <- fossil_phy$tip.label[!(sequence(ape::Ntip(fossil_phy)) %in% descendant_numbers)]
		constraints <- rbind(constraints, data.frame(node=unique_nodes[1], name=fossil_phy$node.label[i], descendants=paste(descendant_names, collapse=" "), nondescendants=paste(non_descendant_names, collapse=" ")))
	}	
	constraints <- constraints[nchar(constraints$nondescendants)>0,]
	return(constraints)
}

PrintConstraints <- function(constraints, output_file="constraints.nex") {
	cat("begin mrbayes;\n", file=output_file, append=FALSE)

	for (i in sequence(nrow(constraints))) {
		cat("\nconstraint backbone partial = ", constraints$descendants[i], " : ", constraints$nondescendants[i], file=output_file, append=TRUE)	
	}
	cat("\nend;\n", file=output_file, append=TRUE)
}

GenerateTipAges <- function(paleodb_ages) {
	species <- subset(paleodb_ages, taxon_rank=="species")
	return(species[, c("taxon_name", "lastapp_min_ma")])
}

PrintTipAges <- function(species_ages) {
	cat("begin mrbayes;\ncalibrate ", file="tip_ages.nex", append=FALSE)
	for (i in sequence(nrow(species_ages))) {
		cat(gsub(" ", "_", species_ages$taxon_name[i]), " = Fixed(", species_ages$lastapp_min_ma[i], ")\n", file="tip_ages.nex", append=TRUE)	
	}
	cat("\n;	prset brlenspr=clock:uniform;
	prset clockvarpr=igr;	
	prset igrvarpr=exp(37.12);
	prset clockratepr = lognorm(-7.08069,2.458582);
	calibrate root=offsetexp(315,0.01234568);
	calibrate holometabola_with_fossils=offsetexp(302,0.0106383);
	prset topologypr=constraints(root, holometabola_with_fossils);
	prset nodeagepr = calibrated;
	end;", file="tip_ages.nex", append=TRUE)

}

PrintMrBayesBatchFile <- function() {
	cat("begin mrbayes;
	execute combined_for_mb.nex;
	execute constraints.nex;
	execute tip_ages.nex;
	mcmcp temp=0.1 nchain=4 samplefreq=1000 printfr=100 nruns=2;
	mcmcp filename=spongeDating;

	mcmc ngen=20000000;	
	end;", file="spongeDating.nex", append=FALSE)
}

RunMrBayes <- function(...) {
	system("mb spongeDating.nex")
}




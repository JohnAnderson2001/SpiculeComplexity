library(targets)
source("functions.R")


tar_option_set(packages=c("phytools", "geiger", "rotl", "phangorn", "paleotree"))



list(
 tar_target(complexity, getComplexity()),
 tar_target(raxml_tree, getRaxmlTree()),
 tar_target(paleodb_tree, getPaleoDBTree()),
 tar_target(charSpecies, getCharSpecies(complexity)),
 tar_target(charGenus, getCharGenus(complexity)),
 tar_target(ott_trees, getOTTTrees(charSpecies, charGenus)),
 tar_target(raxml_treedata, getCleanedRaxmlTreeData(charSpecies=charSpecies, raxml_phy=raxml_tree, charGenus=charGenus, missing_in_ott_tree=ott_trees$missing_in_ott_tree, complexity=complexity)),
 tar_target(constraints, MakeConstraints(paleodb_tree)),
 tar_target(print_constraints, PrintConstraints(constraints))
)
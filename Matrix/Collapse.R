# Create collapsed matrix, merging according to the CollapsibleSpicules.csv sheet

clean_punctuation <- function(x) {
  return(gsub("\\(", '.', gsub(" ", ".", gsub("\\)", ".", gsub("\\,", ".", x)))))
}

collapsible <- read.csv("CollapsibleSpicules.csv")
raw_char <- read.csv("CharacterMatrix.csv")
colnames(raw_char) <- gsub("X2", "2", colnames(raw_char))
collapsible_to_change <- subset(collapsible, Collapse!="No")
for (i in sequence(nrow(collapsible_to_change))) {
  moving_to <-  clean_punctuation(gsub("To ", "", collapsible_to_change$Collapse[i]))
  moving_from <- clean_punctuation(collapsible_to_change$Spicule[i])
  for (row_index in sequence(nrow(raw_char))) {
    if(raw_char[row_index, moving_from]==1) {
      raw_char[row_index, moving_to]<-1
    }
  }
}
unneeded_cols <- colnames(raw_char) %in% clean_punctuation(collapsible_to_change$Spicule)
final_char <- raw_char[,!unneeded_cols]
write.csv(final_char, file="FinalTraits.csv")
final_char_no_pairing <- final_char[, !(colnames(final_char) %in% "Pairing")]
complexity <- apply(final_char_no_pairing[,-1], 1, sum, na.rm=TRUE)
names(complexity) <- final_char_no_pairing$Species
complexity_df <- data.frame(Species=names(complexity), SpiculeTypes=complexity)
write.csv(complexity_df, file="Complexity.csv", row.names=FALSE)

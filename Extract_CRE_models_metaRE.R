library(metaRE)
library(tidyverse)
library(tidyr)
library(dplyr)
#load upstream regions----

load("~/upstream_seq_A.thaliana.Rdata")
promoters_1500 <- setNames(as.character(seq1500$`Flank (Gene)`), seq1500$`Gene stable ID`)
#hexamer
elements <- enumerateOligomers(promoters_1500, k = 6)

#Repeats
TGTC_elements <- c()
for(i in 1:50){
  TGTC_elements[i] <- paste("TGTCNN", paste(rep("N", i), collapse = ""), paste("TGTC", sep = ""), sep = "")
}

TGTC_elements <- c(paste("TGTCNN", "TGTC", sep = ""), TGTC_elements)

TGTC_GACA_elements <- c()
for(i in 1:50){
  TGTC_GACA_elements[i] <- paste("TGTCNN", paste(rep("N", i), collapse = ""), paste("NNGACA", sep = ""), sep = "")
}
TGTC_GACA_elements <- c(paste("TGTCNN", "NNGACA", sep = ""), TGTC_GACA_elements)

GACA_TGTC_elements <- c()
for(i in 1:50){
  GACA_TGTC_elements[i] <- paste("GACA", paste(rep("N", i), collapse = ""), paste("TGTC", sep = ""), sep = "")
}

GACA_TGTC_elements <- c(paste("GACA", "TGTC", sep = ""), GACA_TGTC_elements)

TGTC_TGTC_repeats <- enumeratePatterns(promoters_1500, TGTC_elements)
TGTC_GACA_repeats <- enumeratePatterns(promoters_1500, TGTC_GACA_elements)
GACA_TGTC_repeats <- enumeratePatterns(promoters_1500, GACA_TGTC_elements)

TGTC_repeats <- c(TGTC_TGTC_repeats, TGTC_GACA_repeats, GACA_TGTC_repeats)
TGTC_repeats <- GeneClassificationSparse(TGTC_repeats, geneNames = geneNames(TGTC_TGTC_repeats))

#biparite elements
elements <- enumerateDyadsWithCore(promoters_1500, k = 6, "TGTCNN", 0, 4)

#build matrix of cis-elements----
cis_genes <- lapply(1:length(elements), function(x) 
  geneNames(elements)[elements[[x]]])

motif_df <- data.frame(matrix(ncol = length(elements), nrow = length(geneNames(elements))))
colnames(motif_df) <- names(elements)

row.names(motif_df) <- geneNames(elements)

for(i in 1:length(elements)){
  motif_df[ ,i] <- row.names(motif_df) %in% cis_genes[[i]]
}


#change order of the dataset, so it can be used for scPortrait-----
motif_df <- data.frame(Gene = geneNames(elements),
                   motif_df)


df_long <- motif_df %>%
  gather(key = "Element", value = "Value", -1)

#Filter the rows based on TRUE
df_filtered <- df_long %>%
  filter(Value == TRUE)

#Select only the Gene and Element columns
motif_df <- df_filtered %>%
  select(Gene, Element)

#Arrange in ascending order
motif_df <- motif_df %>%
  arrange(Gene)

#saving to RDS file-----
saveRDS(motif_df, "~/.Rdata")







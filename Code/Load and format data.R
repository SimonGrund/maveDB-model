## Load json file
library(data.table)
library(tidyverse)

### BRCA1 Ring and BRCT domain depletion scores: https://mavedb.org/score-sets/urn:mavedb:00000081-a-1
d = fread("Data/depmap_export_2025-11-26 09_50_16.184205_subsetted.csv") 

clin_data = fread("Data/clin_data.csv")#The depletion score reported here is the number of replicates where the variant was depleted relative to the corresponding control siRNA replicate.
colnames(clin_data) = c("depmap_id", "RNAi")
# Split into WT, Location, and MT
#Make three columns out of hgvs_pro — first, skip the first 2 letters, then the next three letters in first column (WT), three last letters in second column (MT) and finally, whatever was between WT and MT in a columned named Location

d <- d %>%
  dplyr::select(-cell_line_display_name, -lineage_1, -lineage_2, -lineage_3, -lineage_4, -lineage_6)

d = left_join(clin_data, d)

#Remove "Expression Public 25Q3 " from all colnames in d
#colnames(d) = gsub("Expression Public 25Q3 ", "", colnames(d))

write.table(d, "tmp/LATEST_formatted_data.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

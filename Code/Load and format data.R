## Load json file
library(data.table)
library(tidyverse)

### BRCA1 Ring and BRCT domain depletion scores: https://mavedb.org/score-sets/urn:mavedb:00000081-a-1
d = fread("Data/urn_mavedb_00000081-a-1_scores.csv") #The depletion score reported here is the number of replicates where the variant was depleted relative to the corresponding control siRNA replicate.

# Split into WT, Location, and MT
#Make three columns out of hgvs_pro — first, skip the first 2 letters, then the next three letters in first column (WT), three last letters in second column (MT) and finally, whatever was between WT and MT in a columned named Location

d <- d %>%
  mutate(
    hgvs_trim = str_sub(hgvs_pro, 3),  # Remove "p." prefix
    WT = str_sub(hgvs_trim, 1, 3),
    MT = str_sub(hgvs_trim, -3),
    Location = str_extract(hgvs_trim, "(?<=^[A-Za-z]{3})[0-9]+(?=[A-Za-z]{3}$)")
  ) %>%
  select(hgvs_pro, WT, Location, MT, Score = score, Replicates = replicates) %>%
  na.omit() %>%
  mutate(
    Freq = ifelse(Replicates == 0, 0, Score / Replicates), # Calculate frequency
    Freq = ifelse(Freq > 1, 1, Freq) # Cap frequency at 1
  ) 

d$Location = as.numeric(d$Location) # Convert Location to numeric)

write.table(d, "Data/formatted_BRCA1_depletion_scores.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

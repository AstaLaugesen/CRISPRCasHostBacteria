#reset
rm(list=ls())

#read data
setwd("~/data")
clusterCRISPRs <- read.table("/tmp/allClustersCRISPRs.txt", sep="\t")

#renaming
clusterCRISPRs <- clusterCRISPRs %>%
  #extracting rownames:
  rownames_to_column() %>% 
  #substituting 
  select(rowname, matches("w_CRISPR")) 

#removing _w_CRISPR from col name
names(clusterCRISPRs) <- str_replace(names(clusterCRISPRs),"_w_CRISPR","") 

#save data for other analysis
save(clusterCRISPRs, file = "crispr/CRISPRsInClusters.RData", compress = T)

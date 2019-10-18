##Goal1: to generate unique gene lists from the comparison of 2xLPS S100A9 WT vs Mutant
##Goal2: input unique gene lists into NA for unsupervised networks and generate pathways
##Goal3: input the same unique gene list into SIGORA for enriched pathways dependent on the gene list
##Goal4: compare and identify common pathways and query the biological relevance for each condition

#setwd
#setwd("C:/Users/beverlie/Google Drive/Proteomics Analysis on Human Monocytes/Unique gene and pathway lists/Proteomics/PropsedOmics_Proteomics")
#setwd for mac
setwd("~/Google Drive/RNA-Seq Analysis on iPSDM_MGST1 & S100A9/S100A9/DE_genes/S100A9_WTvMUT")


#load tidyverse library the includes dplyr
install.packages ("tidyverse")
library(tidyverse)
library(dplyr)
install.packages("VennDiagram")
library(VennDiagram)
install.packages("sigora")
library(sigora)

#read proteomic dataframes with paired ttest values included
WT <- read_csv("group_Wt_LPSLPS_vs_Vehicle_20190626.csv")
MT <- read_csv("group_Mt_LPSLPS_vs_Vehicle_20190626.csv")

#remove na values
filter.MT <- na.omit(MT)
filter.WT <- na.omit(WT)

#select columns: gene, Log2FC, padj, hgnc_symbol
filter.MT <- filter.MT %>% 
  select(c(1, 3, 7, 10)) 
filter.WT <- filter.WT %>% 
  select(c(1, 3, 7, 10)) 

#cut off log2FoldChange value >1.5 & < -1.5
filtered.MT <- filter(filter.MT, log2FoldChange > 1.5 | log2FoldChange < -1.5)
filtered.WT <- filter(filter.WT, log2FoldChange > 1.5 | log2FoldChange < -1.5)


#find unique genes for each comparison, no FC cutoff 
unique.MT <- anti_join(filter.MT, filter.WT, by = "gene" )
unique.WT <- anti_join(filter.WT, filter.MT, by = "gene")
#find gene count of shared genes, no FC cutoff 
shared.genes <- merge(filter.MT, filter.WT, by= "gene", all = FALSE)


#find unique genes for each comparison, WITH FC cutoff 
FCunique.MT <- anti_join(filtered.MT, filtered.WT, by = "gene" )
FCunique.WT <- anti_join(filtered.WT, filtered.MT, by = "gene")
#find gene count of shared genes, WITH FC cutoff 
FCshared.genes <- merge(filtered.MT, filtered.WT, by= "gene", all = FALSE)



#save csv file and look up pathways on NetworkAnalyst
write.csv(unique.WT, file= "UniqueWildtype_LPSLPSvVehicle_genes.csv")
##on NA 1st order, 
##on NA zero order, 

write.csv(unique.MT, file = "UniqueMutant_LPSLPSvVehicle_genes.csv")
##on NA 1st order, 
##on NA zero order, 
##on minimum order, 


#save FC csv file and look up pathways on NetworkAnalyst
write.csv(FCunique.WT, file= "FCUniqueWildtype_LPSLPSvVehicle_genes.csv")
##on NA 1st order, 170 nodes, 172 edges, 21 seeds
##on NA zero order, 

write.csv(FCunique.MT, file = "FCUniqueMutant_LPSLPSvVehicle_genes.csv")
##on NA 1st order, 1089 nodes, 1455 edges, 96 seeds
##on NA zero order, 
##on minimum order, 


#check number of genes
nrow(unique.MT)
nrow(unique.WT)

#add unique DE protein lists and category names, no FC cutoff 
vennplot <- draw.pairwise.venn(area1 = 780,area2 = 993,cross.area = 642,category = c("Wildtype_LPSLPSvVehicle", "Mutant_LPSLPSvVehicle"))
#add unique DE protein lists and category names, WITH FC cutoff 
FCvennplot <- draw.pairwise.venn(area1 = 236,area2 = 329,cross.area = 201,category = c("Wildtype_LPSLPSvVehicle", "Mutant_LPSLPSvVehicle"))

#save svg and manipulate on inkscape
#ggsave(vennplot, file = "proteomics_vennplot.svg", device = "svg")

########check if there are any shared pathways btwn WT and Mutant, no FC cutoff
FirstWT.pathway <- read_csv("1stNA_LPSLPSvVeh_Wildtype_Pathways.csv")
FirstMT.pathway <- read_csv("1stNA_LPSLPSvVeh_Mutant_Pathways.csv")
First.shared.pathways <- merge(FirstMT.pathway, FirstWT.pathway, by= "Pathway", all = FALSE)
write.csv(First.shared.pathways, file="S100A9_LPSLPSvVehicle_1stSharedPathways.csv")

#ZeroWT.pathway <- read_csv("ZeroNA_LPSLPSvVeh_Wildtype_Pathways.csv")
#ZeroMT.pathway <- read_csv("ZeroNA_LPSLPSvVeh_Mutant_Pathways.csv")
#Zero.sared.pathways <- merge(ZeroMT.pathway, ZeroWT.pathway, by= "Pathway", all = FALSE)
#write.csv(Zero.shared.pathways, file="MGST1_LPSLPSvVehicle_ZeroSharedPathways.csv")

########check if there are any shared pathways btwn WT and Mutant WITH FC cutoff
FC.FirstWT.pathway <- read_csv("1stNA_LPSLPSvVeh_Wildtype_Pathways.csv")
FC.FirstMT.pathway <- read_csv("1stNA_LPSLPSvVeh_Mutant_Pathways.csv")
FC.First.shared.pathways <- merge(FC.FirstMT.pathway, FC.FirstWT.pathway, by= "Pathway", all = FALSE)
write.csv(FC.First.shared.pathways, file="S100A9_LPSLPSvVehicle_FC1stSharedPathways.csv")

#########################################script from Travis to generate pathway list from SIGORA
#sigora with no FC cutoff
sigora_unique.WT <- sigora(reaH, level = 4, queryList = unique.WT$gene) %>%
  pluck("summary_results") %>%
  filter(Bonferroni <= 0.001)
#found 0 enriched LPSLPSvVehicle pathways from Sigora

sigora_unique.MT <- sigora(reaH, level = 4, queryList = unique.MT$gene) %>%
  pluck("summary_results") %>% 
  filter(Bonferroni <= 0.001)
#found 3 enriched pathways

#save SIGORA lists as .csv
write.csv(sigora_unique.WT, file = "Sigora_LPSLPSvVehicle_WT.csv")
write.csv(sigora_unique.MT, file = "Sigora_LPSLPSvVehicle_MT.csv")


#sigora with WITH FC cutoff
sigora_FCunique.WT <- sigora(reaH, level = 4, queryList = FCunique.WT$gene) %>%
  pluck("summary_results") %>%
  filter(Bonferroni <= 0.001)
#found 0 enriched LPSLPSvVehicle pathways from Sigora

sigora_FCunique.MT <- sigora(reaH, level = 4, queryList = FCunique.MT$gene) %>%
  pluck("summary_results") %>% 
  filter(Bonferroni <= 0.001)
#found 3 enriched pathways


#save SIGORA lists as .csv
write.csv(sigora_unique.WT, file = "Sigora_LPSLPSvVehicle_WT.csv")
write.csv(sigora_unique.MT, file = "Sigora_LPSLPSvVehicle_MT.csv")
##################find common pathways btwn NA and Sigora
WT.common.pathways <- merge(FirstWT.pathway, sigora_unique.WT, by.x = "Pathway", by.y = "description", all = FALSE)

MT.common.pathways <- merge(FirstMT.pathway, sigora_unique.MT, by.x = "Pathway", by.y = "description", all = FALSE)

############################################find shared pathway lists from NA and Sigora
NA.LPS <- read_csv("Pairedttest_LPSvVeh_Proteomics.csv")






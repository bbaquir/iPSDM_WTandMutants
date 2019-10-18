##Goal1: to generate unique gene lists from the comparison of 2xLPS MGST1 WT vs Mutant
##Goal2: input unique gene lists into NA for unsupervised networks and generate pathways
##Goal3: input the same unique gene list into SIGORA for enriched pathways dependent on the gene list
##Goal4: compare and identify common pathways and query the biological relevance for each condition

#setwd
setwd("C:/Users/beverlie/Google Drive/RNA-Seq Analysis on iPSDM_MGST1 & S100A9/MGST1/DE_genes/MGST1_WTvMUT")
#setwd for mac
setwd("~/Google Drive/RNA-Seq Analysis on iPSDM_MGST1 & S100A9/MGST1/DE_genes/MGST1_WTvMUT")


#load tidyverse library the includes dplyr
#install.packages ("tidyverse")
library(tidyverse)
library(dplyr)
#install.packages("VennDiagram")
library(VennDiagram)
#install.packages("sigora")
library(sigora)

#read proteomic dataframes with paired ttest values included
MT <- read_csv("group_Mt_LPSLPS_vs_Vehicle_20190626.csv")
WT <- read_csv("group_Wt_LPSLPS_vs_Vehicle_20190626.csv")

#remove na values
filter.MT <- na.omit(MT)
filter.WT <- na.omit(WT)

#select columns: gene, Log2FC, padj, hgnc_symbol
filter.MT <- filter.MT %>% 
  select(c(1, 3, 7, 10)) 
filter.WT <- filter.WT %>% 
  select(c(1, 3, 7, 10)) 


############################DID NOT FILTER ITS NOT WORKING!
###filter.LPS$p.adj <- p.adjust(filter.LPS$d_pairedttest, method = "BH")


###filter.LPSLPS <- filter.LPSLPS %>% select(c(2, 12:13))
####filter.LPSLPS$p.adj <- p.adjust(filter.LPSLPS$dd_pairedttest, method = "BH")

###padj values are already <0.05


#cut off log2FoldChange value >1.5 & < -1.5
filtered.MT <- filter(filter.MT, log2FoldChange > 1.5 | log2FoldChange < -1.5)
filtered.WT <- filter(filter.WT, log2FoldChange > 1.5 | log2FoldChange < -1.5)

#df <- filter.MT %>% filter (log2FoldChange > 1.5 & log2FoldChange < -1.5)


#filtered.MT <- filter.MT %>% 
  filter(log2FoldChange < -1.5) 

#filter(log2FoldChange > 1.5)
#filtered.WT <- filter.WT %>%
  filter(padj < 0.5)

#find unique genes for each comparison, no FC cutoff 
unique.MT <- anti_join(filter.MT, filter.WT, by = "gene" )
unique.WT <- anti_join(filter.WT, filter.MT, by = "gene")

#find gene count of shared genes, no FC cutoff 
shared.genes <- merge(filter.MT, filter.WT, by= "gene", all = FALSE)

#save csv file and look up pathways on NetworkAnalyst
write.csv(unique.WT, file= "UniqueWildtype_LPSLPSvVehicle_genes.csv")
##on NA 1st order, nodes=4217, edges=8083, seeds=423
##on NA zero order, nodes=102, edges=127, seeds=102

write.csv(unique.MT, file = "UniqueMutant_LPSLPSvVehicle_genes.csv")
##on NA 1st order, nodes=5548, edges=12508, seeds=657
##on NA zero order, nodes=199, edges= 261, seeds=199
##on minimum order, nodes= 1417, edges=5335 , seeds=657 



#check number of genes
nrow(unique.MT)
nrow(unique.WT)

#add unique DE protein lists and category names, no FC cutoff 
vennplot <- draw.pairwise.venn(area1 = 2077,area2 = 2320,cross.area = 1568,category = c("Wildtype_LPSLPSvVehicle", "Mutant_LPSLPSvVehicle"))

#save svg and manipulate on inkscape
#ggsave(vennplot, file = "proteomics_vennplot.svg", device = "svg")

########check if there are any shared pathways btwn WT and Mutant
FirstWT.pathway <- read_csv("1stNA_LPSLPSvVeh_Wildtype_Pathways.csv")
FirstMT.pathway <- read_csv("1stNA_LPSLPSvVeh_Mutant_Pathways.csv")
First.shared.pathways <- merge(FirstMT.pathway, FirstWT.pathway, by= "Pathway", all = FALSE)
write.csv(First.shared.pathways, file="MGST1_LPSLPSvVehicle_1stSharedPathways.csv")

ZeroWT.pathway <- read_csv("ZeroNA_LPSLPSvVeh_Wildtype_Pathways.csv")
ZeroMT.pathway <- read_csv("ZeroNA_LPSLPSvVeh_Mutant_Pathways.csv")
Zero.shared.pathways <- merge(ZeroMT.pathway, ZeroWT.pathway, by= "Pathway", all = FALSE)
write.csv(Zero.shared.pathways, file="MGST1_LPSLPSvVehicle_ZeroSharedPathways.csv")


#########################################script from Travis to generate pathway list from SIGORA
#sigora with no FC cutoff
sigora_unique.WT <- sigora(reaH, level = 4, queryList = unique.WT$gene) %>%
  pluck("summary_results") %>%
  filter(Bonferroni <= 0.001)
#found 6 enriched LPSLPSvVehicle pathways from Sigora
  
sigora_unique.MT <- sigora(reaH, level = 4, queryList = unique.MT$gene) %>%
  pluck("summary_results") %>% 
  filter(Bonferroni <= 0.001)
#found 10 enriched pathways


#save SIGORA lists as .csv
write.csv(sigora_unique.WT, file = "Sigora_LPSLPSvVehicle_WT.csv")
write.csv(sigora_unique.MT, file = "Sigora_LPSLPSvVehicle_MT.csv")

##################find common pathways btwn NA and Sigora
WT.common.pathways <- merge(FirstWT.pathway, sigora_unique.WT, by.x = "Pathway", by.y = "description", all = FALSE)

MT.common.pathways <- merge(FirstMT.pathway, sigora_unique.MT, by.x = "Pathway", by.y = "description", all = FALSE)

############################################find shared pathway lists from NA and Sigora
NA.LPS <- read_csv("Pairedttest_LPSvVeh_Proteomics.csv")






# GO Enrichment Test
library(topGO)
library(dplyr)

#import all go terms as a dataframe 
ALL_go_terms <- read.delim("../transcriptome/GO_terms/HF4_all_go_terms_topgo.txt")

#import high fst contigs and annotation
high_fst_contigs <- read.table("../transcriptome/high_fst_contigs")
high_fst_fsts <- read.table("out_oct6.weir.fst", header = T)
HF4_annotation <- read.csv("../transcriptome/HF4_final_annotation.csv")

#filter to select "subjects" that are high fst
high_fst_subjects <- HF4_annotation %>% 
  filter(Query %in% high_fst_contigs$V1)

#add fst values to subjects for high fst snps

highest_fst <- filter(high_fst_fsts, WEIR_AND_COCKERHAM_FST > 0.8)

high_fst_with_fst <- HF4_annotation %>% 
 filter(Query %in% highest_fst$CHROM)
  
write.csv(x = high_fst_with_fst, file = "high_fst_annotation.csv")

#filter to pull GO terms for the high fst contigs
high_fst_GOs <- ALL_go_terms %>% 
  filter(Entry %in% high_fst_subjects$Subject)

# read in GO annotations with TopGO 
geneID2GO <- readMappings(file = "../transcriptome/GO_terms/HF4_all_go_terms_topgo.txt")

# define "gene universe" or all the GO terms in your dataset 
all_genes <- names(geneID2GO)

#read in GO annotations for the genes of interest (high fst) with TopGO
genesOfInterest <- as.character(high_fst_GOs$Entry)
# define "genes of interest" as a subset of the total dataset 
geneList <- factor(as.integer(all_genes %in% genesOfInterest))
names(geneList) <- all_genes

# Create TopGOData object for Biological Processes 
GO_data_BP <- new("topGOdata", description = "Pleaxuara homomalla SNPs",
               ontology = "BP", allGenes= geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize =10)
GO_data_BP

#Review significantly enriched genes  
sig_genes_BP <- sigGenes(GO_data_BP)
str(sig_genes_BP)
#Number of significantly enriched genes for Biological Processes 
numSigGenes(GO_data_BP)
#97

#Run Fisher's Exact Test for Biological Processes 
#Note that this does not correct for multiple comparisons 
BP_result_Fisher <- runTest(GO_data_BP, algorithm = "classic", statistic = "fisher", )
BP_result_Fisher
#9 terms have raw p<0.01 

#BP_list <- usedGO(object = GO_data_BP)

BP_fisher_only <- GenTable(GO_data_BP, classicFisher = BP_result_Fisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(BP_list))
write.table(BP_fisher_only, file = "BP_fisher_only_nodesize_10.txt") 

p.adjust(BP_fisher_only$classicFisher, method = "fdr")
#no terms are significantly enriched following FDR adjustment

#Run Kolmogorov-Smirnov Test
#BP_result_KS <- runTest(GO_data_BP, algorithm = "classic", statistic = "ks")
#BP_result_KS

# Analysis of BP results
all_BP_results <- GenTable(GO_data_BP, classicFisher = BP_result_Fisher, 
                           orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 500)
all_BP_results

#BP_ks_pValue.classic <- score(BP_result_KS) 
#BP_fisher_pValue.classic <- score(BP_result_Fisher)
#BP_fisher <- as.data.frame(BP_fisher_pValue.classic)
#write.csv(BP_fisher, "../transcriptome/GO_analyses/BP_fisher_10nodes.csv", quote = F)
#BP_fisher <- read.csv("../transcriptome/GO_analyses/BP_fisher_10nodes.csv", header = T)
#BP_KS <- data.frame(BP_ks_pValue.classic)
#write.csv(BP_KS, "../transcriptome/GO_analyses/BP_KS.csv", quote = F)
#BP_KS <- read.csv("../transcriptome/GO_analyses/BP_KS.csv", header = T)
#BP_all <- merge.data.frame(BP_KS, BP_fisher, by = "GO_term")
#BP_sig_pvals <- BP_fisher %>% 
#  filter(BP_fisher_pValue.classic < 0.01)
#write.csv(BP_sig_pvals, "../transcriptome/GO_analyses/BP_signficant_GOs_10nodes.csv", quote = F)

showSigOfNodes(GO_data_BP, score(BP_result_Fisher), firstSigNodes = 10, useInfo = 'all')
printGraph(GO_data_BP, BP_result_Fisher, firstSigNodes = 10, fn.prefix = "tGO", useInfo = 'all')

# Analyze each GO term to see what genes are attached
BP_terms_df <- all_BP_results %>% 
  filter(classicFisher < 0.01)

BP_terms <- BP_terms_df$GO.ID
BPgenes <- genesInTerm(GO_data_BP, BP_terms)

GO_0000027 <- as.vector(BPgenes$`GO:0000027`)
as.data.frame(GO_0000027)
GO_0000027_subjects <- high_fst_subjects %>% 
  filter(Subject %in% GO_0000027) %>% 
  select(Query, Subject) %>% 
  mutate(GO_ID = c("GO:0000027"))
write.csv(GO_0000027_subjects, "GO_0000027.csv")

GO_0006511 <- as.vector(BPgenes$`GO:0006511`)
as.data.frame(GO_0006511)
GO_0006511_subjects <- high_fst_subjects %>% 
  filter(Subject %in% GO_0006511) %>% 
  select(Query, Subject) %>% 
  mutate(GO_ID = c("GO:0006511"))
write.csv(GO_0006511_subjects, "GO_0006511.csv")

GO_0006739 <- as.vector(BPgenes$`GO:0006739`)
as.data.frame(GO_0006739)
GO_0006739_subjects <- high_fst_subjects %>% 
  filter(Subject %in% GO_0006739) %>% 
  select(Query, Subject) %>% 
  mutate(GO_ID = c("GO:0006739"))
write.csv(GO_0006739_subjects, "GO_0006739.csv")

GO_0019941 <- as.vector(BPgenes$`GO:0019941`)
as.data.frame(GO_0019941)
GO_0019941_subjects <- high_fst_subjects %>% 
  filter(Subject %in% GO_0019941) %>% 
  select(Query, Subject) %>% 
  mutate(GO_ID = c("GO:0019941"))
write.csv(GO_0006511_subjects, "GO_0019941.csv")

GO_0043632 <- as.vector(BPgenes$`GO:0043632`)
as.data.frame(GO_0043632)
GO_0043632_subjects <- high_fst_subjects %>% 
  filter(Subject %in% GO_0043632) %>% 
  select(Query, Subject) %>% 
  mutate(GO_ID = c("GO:0043632"))
write.csv(GO_0043632_subjects, "GO_0043632.csv")

GO_0044262 <- as.vector(BPgenes$`GO:0044262`)
as.data.frame(GO_0044262) 
GO_0044262_subjects <- high_fst_subjects %>% 
  filter(Subject %in% GO_0044262) %>% 
  select(Query, Subject) %>% 
  mutate(GO_ID = c("GO:0044262"))
write.csv(GO_0044262_subjects, "GO_0044262.csv")

GO_0044265 <- as.vector(BPgenes$`GO:0044265`)
as.data.frame(GO_0044265)
GO_0044265_subjects <- high_fst_subjects %>% 
  filter(Subject %in% GO_0044265) %>% 
  select(Query, Subject) %>% 
  mutate(GO_ID = c("GO:0044265"))
write.csv(GO_0044265_subjects, "GO_0044265.csv")

GO_0044271 <- as.vector(BPgenes$`GO:0044271`)
as.data.frame <- GO_0044271 
GO_0044271_subjects <- high_fst_subjects %>% 
  filter(Subject %in% GO_0044271) %>% 
  select(Query, Subject) %>% 
  mutate(GO_ID = c("GO:0044271"))
write.csv(GO_0044271_subjects, "GO_0044271.csv")

GO_1901575 <- as.vector(BPgenes$`GO:1901575`)
as.data.frame <- GO_1901575 
GO_1901575_subjects <- high_fst_subjects %>% 
  filter(Subject %in% GO_1901575) %>% 
  select(Query, Subject) %>% 
  mutate(GO_ID = c("GO:1901575"))
write.csv(GO_1901575_subjects, "GO_1901575.csv")

#########################################################################################
# Cellular Component 
# Create TopGOData object for Cellular Component 
GO_data_CC <- new("topGOdata", description = "Pleaxuara homomalla SNPs",
                  ontology = "CC", allGenes= geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
GO_data_CC

#Review significantly enriched genes  
sig_genes_CC <- sigGenes(GO_data_CC)
str(sig_genes_CC)
#Number of significantly enriched genes for Biological Processes 
numSigGenes(GO_data_CC)
#114

#Run Fisher's Exact Test for Biological Processes 
#Note that this does not correct for multiple comparisons 
CC_result_Fisher <- runTest(GO_data_CC, algorithm = "classic", statistic = "fisher")
CC_result_Fisher
# 5 terms with p < 0.01

CC_all <- usedGO(GO_data_CC)
CC_fisher_only <- GenTable(GO_data_CC, classicFisher = CC_result_Fisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(CC_all))
write.table(CC_fisher_only, file = "CC_fisher_only_10nodes.txt") 
p.adjust(CC_fisher_only$classicFisher, method = "fdr")
# no terms are signifcant following FDR correction

#Run Kolmogorov-Smirnov Test
#CC_result_KS <- runTest(GO_data_CC, algorithm = "classic", statistic = "ks")
#CC_result_KS

# Analysis of CC results
all_CC_results <- GenTable(GO_data_CC, classicFisher = CC_result_Fisher,
                           orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 150)
all_CC_results

#CC_ks_pValue.classic <- score(CC_result_KS) 
#CC_fisher_pValue.classic <- score(CC_result_Fisher)
#CC_fisher <- as.data.frame(CC_fisher_pValue.classic)
#write.csv(CC_fisher, "../transcriptome/GO_analyses/CC_fisher.csv", quote = F)
#CC_fisher <- read.csv("../transcriptome/GO_analyses/CC_fisher.csv", header = T)
#CC_KS <- data.frame(CC_ks_pValue.classic)
#write.csv(CC_KS, "../transcriptome/GO_analyses/CC_KS.csv", quote = F)
#CC_KS <- read.csv("../transcriptome/GO_analyses/CC_KS.csv", header = T)
#CC_all <- merge.data.frame(CC_KS, CC_fisher, by = "GO_terms")
#CC_sig_pvals <- CC_all %>% 
#  filter(KS_pvals < 0.05 & Fisher_pvals < 0.05)
#write.csv(CC_sig_pvals, "../transcriptome/GO_analyses/CC_signficant_GOs.csv", quote = F)

showSigOfNodes(GO_data_CC, score(CC_result_Fisher), firstSigNodes = 10, useInfo = 'all')
printGraph(GO_data_CC, CC_result_Fisher, firstSigNodes = 10, fn.prefix = "tGO", useInfo = 'all')

# Analyze each GO term to see what genes are attached
CC_terms_df <- all_CC_results %>% 
  filter(classicFisher < 0.01)

CC_terms <- CC_terms_df$GO.ID
CCgenes <- genesInTerm(GO_data_CC, CC_terms)

GO_0005622 <- as.vector(CCgenes$`GO:0005622`)
as.data.frame(GO_0005622)
GO_0005622_subjects <- high_fst_subjects %>% 
  filter(Subject %in% GO_0005622) %>% 
  select(Query, Subject) %>% 
  mutate(GO_ID = c("GO:0005622"))
write.csv(GO_0005622_subjects, "GO_0005622.csv")

GO_0043227 <- as.vector(CCgenes$`GO:0043227`)
as.data.frame(GO_0043227)
GO_0043227_subjects <- high_fst_subjects %>% 
  filter(Subject %in% GO_0043227) %>% 
  select(Query, Subject) %>% 
  mutate(GO_ID = c("GO:0043227"))
write.csv(GO_0043227_subjects, "GO_0043227.csv")

GO_0043226 <- as.vector(CCgenes$`GO:0043226`)
as.data.frame(GO_0043226)
GO_0043226_subjects <- high_fst_subjects %>% 
  filter(Subject %in% GO_0043226) %>% 
  select(Query, Subject) %>% 
  mutate(GO_ID = c("GO:0043226"))
write.csv(GO_0043226_subjects, "GO_0043226.csv")

GO_0043229 <- as.vector(CCgenes$`GO:0043229`)
as.data.frame(GO_0043229)
GO_0043229_subjects <- high_fst_subjects %>% 
  filter(Subject %in% GO_0043229) %>% 
  select(Query, Subject) %>% 
  mutate(GO_ID = c("GO:0043229"))
write.csv(GO_0043229_subjects, "GO_0043229.csv")

GO_0044424 <- as.vector(CCgenes$`GO:0044424`)
as.data.frame(GO_0044424)
GO_0044424_subjects <- high_fst_subjects %>% 
  filter(Subject %in% GO_0044424) %>% 
  select(Query, Subject) %>% 
  mutate(GO_ID = c("GO:0044424"))
write.csv(GO_0044424_subjects, "GO_0044424.csv")

###################################################################################################
# Molecular Function 
# Create TopGOData object for Molecular Function
GO_data_MF <- new("topGOdata", description = "Pleaxuara homomalla SNPs",
                  ontology = "MF", allGenes= geneList, annot = annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
GO_data_MF

#Review significantly enriched genes  
sig_genes_MF <- sigGenes(GO_data_MF)
str(sig_genes_MF)
#Number of significantly enriched genes for Molecular Functions 
numSigGenes(GO_data_MF)
#121

#Run Fisher's Exact Test for Molecular Functions 
#Note that this does not correct for multiple comparisons 
MF_result_Fisher <- runTest(GO_data_MF, algorithm = "classic", statistic = "fisher")
MF_result_Fisher
#1 term with p < 0.01

MF_all <- usedGO(GO_data_MF)

MF_fisher_only <- GenTable(GO_data_MF, classicFisher = MF_result_Fisher, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(MF_all))
write.table(MF_fisher_only, file = "MF_fisher_only_10_nodes.txt") 
p.adjust(MF_fisher_only$classicFisher, method = "fdr")

#Run Kolmogorov-Smirnov Test
#MF_result_KS <- runTest(GO_data_MF, algorithm = "classic", statistic = "ks")
#MF_result_KS

# Analysis of MF results
all_MF_results <- GenTable(GO_data_MF, classicFisher = MF_result_Fisher,
                           orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 100)
all_MF_results

#MF_ks_pValue.classic <- score(MF_result_KS) 
#MF_fisher_pValue.classic <- score(MF_result_Fisher)
#MF_fisher <- as.data.frame(MF_fisher_pValue.classic)
#write.csv(MF_fisher, "../transcriptome/GO_analyses/MF_fisher.csv", quote = F)
#MF_fisher <- read.csv("../transcriptome/GO_analyses/MF_fisher.csv", header = T)
#MF_KS <- data.frame(MF_ks_pValue.classic)
#write.csv(MF_KS, "../transcriptome/GO_analyses/MF_KS.csv", quote = F)
#MF_KS <- read.csv("../transcriptome/GO_analyses/MF_KS.csv", header = T)
#MF_all <- merge.data.frame(MF_KS, MF_fisher, by = "GO_terms")
#MF_sig_pvals <- MF_all %>% 
#  filter(KS_pvals < 0.05 & Fisher_pvals < 0.05)
#write.csv(MF_sig_pvals, "../transcriptome/GO_analyses/MF_signficant_GOs.csv", quote = F)


showSigOfNodes(GO_data_MF, score(MF_result_Fisher), firstSigNodes = 10, useInfo = 'all')
printGraph(GO_data_MF, MF_result_Fisher, firstSigNodes = 10, fn.prefix = "tGO", useInfo = 'all')

# Analyze each GO term to see what genes are attached
MF_terms_df <- all_MF_results %>% 
  filter(classicFisher < 0.01)

MF_terms <- MF_terms_df$GO.ID
MFgenes <- genesInTerm(GO_data_MF, MF_terms)

GO_0015103 <- as.vector(MFgenes$`GO:0015103`)
as.data.frame(GO_0015103)
GO_0015103_subjects <- high_fst_subjects %>% 
  filter(Subject %in% GO_0015103) %>% 
  select(Query, Subject) %>% 
  mutate(GO_ID = c("GO:0015103"))
write.csv(GO_0015103_subjects, "GO_0015103.csv")
##################################################################################################
#Assign genes to subject, species, and protein names 
go_genes <- read.delim("GO_genes.tab")
go_enriched <- read.csv("../transcriptome/GO_analyses/enriched_genes.csv")
go_total <- merge(y=go_genes, x=go_enriched, by = "Blast.Subject")
write.csv(go_total, "enriched_genes_names.csv")

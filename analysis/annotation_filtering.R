# transcriptome annotation cleaning
library(dplyr)

# Step 1. Filter annotations from the ensembl cnidarians 
HF4_ensembl_cnidaria <- read.table("../transcriptome/Annotations/HF4_ensembl_annotation_april29")
HF4_ensembl_cnidaria_sorted <- HF4_ensembl_cnidaria %>% 
  group_by(V1) %>% #group by the transcript
  filter(V4 > 100) %>% #only keep hits that overlap the query by at least 100 bp
  filter(row_number() == 1) #only keep the hit with the best e-value 
# 17,911 transcripts annotated with Ensebml Cnidaria  


KG1_ensebml_cnidaria <- read.table("../transcriptome/Annotations/KG1_ensembl_annotation_april29")
KG1_ensembl_cnidaria_sorted <- KG1_ensebml_cnidaria %>% 
  group_by(V1) %>% 
  filter(V4 > 100) %>% 
  filter(row_number() == 1)
# 25,645 transcripts annotated with Ensebml Cnidaria

# Step 2. Filter annotations from the ensembl metazoa data set
HF4_ensembl_metazoa <- read.table("../transcriptome/Annotations/HF4_ensembl_metazoa_annotation_may6")

HF4_metazoa_unique <- HF4_ensembl_metazoa[!HF4_ensembl_metazoa$V1 %in% HF4_ensembl_cnidaria_sorted$V1,]

HF4_ensembl_metazoa_sorted <- HF4_metazoa_unique %>% 
  group_by(V1) %>% 
  filter(V4 > 100) %>% 
  filter(row_number() == 1)
#4,137 transcripts annotated with Ensebml Metazoa

KG1_ensembl_metazoa <- read.table("../transcriptome/Annotations/KG1_ensembl_metazoa_annotation_may9") 

KG1_metazoa_unique <- KG1_ensembl_metazoa[!KG1_ensembl_metazoa$V1 %in% KG1_ensembl_cnidaria_sorted$V1,]

KG1_ensembl_metazoa_sorted <- KG1_metazoa_unique %>% 
  group_by(V1) %>% 
  filter(V4 > 100) %>% 
  filter(row_number() == 1)
# 7,433 transcripts annotated with Ensembl Metazoa

# Step 3. Filter annotations from the Uniprot cnidarian data set 

HF4_uniprot_cnidaria <- read.table("../transcriptome/Annotations/HF4_uniprot_annotation_april29")

HF4_uniprot_unique <- HF4_uniprot_cnidaria[!HF4_uniprot_cnidaria$V1 %in% HF4_ensembl_cnidaria_sorted$V1,]

HF4_uniprot_unique <- HF4_uniprot_unique[!HF4_uniprot_unique$V1 %in% HF4_ensembl_metazoa_sorted$V1,]

HF4_uniprot_cnidaria_sorted <- HF4_uniprot_unique %>% 
  group_by(V1) %>% 
  filter(V4 > 100) %>% 
  filter(row_number() == 1)
#1,487 transcripts annotated with Uniprot Cnidaria 
# We will use this later, write out to .txt
write.table(HF4_uniprot_cnidaria_sorted, "../transcriptome/HF4_uniprot_cnidaria_sorted.txt",
            sep = '\t', quote = F)


KG1_uniprot_cnidaria <- read.table("../transcriptome/Annotations/KG1_uniprot_annotation_april29")

KG1_uniprot_unique <- KG1_uniprot_cnidaria[!KG1_uniprot_cnidaria$V1 %in% KG1_ensembl_cnidaria_sorted$V1,]

KG1_uniprot_unique <- KG1_uniprot_cnidaria[!KG1_uniprot_cnidaria$V1 %in% KG1_ensembl_cnidaria_sorted$V1,]

KG1_uniprot_cnidaria_sorted <- KG1_uniprot_unique %>% 
  group_by(V1) %>% 
  filter(V4 > 100) %>% 
  filter(row_number() == 1)
#8,150 transcripts annotated with Uniprot Cnidaria 
# We will use this later, write out to .txt
write.table(KG1_uniprot_cnidaria_sorted, "../transcriptome/KG1_uniprot_cnidaria_sorted.txt",
            sep = '\t', quote = F)

# Step 4. Concatenate the three annotation files into one
HF4_ensembl_annotation <- rbind(HF4_ensembl_cnidaria_sorted, HF4_ensembl_metazoa_sorted)
HF4_final_annotation <- rbind(HF4_ensembl_annotation, HF4_uniprot_cnidaria_sorted)
write.csv(HF4_final_annotation, "../HF4_final_annotation.csv")

KG1_ensebml_annotation <- rbind(KG1_ensembl_cnidaria_sorted, KG1_ensembl_metazoa_sorted)
KG1_final_annotation <- rbind(KG1_ensebml_annotation, KG1_uniprot_cnidaria_sorted)
write.csv(KG1_final_annotation, "../transcriptome/KG1_final_annotation.csv")



# GO terms for transcripts

# Step 1. Pull GO terms for the Ensembl Cnidaria Annotation
GO_ensembl_cnidaria <- read.delim("../transcriptome/GO_terms/Ensembl_cnidaria_GO.txt", sep = '\t', header = T)

HF4_ensembl_cnidaria_GO <- GO_ensembl_cnidaria %>% 
  filter(Transcript.stable.ID %in% HF4_ensembl_cnidaria_sorted$V2)

KG1_ensembl_cnidaria_GO <- GO_ensembl_cnidaria %>% 
  filter(Transcript.stable.ID %in% KG1_ensembl_cnidaria_sorted$V2)

# Step 2. Pull GO terms for the Ensembl Metazoa Annotation 
GO_ensembl_metazoa <- read.delim("../transcriptome/GO_terms/Ensembl_metazoa_GO.txt", sep = '\t', header = T)

HF4_ensembl_metazoa_GO <- GO_ensembl_metazoa %>% 
  filter(Transcript.stable.ID %in% HF4_ensembl_metazoa_sorted$V2)

KG1_ensembl_metazoa_GO <- GO_ensembl_metazoa %>% 
  filter(Transcript.stable.ID %in% KG1_ensembl_metazoa_sorted$V2)

# Step 3. Pull Go Terms from Uniprot dataset 
Uniprot_GO_terms <- read.delim("../transcriptome/GO_terms/uniprot-cnidaria-swiss-trembl-goterms.txt", sep ='\t')

# Note that here we read in the Unirpot cnidarian annotations, but after we have replaced 
# "|" which are typically used in Uniprot datasets with tabs 

HF4_uniprot_cnidaria_annotation <- read.table("../transcriptome/HF4_uniprot_cnidaria_sorted.txt",
                                              header = T)
HF4_uniprot_cnidaria_GO <- Uniprot_GO_terms %>% 
  filter(Entry.name %in% HF4_uniprot_cnidaria_annotation$V5)

KG1_uniprot_cnidaria_annotation <- read.table("../transcriptome/KG1_uniprot_cnidaria_sorted.txt",
                                              header = T)

KG1_uniprot_cnidaria_GO <- Uniprot_GO_terms %>% 
  filter(Entry.name %in% KG1_uniprot_cnidaria_annotation$V5)

# Step 4. Concatenate GO term lists and write out 
HF4_ensembl_GO_terms <- rbind(HF4_ensembl_cnidaria_GO, HF4_ensembl_metazoa_GO)
write.csv(HF4_ensembl_GO_terms, "../transcriptome/HF4_ensembl_GO_terms.csv")

KG1_ensembl_GO_terms <- rbind(KG1_ensembl_cnidaria_GO, KG1_ensembl_metazoa_GO)
write.csv(KG1_ensembl_GO_terms, "../transcriptome/KG1_ensembl_GO_terms.csv")

write.table(HF4_uniprot_cnidaria_GO, "../transcriptome/HF4_uniprot_cnidaria_GO.txt", sep = '\t',
            quote = F)
write.table(KG1_uniprot_cnidaria_GO, "../transcriptome/KG1_uniprot_cnidaria_GO.txt", sep = '\t',
            quote = F)

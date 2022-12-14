

#############GENE ONTOLOGY ENRICHMENT ANALYSIS various


#how many genes are in this group considered a pathway
#how many of those genes are DE?
#From those two ratios, you can determine how likely it is that the pathway is affected

#Can use general enrichment approaches or statistical GSEA
#general hypergeometric enrichment: significant genes, binary, dependent on how significance is defined, compares 2+ states
#GSEA: all genes, numeric fold change values, mean fold per change, compares 2 states only


BiocManager::install("clusterProfiler")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("pathview")

library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(AnnotationHub)
library(ensembldb)
library(tidyverse)
library(ggplot2)

#dropping NA values
sigs <- na.omit(res)

#filtering only the most significantly deregulated genes
sigs <- sigs[sigs$padj <0.05 & sigs$log2FoldChange < -0.5,]

#gene IDs are genes object
genes <- rownames(sigs)

#perform enrichGO analysis
GO_res <- enrichGO(gene = genes, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

#save results as data frame
GO_results <- as.data.frame(GO_res)

#stricter p value filtering
gres0.01 <- GO_results[GO_results$p.adjust <0.02,]

#create plot of filtered results
plot <- ggplot(gres0.01, aes(Description, -log(pvalue)))+
  geom_bar(stat='identity')+
  xlab(NULL)+
  ylab("-log(pvalues)")+
  ggtitle("Downregulated Global Pathways\n GOenrich (FDR < 0.02)")

#rotate plot to visualize descriptions more easily 
plot + coord_flip()


# general enrichment analysis
#https://www.youtube.com/watch?v=Bzu4_yDcBLY

BiocManager::install("msigdbr")
BiocManager::install("fgsea")
library(msigdbr)
library(fgsea)

#viewing human database of pathways, concordant expression, not just related
H <- msigdbr(species="Homo sapiens", category = "H")
class(H)


################OVERREPRESENTATION enrichment GO analysis
#HALLMARK gene set comparisons


#converting DESeq2 results to dataframe
resdf <- as.data.frame(res.sp)

#filtering by log2foldchange
sig.sp <- res.p.df %>% dplyr::filter(log2FoldChange < -0.5)

#creating ENSEMBL IDs from rownames of signif DF
sig.ens <- unique(rownames(signif))
H.ens <- dplyr::select(H, gs_name, ensembl_gene)

#enrichment
enrich.h <- enricher(gene = sig.ens, TERM2GENE = H.ens)

#extract results
class(enrich.h) #enrich result, but is S4 datatype

#use @ symbol instead of $ to extract data from object
head(enrich.h@result)
#outputs 2 ratios, # of sig genes/# of all hallmark genes
#BgRatio = # of genes of geneset/genes in hallmark

class(enrich.h@result$GeneRatio)
#format results

enrich.H.df <- enrich.h@result %>% 
  #separate ratios into 2 columns of data
  separate(BgRatio, into=c("size.term","size.category"), sep="/") %>% 
  separate(GeneRatio, into=c("size.overlap.term", "size.overlap.category"),
           sep="/") %>% 
  #convert to numeric
  mutate_at(vars("size.term","size.category",
                 "size.overlap.term","size.overlap.category"),
            as.numeric) %>% 
  #Calculate k/K
  mutate("k.K"=size.overlap.term/size.term)

#visualizing
enrich.H.df %>% 
  dplyr::filter(p.adjust <= 0.1) %>% 
  #removing _ and HALLMARK
  mutate(Description = gsub("HALLMARK_","", Description),
         Description = gsub("_"," ", Description)) %>% 
  
  ggplot(aes(x=reorder(Description, k.K), #Reorder gene sets by k/K values
             y=k.K)) +
  geom_col() +
  theme_classic() +
  #Some more customization to pretty it up
  #Flip x and y so long labels can be read
  coord_flip() +
  #fix labels
  labs(y="Sig. genes in set / Total genes in set \nk/K",
       x="Gene set",
       title = "Downregulated Global Pathways \nHallmark (FDR < 0.15)")



##############spinal cord only


#converting DESeq2 results to dataframe
res.p.df <- as.data.frame(res.sp)

#filtering by log2foldchange
sig.sp <- res.p.df %>% dplyr::filter(log2FoldChange < -0.5)

#creating ENSEMBL IDs from rownames of signif DF
sig.ens <- unique(rownames(sig.sp))
H.ens <- dplyr::select(H, gs_name, ensembl_gene)

#enrichment
enrich.h <- enricher(gene = sig.ens, TERM2GENE = H.ens)

#extract results
class(enrich.h) #enrich result, but is S4 datatype

#use @ symbol instead of $ to extract data from object
head(enrich.h@result)
#outputs 2 ratios, # of sig genes/# of all hallmark genes
#BgRatio = # of genes of geneset/genes in hallmark

class(enrich.h@result$GeneRatio)
#format results

enrich.H.df <- enrich.h@result %>% 
  #separate ratios into 2 columns of data
  separate(BgRatio, into=c("size.term","size.category"), sep="/") %>% 
  separate(GeneRatio, into=c("size.overlap.term", "size.overlap.category"),
           sep="/") %>% 
  #convert to numeric
  mutate_at(vars("size.term","size.category",
                 "size.overlap.term","size.overlap.category"),
            as.numeric) %>% 
  #Calculate k/K
  mutate("k.K"=size.overlap.term/size.term)

#visualizing
enrich.H.df %>% 
  dplyr::filter(p.adjust <= 0.1) %>% 
  #Beautify descriptions by removing _ and HALLMARK
  mutate(Description = gsub("HALLMARK_","", Description),
         Description = gsub("_"," ", Description)) %>% 
  
  ggplot(aes(x=reorder(Description, k.K, col= Count), #Reorder gene sets by k/K values
             y=k.K)) +
  geom_col() +
  theme_classic() +
  #Some more customization to pretty it up
  #Flip x and y so long labels can be read
  coord_flip() +
  #fix labels
  labs(y="Sig. genes in set / Total genes in set \nk/K",
       x="Gene set",
       title = "Downregulated Spinal Cord Pathways \nHallmark (FDR < 0.15)")



#######################


#Gene set enrichment analysis, GSEA
#takes in fold change
#taking ES(S) value and normalizing to gene set size, bc its easier to find sig in small gene sets

#want to run GSEA in tissue specific context, group samples that make sense

####Format GSEA database ####

#format into proper structure
H.ens.ls <- H %>%
  dplyr::select(gs_name, ensembl_gene) %>%
  group_by(gs_name) %>%
  summarize(all.genes = unique(list(ensembl_gene))) %>%
  deframe()

#calculate fold change
#Extract expression data
FC <- as.data.frame(dat$E) %>% 
  #Move gene IDs from rownames to a column
  rownames_to_column("ensembl_gene_id") %>% 
  #Make long format with all expression values in 1 column
  pivot_longer(-ensembl_gene_id, 
               names_to = "libID", values_to = "expression") %>% 
  #Extract RSID and TB condition from libID
  #If this info was not in the libID, we could get it by joining
  # with dat$targets
  separate(libID, into = c("RSID","condition"), sep="_") %>% 
  #Make wide with media and tb expression in separate columns
  pivot_wider(names_from = condition, values_from = expression) %>% 
  #Calculate tb minus media fold change (delta for change)
  #Because expression values are log2, subtraction is the same as division
  mutate(delta = TB-MEDIA) %>% 
  #Calculate mean fold change per gene
  group_by(ensembl_gene_id) %>% 
  summarise(mean.delta = mean(delta, na.rm=TRUE)) %>% 
  #Arrange by descending fold change
  arrange(desc(mean.delta))


#HBC Training tutorial

## Extract the foldchanges
gs_fc <- res$log2FoldChange

## Name each fold change with the corresponding Entrez ID
names(gs_fc) <- res_entrez$entrez

## Sort fold changes in decreasing order
gs_fc <- sort(gs_fc, decreasing = TRUE)

head(gs_fc)

## GSEA using gene sets from KEGG pathways
gseaKEGG <- gseKEGG(geneList = gs_fc, # ordered named vector of fold changes (Entrez IDs are the associated names)
                    organism = "hsa", # supported organisms listed below
                    nPerm = 1000, # default number permutations
                    minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                    pvalueCutoff = 0.05, # padj cutoff value
                    verbose = FALSE)

## Extract the GSEA results
gseaKEGG_results <- gseaKEGG@result

## Write GSEA results to file
View(gseaKEGG_results)

write.csv(gseaKEGG_results, "results/gseaOE_kegg.csv", quote=F)

## Plot the GSEA plot for a single enriched pathway, `hsa03040`
gseaplot(gseaKEGG, geneSetID = 'hsa03040')

#importing meta data labels
#from GEO -> SRA run selector -> metadata
meta <- read.delim("SraRunTable.txt",header=TRUE,sep=',')

#finding possible classifications in the metadata GROUP column
unique(meta$Group)

#creating new table with alt_ID, source, and group
library(dplyr)
library(tidyverse)

source_key <- meta %>% select(sample_id_alt,source_name,Group)
#'.' is separator in dataset, '-' in in metadata

#substituting out '.' character for '-'
source_key$sample_id_alt <- gsub('-', '.', source_key$sample_id_alt)

#gene ontology analysis packages required
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)

#filtering (+) expressed >0.5 fold change

genes_to_test <- rownames(sigs[sigs$log2FoldChange > 0.5,])

#example command for creating GO terms
GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont ="BP")

#convert results to data frame
as.data.frame(GO_results)
#plot
fit <- plot(barplot(GO_results, showCategory =20))

ggplot(GO_results, aes(x=))

#cleaning data
#filter low expressed genes, e.g. remove genes with less than 10 reads across N samples

#Finding any possible class imbalance
#counts of Tissue libraries
library(stringr)

#populating categories into table
tissue.count <- data.frame(table(coldata$Tissue))

#plotting tissue count
ggplot(tissue.count, aes(Var1, Freq))+
  geom_bar(stat='identity')+
  xlab(NULL)+
  theme(axis.text.x = element_text(angle = 70,vjust = 0.5))+
  ggtitle("Tissue Sample Count")

#finding counts of Groups, total
count(coldata, "Group")

#into table
group.count <- data.frame(table(coldata$Group))

#plotting
ggplot(group.count, aes(Var1, Freq))+
  geom_bar(stat='identity')+
  xlab(NULL)+
  theme(axis.text.x = element_text(angle = 70,vjust = 0.5))+
  ggtitle("Tissue Sample Count")


############additional analyses
#only ALS spectrum MND samples

#make dataset isolates for easier compare b/t ALS and control groups
col.als <- coldata %>% filter(grepl("ALS Spectrum MND", Group))
col.con <- coldata %>% filter(grepl("Non-Neurological Control", Group))
col.sp <- coldata_comp %>% filter(grepl("Spinal", Tissue))
col.cor <- coldata_comp %>% filter(grepl("Cortex", Tissue))

#########SPINAL CORD only analysis

#naming spinal cord alt IDs
col.sp.id <- col.sp$sample_id_alt

#Combine the two filtered sets
coldata_comp <- rbind(col.als, col.con)

#remove groups that contain Other
coldata_comp <- coldata_comp %>% filter(!grepl("Other", Group))

#only 2 groups now
unique(coldata_comp$Group)

#selecting specific columns that match spinal cord libraries. Rest of IDs are saved in desktop folder
counts.sp <- counts[,c("CGND.HRA.00132", "MANY")]


#making object to store tabulated frequency of each category in group
sp.group.count <- data.frame(table(col.sp$Group))

#plot frequency
ggplot(sp.group.count, aes(Var1, Freq))+
  geom_bar(stat='identity')+
  ggtitle("ALS v Ctrl Group Size")+
  ylab("Sample count")+
  xlab("Group")


#################enrichGO package pathway enrichment
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)

#remove NA values
sig.sp <- na.omit(res.sp)

#filter results for the most significant
sig.sp <- sig.sp[sig.sp$padj <0.05 & sig.sp$log2FoldChange < -1,]

#store gene names as a new object
gene.sp <- rownames(sig.sp)

#GO analysis, supply library, gene ID tag type, and type of GO terms (BP, biological process)
GO_res.sp <- enrichGO(gene = gene.sp, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

#make results into data frame
GO_results.sp <- as.data.frame(GO_res.sp)

#create plot base
plot <- ggplot(GO_results.sp, aes(Description, -log(pvalue)))+
  geom_bar(stat='identity')

#flip axis of plot so bars are horizontal
plot + coord_flip()


####################CORTEX ONLY analysis


#naming cortex alt IDs to output so I can filter by alt ID
col.cor.id <- col.cor$sample_id_alt

#should only be 2 group, ALS and control
unique(coldata_comp$Group)

#selecting specific columns that match spinal cord libraries. Same as col.cor.id output
#Rest of IDs are stored in desktop folder
counts.cor <- counts[,c("CGND.HRA.00115", "MANY")]

#ensure the experiments match
dim(col.cor) #575 x 30
dim(counts.cor) #20487 x 575

#do row names in coldata match column names in counts data
all(colnames(counts.cor) %in% rownames(col.cor)) #TRUE
#are they in the same order?
all(colnames(counts.cor) == rownames(col.cor)) # TRUE

#################enrichGO package pathway enrichment
library(clusterProfiler)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ggplot2)

#remove NA values
sig.cor <- na.omit(res.cor)

#filter results for the most significant
sig.cor <- sig.cor[sig.cor$padj <0.01 & sig.cor$log2FoldChange < -1]

#store gene names as a new object
gene.cor <- rownames(sig.cor)

#GO analysis, supply library, gene ID tag type, and type of GO terms (BP, biological process)
GO_res.cor <- enrichGO(gene = gene.cor, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

#make results into data frame
GO_results.cor <- as.data.frame(GO_res.cor)

#create plot base
plot <- ggplot(GO_results.cor, aes(Description, -log(pvalue)))+
  geom_bar(stat='identity')

#flip axis of plot so bars are horizontal
plot + coord_flip()


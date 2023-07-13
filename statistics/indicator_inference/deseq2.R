# test 2 DESeq2 with squeezemeta count table
# test differential gene expression -> "differential abundance" of some samples using DESeq2

library(tibble)
library(tidyr)
library(ggplot2)
library(dplyr)
library(DESeq2)

# load metadata containing sample grouping
env_eq <- read.csv("../rf_tests_sat4/metadata_batch123_dna.csv", sep=",") %>% select(-sample) %>% dplyr::rename(sample=origin) %>%
  select(sample, Rf_3) %>% filter(Rf_3=="g" | Rf_3=="b") %>% droplevels()

# load ASV table containing read counts
dna_tab <- read.csv("../stats_sat4/DNA_salmon_rare90.csv") %>% column_to_rownames("X")

# keep only metadata of samples which are in the ASV table
env_eq <- env_eq %>% filter(sample %in% rownames(dna_tab))

# load taxonomic assignment to filter only for Bacterial reads
taxo_key <- read.csv("../stats_sat4/taxo_key.csv") %>% select(-X) %>% filter(kingdom=="Bacteria")
taxojoin_bac <- dna_tab[,which(colnames(dna_tab) %in% taxo_key$TAX)]

# keep only ASVs which have at least 100 reads
taxojoin_bac <- taxojoin_bac[,which(colSums(taxojoin_bac) >99)]
taxojoin_bac <- taxojoin_bac[which(rownames(taxojoin_bac) %in% env_eq$sample),]
ints <- taxojoin_bac %>% t() %>% as.data.frame()
tab_in <- ints %>% rownames_to_column("TAX")

#DESeq2
dds <- DESeqDataSetFromMatrix(countData=tab_in, 
                              colData=env_eq, 
                              design=~Rf_3, tidy = TRUE)

# if zeros: option=  "poscounts" (calculating a modified geometric mean by taking the n-th root of the product of the non-zero counts)
dds <- DESeq(dds, sfType = "poscounts")
res <- results(dds)
head(results(dds, tidy=TRUE))
summary(res)
res <- res[order(res$padj),]
saveRDS(res, "DESeq2_results.rds")

# explore results

plotMA(res, ylim=c(-2,2))
plotCounts(dds, gene="TAX_165", intgroup="Rf_3")
head(res)
head(res[6])

# make and filter result data frame
padj_df <- as.data.frame(res[6]) 
log2fc_df <- as.data.frame(res[2])
res_df <- cbind(padj_df, log2fc_df)
sig_df <- res_df %>% filter(padj<0.05) %>% filter(log2FoldChange>0 | log2FoldChange<(-0))

sig_ups <- res_df %>% filter(padj<0.05) %>% filter(log2FoldChange<0)
write.csv(sig_ups, "DESeq2_upregulated_genes.csv")

sig_downs <- res_df %>% filter(padj<0.05) %>% filter(log2FoldChange>0)
write.csv(sig_downs, "DESeq2_downregulated_genes.csv")

saveRDS(dds, "dds.rds")

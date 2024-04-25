# indicators_deseq
# Verena Rubel 
# RPTU Kaiserslautern Landau
# original 30.10.2023 
# last revision 25.04.2024
# DESeq2 with gene input count tables
# test "differential gene expression" using DESeq2

library(tibble)
library(tidyr)
library(ggplot2)
library(dplyr)
library(DESeq2)
library(vegan)
'%!in%' <- function(x,y)!('%in%'(x,y))

# load env with info containing grouping (good/bad)
env <- read.csv2("analysis_cleaned/metadata_RNA_Table_Batch123_cleaned.csv")

# load ASV table -> here: gene table
tab_MT <- read.csv2("analysis_cleaned/RNA_Table_Batch123_cleaned.csv")

# transform first column into kegg id only and remove gene description for anaysis
tab_MT2 <- tab_MT %>% separate(X, into = c("KEGG_id", "KEGG_gene"), sep = "_", remove = TRUE)
tab_in <- tab_MT2 %>% select(-KEGG_gene) %>% column_to_rownames("KEGG_id") %>% round(0) %>%
  rownames_to_column("KEGG_id")

#DESeq2

dds <- DESeqDataSetFromMatrix(countData=tab_in, 
                              colData=env, 
                              design=~env_state, tidy = TRUE)
# if zeros: option=  "poscounts" (calculating a modified geometric mean by taking the 
# n-th root of the product of the non-zero counts)

dds <- DESeq(dds, sfType = "poscounts")
res <- results(dds)
res <- res[order(res$padj),]
write.csv(res, "results/indicators_deseq/DESeq2_output.csv")

png("results/indicators_deseq/DESeq2_results_plot.png")
plotMA(res, ylim=c(-2,2))
dev.off()

# weil reihenfolge b - g : up/down regulated umgedreht
# wir wollen wissen: transition 2 - 5 , g - b up or down?

padj_df <- as.data.frame(res[6]) 
log2fc_df <- as.data.frame(res[2])
res_df <- cbind(padj_df, log2fc_df)
sig_df <- res_df %>% filter(padj<0.05) %>% filter(log2FoldChange>0 | log2FoldChange<(-0))

# for upregulation (which is "downregulation" and therefore negative, transfrom to positive values
sig_ups <- res_df %>% filter(padj<0.05) %>% filter(log2FoldChange<0) %>% 
  abs() %>% filter(log2FoldChange>=2)
write.csv(sig_ups, "results/indicators_deseq/DESeq2_upregulated_genes.csv")
#
sig_downs <- res_df %>% filter(padj<0.05) %>% filter(log2FoldChange>0) %>% 
  abs() %>% filter(log2FoldChange>=2)
write.csv(sig_downs, "results/indicators_deseq/DESeq2_downregulated_genes.csv")

# Plot the data
p <- ggplot(sig_ups, aes(x = padj, y = log2FoldChange)) +
  geom_point() +
  labs(x = "adj. p-value", y = "Log Fold Change", title="698 significantly upregulated genes") +
  theme_light()+
  theme(panel.grid = element_blank())
p
ggsave("results/indicators_deseq/plot_upreg.pdf",p,height=4,width=5)

p2 <- ggplot(sig_downs, aes(x = padj, y = log2FoldChange)) +
  geom_point() +
  labs(x = "adj. p-value", y = "Log Fold Change", title="2443 significantly downregulated genes") +
  theme_light()+
  theme(panel.grid = element_blank())
p2
ggsave("results/indicators_deseq/plot_downreg.pdf",p2,height=4,width=5)

# family stats
# Verena Rubel 
# last edit 28.06.2024
library(tibble)
library(tidyr)
library(dplyr)
library(ggplot2)
library(vegan)
library("VennDiagram", lib.loc = "/home/vdully/Libraries")
library("gplots", lib.loc = "/home/vdully/Libraries")
'%!in%' <- function(x,y)!('%in%'(x,y))

# read number of indicators
n_inds_fam <- read.csv("inds_family_without_aldex_without_deseq_RF_raw.csv", sep=",") 
fam_res_long <- n_inds_fam %>% gather("algo", "value", -type)

fam_res_param_long <- fam_res_long  %>% mutate(MBMG = ifelse(grepl("MB", type), "MB", ifelse(grepl("MG", type), "MG", NA))) %>%
  mutate(dataset = ifelse(grepl("SCO", type), "SCO", ifelse(grepl("NOR", type), "NOR", "both")))

kruskal_test <- kruskal.test(value ~ MBMG, data = fam_res_param_long)
print(kruskal_test)

fam_res_param_long %>% select(value, MBMG) %>% filter(MBMG=="MB") %>% summary()
fam_res_param_long %>% select(value, MBMG) %>% filter(MBMG=="MG") %>% summary()

plot1<-ggplot(fam_res_param_long,aes(MBMG,value)) +
  geom_boxplot() +
  theme_light()+
  theme(panel.grid = element_blank())+
  labs(y="Number of indicator candidates", x="Sequencing strategy")+
  facet_grid(.~algo)+
  theme(
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(color = "black"))
plot1
ggsave("paper_boxplot_1.pdf",plot1,height=3,width=5)

sco_tmp <- fam_res_param_long %>% filter(dataset=="NOR")
kruskal_test <- kruskal.test(value ~ MBMG, data = sco_tmp)
print(kruskal_test)

plot1<-ggplot(fam_res_param_long,aes(MBMG,value)) +
  geom_point() +
  theme_light()+
  theme(panel.grid = element_blank())+
  labs(y="Number of indicator candidates", x="Sequencing strategy")+
  facet_grid(dataset~algo)+
  theme(
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(color = "black"))
plot1
ggsave("paper_boxplot_1B.pdf",plot1,height=3,width=5)

# read prediciton accuracies
res_fam_raw <- read.csv("results_family_typed_RF_raw.csv", sep=",") %>% select(-aldex) %>% select(-deseq2)
res_fam_raw_long <- res_fam_raw %>% gather("algo", "value", -type)
res_fam_raw_long_param <- res_fam_raw_long  %>% mutate(MBMG = ifelse(grepl("MB", type), "MB", ifelse(grepl("MG", type), "MG", NA))) %>%
  mutate(dataset = ifelse(grepl("SCO", type), "SCO", ifelse(grepl("NOR", type), "NOR", "both")))

res_fam_raw_long_param %>% select(value, MBMG) %>% filter(MBMG=="MB") %>% summary()
res_fam_raw_long_param %>% select(value, MBMG) %>% filter(MBMG=="MG") %>% summary()

kruskal_test <- kruskal.test(value ~ MBMG, data = res_fam_raw_long_param)
print(kruskal_test)

plot2<-ggplot(res_fam_raw_long_param,aes(MBMG,value)) +
  geom_boxplot() +
  theme_light()+
  theme(panel.grid = element_blank())+
  labs(y="Prediction accuracy", x="Sequencing strategy")+
  facet_grid(.~algo)+
  theme(
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(color = "black"))
plot2
ggsave("paper_boxplot_2A.pdf",plot2,height=3,width=5)

plot2<-ggplot(res_fam_raw_long_param,aes(MBMG,value)) +
  geom_point() +
  theme_light()+
  theme(panel.grid = element_blank())+
  labs(y="Prediction accuracy", x="Sequencing strategy")+
  facet_grid(dataset~algo)+
  theme(
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(color = "black"))
plot2
ggsave("paper_boxplot_2B.pdf",plot2,height=3,width=5)


# comparison of indicator candidates

# indispec

path_indispec <- "/work/TUK-CoastMon/MB_vs_MG/combined_final/results/indicators_indispec/"
indicator_set_indispec_mb <- read.csv(paste0(path_indispec, "IndicSpec_adj_pval_MB.csv")) %>% filter(p.value.bh < 0.05) %>% select(rn) %>% rename(X=rn)
indicator_set_indispec_mg <- read.csv(paste0(path_indispec, "IndicSpec_adj_pval_MG.csv")) %>% filter(p.value.bh < 0.05) %>% select(rn) %>% rename(X=rn)
indicator_set_indispec_mb_sco <- read.csv(paste0(path_indispec, "IndicSpec_adj_pval_MB_SCO2.csv")) %>% filter(p.value.bh < 0.05) %>% select(rn) %>% rename(X=rn)
indicator_set_indispec_mb_nor <- read.csv(paste0(path_indispec, "IndicSpec_adj_pval_MB_NOR2.csv")) %>% filter(p.value.bh < 0.05) %>% select(rn) %>% rename(X=rn)
indicator_set_indispec_mg_sco <- read.csv(paste0(path_indispec, "IndicSpec_adj_pval_MG_SCO2.csv")) %>% filter(p.value.bh < 0.05) %>% select(rn) %>% rename(X=rn)
indicator_set_indispec_mg_nor <- read.csv(paste0(path_indispec, "IndicSpec_adj_pval_MG_NOR2.csv")) %>% filter(p.value.bh < 0.05) %>% select(rn) %>% rename(X=rn)


path_RF <- "/work/TUK-CoastMon/MB_vs_MG/combined_final/results/indicators_RF/"

indicator_set_RF_mb <- read.csv(paste0(path_RF, "MB_classification/varimp_mean_allmodels_MB_classification_MINIMUM.csv")) %>% filter(min > 0) %>% select(X)
indicator_set_RF_mg <- read.csv(paste0(path_RF, "MG_classification/varimp_mean_allmodels_MG_classification_MINIMUM.csv")) %>% filter(min > 0) %>% select(X)
indicator_set_RF_mb_sco <- read.csv(paste0(path_RF, "MB_classification_SCO/varimp_mean_allmodels_MB_classification_SCO_MINIMUM.csv")) %>% filter(min > 0) %>% select(X)
indicator_set_RF_mb_nor <- read.csv(paste0(path_RF, "MB_classification_NOR/varimp_mean_allmodels_MB_classification_NOR_MINIMUM.csv")) %>% filter(min > 0) %>% select(X)
indicator_set_RF_mg_sco <- read.csv(paste0(path_RF, "MG_classification_SCO/varimp_mean_allmodels_MG_classification_SCO_MINIMUM.csv")) %>% filter(min > 0) %>% select(X)
indicator_set_RF_mg_nor <- read.csv(paste0(path_RF, "MG_classification_NOR/varimp_mean_allmodels_MG_classification_NOR_MINIMUM.csv")) %>% filter(min > 0) %>% select(X)

mb_tmp <- as.data.frame(all_ind2_num) %>% select(contains("MB")) %>% select(-contains("NOR")) %>% select(-contains("SCO")) %>% mutate(total=rowSums(.)) %>% filter(total>0) %>% select(-total) 
mg_tmp <- as.data.frame(all_ind2_num) %>% select(contains("MG")) %>% select(-contains("NOR")) %>% select(-contains("SCO")) %>% mutate(total=rowSums(.)) %>% filter(total>0) %>% select(-total) 
MB <- rownames(mb_tmp)
MG <- rownames(mg_tmp)

MB <- indicator_set_RF_reg_mb$X
MG <- indicator_set_RF_reg_mg$X

# Calculate intersections and unique elements
intersection <- length(dplyr::intersect(MB, MG))
only_MB <- length(setdiff(MB, MG))
only_MG <- length(setdiff(MG, MB))
#only_MB2 <- setdiff(MB, MG)
#write.csv(only_MB2, "only_MB2.csv")
#only_MG2 <- setdiff(MG, MB)
#write.csv(only_MG2, "only_MG2.csv")
total_MB <- length(MB)
total_MG <- length(MG)

# Calculate the total union
total_union <- length(union(MB, MG))

# Calculate percentages relative to the total union
percentage_only_MB <- only_MB / total_union * 100
percentage_only_MG <- only_MG / total_union * 100
percentage_intersection <- intersection / total_union * 100

# Plot the Venn diagram manually
plot(0, 0, type = "n", xlim = c(-1, 3), ylim = c(-1, 3), xlab = "", ylab = "", xaxt = 'n', yaxt = 'n')

# Draw circles
symbols(0, 1, circles = 1, inches = FALSE, add = TRUE, bg = alpha("blue", 0.1), fg = NULL)
symbols(1, 1, circles = 1, inches = FALSE, add = TRUE, bg = alpha("yellow", 0.1), fg = NULL)

# Add text labels with correct percentages
text(-0.5, 1, paste("MB", only_MB, "\n", round(percentage_only_MB, 2), "%"))
text(1.5, 1, paste("MG", only_MG, "\n", round(percentage_only_MG, 2), "%"))
text(0.5, 1, paste("", intersection, "\n", round(percentage_intersection, 2), "%"))


# heatmap for identity of indicators: only indispec RF datasets MB MG
# get all indicators and mapping
# make binary

# make list of all families
mb_tab <- rownames(read.csv("../data_created/reads_fam_MB_clean.csv", row.names = 1))
mg_tab <- rownames(read.csv("../data_created/reads_fam_MG_clean.csv", row.names = 1))
fams <- as.data.frame(unique(c(mb_tab, mg_tab)))
colnames(fams) <- "Family"
print(fams)

# indispec

path_indispec <- "/work/TUK-CoastMon/MB_vs_MG/combined_final/results/indicators_indispec/"

indicator_set_indispec_mb <- read.csv(paste0(path_indispec, "IndicSpec_adj_pval_MB.csv")) %>% filter(p.value.bh < 0.05) %>% select(rn) %>% rename(X=rn) %>% mutate(indic_MB=1)%>% rename(Family=X)
indicator_set_indispec_mg <- read.csv(paste0(path_indispec, "IndicSpec_adj_pval_MG.csv")) %>% filter(p.value.bh < 0.05) %>% select(rn) %>% rename(X=rn) %>% mutate(indic_MG=1)%>% rename(Family=X)
indicator_set_indispec_mb_sco <- read.csv(paste0(path_indispec, "IndicSpec_adj_pval_MB_SCO2.csv")) %>% filter(p.value.bh < 0.05) %>% select(rn) %>% rename(X=rn) %>% mutate(indic_MB_SCO=1)%>% rename(Family=X)
indicator_set_indispec_mb_nor <- read.csv(paste0(path_indispec, "IndicSpec_adj_pval_MB_NOR2.csv")) %>% filter(p.value.bh < 0.05) %>% select(rn) %>% rename(X=rn) %>% mutate(indic_MB_NOR=1)%>% rename(Family=X)
indicator_set_indispec_mg_sco <- read.csv(paste0(path_indispec, "IndicSpec_adj_pval_MG_SCO2.csv")) %>% filter(p.value.bh < 0.05) %>% select(rn) %>% rename(X=rn) %>% mutate(indic_MG_SCO=1)%>% rename(Family=X)
indicator_set_indispec_mg_nor <- read.csv(paste0(path_indispec, "IndicSpec_adj_pval_MG_NOR2.csv")) %>% filter(p.value.bh < 0.05) %>% select(rn) %>% rename(X=rn) %>% mutate(indic_MG_NOR=1)%>% rename(Family=X)


path_RF <- "/work/TUK-CoastMon/MB_vs_MG/combined_final/results/indicators_RF/"

indicator_set_RF_mb <- read.csv(paste0(path_RF, "MB_classification/varimp_mean_allmodels_MB_classification_MINIMUM.csv")) %>% filter(min > 0) %>% select(X) %>% mutate(rf_class_MB=1)%>% rename(Family=X)
indicator_set_RF_mg <- read.csv(paste0(path_RF, "MG_classification/varimp_mean_allmodels_MG_classification_MINIMUM.csv")) %>% filter(min > 0) %>% select(X) %>% mutate(rf_class_MG=1)%>% rename(Family=X)
indicator_set_RF_mb_sco <- read.csv(paste0(path_RF, "MB_classification_SCO/varimp_mean_allmodels_MB_classification_SCO_MINIMUM.csv")) %>% filter(min > 0) %>% select(X)%>% mutate(rf_class_MB_SCO=1)%>% rename(Family=X)
indicator_set_RF_mb_nor <- read.csv(paste0(path_RF, "MB_classification_NOR/varimp_mean_allmodels_MB_classification_NOR_MINIMUM.csv")) %>% filter(min > 0) %>% select(X)%>% mutate(rf_class_MB_NOR=1)%>% rename(Family=X)
indicator_set_RF_mg_sco <- read.csv(paste0(path_RF, "MG_classification_SCO/varimp_mean_allmodels_MG_classification_SCO_MINIMUM.csv")) %>% filter(min > 0) %>% select(X)%>% mutate(rf_class_MG_SCO=1)%>% rename(Family=X)
indicator_set_RF_mg_nor <- read.csv(paste0(path_RF, "MG_classification_NOR/varimp_mean_allmodels_MG_classification_NOR_MINIMUM.csv")) %>% filter(min > 0) %>% select(X)%>% mutate(rf_class_MG_NOR=1)%>% rename(Family=X)


inidc_ind <- fams %>% left_join(indicator_set_indispec_mb) %>% left_join(indicator_set_indispec_mg) %>% left_join(indicator_set_indispec_mb_sco) %>% left_join(indicator_set_indispec_mb_nor) %>% left_join(indicator_set_indispec_mg_sco) %>% left_join(indicator_set_indispec_mg_nor)
rf_class_ind <- fams %>% left_join(indicator_set_RF_mb) %>% left_join(indicator_set_RF_mg) %>% left_join(indicator_set_RF_mb_sco) %>% left_join(indicator_set_RF_mb_nor) %>% left_join(indicator_set_RF_mg_sco) %>% left_join(indicator_set_RF_mg_nor)
rf_regr_ind <- fams %>% left_join(indicator_set_RF_reg_mb) %>% left_join(indicator_set_RF_reg_mg) %>% left_join(indicator_set_RF_reg_mb_sco) %>% left_join(indicator_set_RF_reg_mb_nor) %>% left_join(indicator_set_RF_reg_mg_sco) %>% left_join(indicator_set_RF_reg_mg_nor)

all_ind <- inidc_ind %>%  left_join(rf_class_ind) %>% left_join(rf_regr_ind)
all_ind[is.na(all_ind)] <- 0
all_ind2 <- all_ind %>% column_to_rownames("Family") %>%  mutate(total=rowSums(.)) %>% arrange(-total)
write.csv(all_ind2, "all_indicators_binary_3algos.csv")

mapping <- as.data.frame(table(colnames(all_ind2), dnn = list("origin"))) %>% filter(origin!="total") %>% select(-Freq)  %>%
  mutate(method = ifelse(grepl("MB", origin), "MB", "MG")) %>% 
  mutate(
    method = ifelse(grepl("MB", origin), "MB", "MG"),
    cty = case_when(
      grepl("NOR", origin) ~ "NOR",
      grepl("SCO", origin) ~ "SCO",
      TRUE ~ "NOR+SCO"
    ),
    indic_method = case_when(
      grepl("indic", origin) ~ "indic",
      grepl("rf_class", origin) ~ "rf_class",
      grepl("rf_regr", origin) ~ "rf_regr",
      TRUE ~ NA_character_  # or specify another default value if needed
    )
  )
write.csv(mapping, "mapping_methods_3algos.csv")

mapping <- read.csv("mapping_methods_3algos.csv")
all_ind <- read.csv("all_indicators_binary_3algos.csv", row.names = 1) %>% select(-total)
all_ind_filt <- all_ind %>% select(-contains("rf_regr")) %>% droplevels()
all_ind2 <- all_ind_filt[rowSums(all_ind_filt)>0, ]
all_ind2_num <- as.matrix(all_ind2[sapply(all_ind2, is.numeric)])
colnames_rf_new <- gsub("rf_class_", "Variable Importance ", colnames(all_ind2_num))
colnames_rf_new <- gsub("indic_", "IndicSpecies ", colnames_rf_new)
colnames_rf_new <- gsub("MG_NOR", "NOR MG", colnames_rf_new)
colnames_rf_new <- gsub("MG_SCO", "SCO MG", colnames_rf_new)
colnames_rf_new <- gsub("MB_NOR", "NOR MB", colnames_rf_new)
colnames_rf_new <- gsub("MB_SCO", "SCO MB", colnames_rf_new)
print(colnames_rf_new)
colnames_rf_new[1] <- "IndicSpecies NOR/SCO MB"
colnames_rf_new[2] <- "IndicSpecies NOR/SCO MG"
colnames_rf_new[7] <- "Variable Importance NOR/SCO MB"
colnames_rf_new[8] <- "Variable Importance NOR/SCO MG"

colnames(all_ind2_num) <- colnames_rf_new



# Specify the custom column order
custom_col_order <- c("IndicSpecies SCO MB", "Variable Importance SCO MB", "IndicSpecies NOR MB", "Variable Importance NOR MB",
                      "IndicSpecies NOR/SCO MB", "Variable Importance NOR/SCO MB", 
                      "IndicSpecies SCO MG", "Variable Importance SCO MG", "IndicSpecies NOR MG", "Variable Importance NOR MG",
                      "IndicSpecies NOR/SCO MG", "Variable Importance NOR/SCO MG")

# Reorder the columns in the matrix
all_ind2_num <- all_ind2_num[, custom_col_order]

# Calculate the row order based on the frequency of values in the first column
row_order <- rownames(all_ind2_num)

# Use the heatmap.2 function with the specified row order
pdf("heatmap_output_z.pdf", width=8, height=12)
heatmap.2(all_ind2_num[row_order, ], 
          col = c("white", "black"), 
          colsep=0:ncol(all_ind2_num),
          rowsep=0:nrow(all_ind2_num),
          sepcolor="grey",        # Grid color
          sepwidth=c(0.01, 0.01), # Thickness of the separators  
          Rowv=NA, 
          Colv=NA,                # Prevent clustering and enforce order
          dendrogram = "none", 
          trace="none", 
          key="none", 
          margins = c(10, 10),    # Adjust margins to give space at the top
          cexCol=0.8,             # Control size of column labels
          srtCol=45,              # Rotate the column labels (adjust angle as needed)
          adjCol = c(0.5, 7),     # Adjust position of column labels
          offsetCol = 5        # Offset labels to move them to the top
)
dev.off()


# found families (all,. not only indicators)
# read table and order columns alphabetically!
# load env with info containing grouping (good/bad)
out_sams <- read.csv("../samples_out_NMDS_MB_MG.csv")
env0 <- read.csv("../data/metadata_MB_MG.csv")
env0 <- env0[order(env0$Sample),]
env <- env0[which(env0$Sample %!in% out_sams$S3), ]
env_cty <- env %>% dplyr::select(Sample, Country)
env_ambi <- env %>% dplyr::select(Sample, AMBI)
env_eqform <- env %>% mutate(EQ_class = ifelse(EQ %in% c(1, 2), "good", "bad")) %>% dplyr::select(Sample, EQ_class)
sams_nor <- env %>% filter(Country=="NOR")
sams_sco <- env %>% filter(Country=="SCO")

tab_MB <- read.csv("../data_created/reads_fam_MB_clean.csv", row.names = 1)
tab_MB <- tab_MB %>% dplyr::select(order(colnames(tab_MB)))
tab_MG <- read.csv("../data_created/reads_fam_MG_clean.csv", row.names = 1) 
tab_MG <- tab_MG %>% dplyr::select(order(colnames(tab_MG)))


mb_tmp <- tab_MB[,which(colnames(tab_MB) %in% sams_sco$Sample)] %>%  mutate(total=rowSums(.)) %>% filter(total>0) %>% select(-total)
mg_tmp <- tab_MG[,which(colnames(tab_MG) %in% sams_sco$Sample)] %>%  mutate(total=rowSums(.)) %>% filter(total>0) %>% select(-total)


## next task: find MG indicators only
# are they strong?
# are they universal?

# start: 43 universal indicators of MG/MB
only_MG2 <- read.csv("only_MG2.csv", row.names = 1)
only_MB2 <- read.csv("only_MB2.csv", row.names = 1)
# question: are they strong?
# maybe answer per dataset?

# rf repred results. need: column in df mapping the dataset. 
# make boxplot MB indicators MG indicators per data set
# read indicators

#indispec

path <- "/work/TUK-CoastMon/MB_vs_MG/combined_final/reprediction_RF/results/family/indispec/MG/"
temp = list.files(path, pattern="*rds")
temp2 <- paste0(path, temp)
sams <- sapply(strsplit(basename(temp), ".rds"), `[`, 1)
sams2 <- sapply(strsplit(basename(sams), "model_acc_"), `[`, 2)
for (i in 1:length(temp2)) assign(sams2[i], readRDS(temp2[i]))
path <- "/work/TUK-CoastMon/MB_vs_MG/combined_final/reprediction_RF/results/family/indispec/MB/"
temp = list.files(path, pattern="*rds")
temp2 <- paste0(path, temp)
sams <- sapply(strsplit(basename(temp), ".rds"), `[`, 1)
sams2 <- sapply(strsplit(basename(sams), "model_acc_"), `[`, 2)
for (i in 1:length(temp2)) assign(sams2[i], readRDS(temp2[i]))

tmp_1 <- importance(indicator_set_indispec_mg_nor) %>% as.data.frame() %>% rownames_to_column("Family") %>% 
  select(Family, MeanDecreaseAccuracy) %>% mutate(algo="indispec") %>% mutate(MBMG="MG") %>% mutate(dataset="NOR")
tmp_2 <- importance(indicator_set_indispec_mg_sco) %>% as.data.frame() %>% rownames_to_column("Family") %>% 
  select(Family, MeanDecreaseAccuracy) %>% mutate(algo="indispec") %>% mutate(MBMG="MG") %>% mutate(dataset="SCO")
tmp_3 <- importance(indicator_set_indispec_mg) %>% as.data.frame() %>% rownames_to_column("Family") %>% 
  select(Family, MeanDecreaseAccuracy) %>% mutate(algo="indispec") %>% mutate(MBMG="MG") %>% mutate(dataset="both")

tmp_4 <- importance(indicator_set_indispec_mb_nor) %>% as.data.frame() %>% rownames_to_column("Family") %>% 
  select(Family, MeanDecreaseAccuracy) %>% mutate(algo="indispec") %>% mutate(MBMG="MB") %>% mutate(dataset="NOR")
tmp_5 <- importance(indicator_set_indispec_mb_sco) %>% as.data.frame() %>% rownames_to_column("Family") %>% 
  select(Family, MeanDecreaseAccuracy) %>% mutate(algo="indispec") %>% mutate(MBMG="MB") %>% mutate(dataset="SCO")
tmp_6 <- importance(indicator_set_indispec_mb) %>% as.data.frame() %>% rownames_to_column("Family") %>% 
  select(Family, MeanDecreaseAccuracy) %>% mutate(algo="indispec") %>% mutate(MBMG="MB") %>% mutate(dataset="both")

indispec_varimps <- rbind(tmp_1, tmp_2, tmp_3, tmp_4, tmp_5, tmp_6)
write.csv(indispec_varimps, "results/asvs/indispec/indispec_varimps_combined.csv")
indispec_varimps <- read.csv("results/asvs/indispec/indispec_varimps_combined.csv", row.names = 1)

indispec_varimps_set <- indispec_varimps %>% 
  mutate(set=ifelse(Family %in% only_MG2$x, "MG only", ifelse(Family %in% only_MB2$x, "MB only", "both")))

# RF classification

path <- "/work/TUK-CoastMon/MB_vs_MG/combined_final/reprediction_RF/results/family/RF_classification/MG/"
temp = list.files(path, pattern="*rds")
temp2 <- paste0(path, temp)
sams <- sapply(strsplit(basename(temp), ".rds"), `[`, 1)
sams2 <- sapply(strsplit(basename(sams), "model_acc_"), `[`, 2)
for (i in 1:length(temp2)) assign(sams2[i], readRDS(temp2[i]))
path <- "/work/TUK-CoastMon/MB_vs_MG/combined_final/reprediction_RF/results/family/RF_classification/MB/"
temp = list.files(path, pattern="*rds")
temp2 <- paste0(path, temp)
sams <- sapply(strsplit(basename(temp), ".rds"), `[`, 1)
sams2 <- sapply(strsplit(basename(sams), "model_acc_"), `[`, 2)
for (i in 1:length(temp2)) assign(sams2[i], readRDS(temp2[i]))

tmp_1 <- importance(indicator_set_indispec_mg_nor) %>% as.data.frame() %>% rownames_to_column("Family") %>% 
  select(Family, MeanDecreaseAccuracy) %>% mutate(algo="RF_classification") %>% mutate(MBMG="MG") %>% mutate(dataset="NOR")
tmp_2 <- importance(indicator_set_indispec_mg_sco) %>% as.data.frame() %>% rownames_to_column("Family") %>% 
  select(Family, MeanDecreaseAccuracy) %>% mutate(algo="RF_classification") %>% mutate(MBMG="MG") %>% mutate(dataset="SCO")
tmp_3 <- importance(indicator_set_indispec_mg) %>% as.data.frame() %>% rownames_to_column("Family") %>% 
  select(Family, MeanDecreaseAccuracy) %>% mutate(algo="RF_classification") %>% mutate(MBMG="MG") %>% mutate(dataset="both")

tmp_4 <- importance(indicator_set_indispec_mb_nor) %>% as.data.frame() %>% rownames_to_column("Family") %>% 
  select(Family, MeanDecreaseAccuracy) %>% mutate(algo="RF_classification") %>% mutate(MBMG="MB") %>% mutate(dataset="NOR")
tmp_5 <- importance(indicator_set_indispec_mb_sco) %>% as.data.frame() %>% rownames_to_column("Family") %>% 
  select(Family, MeanDecreaseAccuracy) %>% mutate(algo="RF_classification") %>% mutate(MBMG="MB") %>% mutate(dataset="SCO")
tmp_6 <- importance(indicator_set_indispec_mb) %>% as.data.frame() %>% rownames_to_column("Family") %>% 
  select(Family, MeanDecreaseAccuracy) %>% mutate(algo="RF_classification") %>% mutate(MBMG="MB") %>% mutate(dataset="both")

rf_class_varimps <- rbind(tmp_1, tmp_2, tmp_3, tmp_4, tmp_5, tmp_6)
write.csv(rf_class_varimps, "results/asvs/RF_classification/rf_class_varimps_combined.csv")

rf_class_varimps_set <- rf_class_varimps %>% 
  mutate(set=ifelse(Family %in% only_MG2$x, "MG only", ifelse(Family %in% only_MB2$x, "MB only", "both")))

# plot

varimp_set <- rbind(indispec_varimps_set, rf_class_varimps_set, rf_regr_varimps_set)

plot_varimp_indic<-ggplot(varimp_set,aes(set,MeanDecreaseAccuracy)) +
  geom_boxplot() +
  theme_light()+
  theme(panel.grid = element_blank())+
  labs(y="mean variable importance", x="indicator occurence")+
  facet_grid(.~algo)+
  theme(
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(color = "black"))
plot_varimp_indic
ggsave("paper_boxplot_varimp_indic_A.pdf",plot_varimp_indic,height=3,width=6)

plot_varimp_indic<-ggplot(varimp_set,aes(set,MeanDecreaseAccuracy)) +
  geom_boxplot() +
  theme_light()+
  theme(panel.grid = element_blank())+
  labs(y="mean variable importance", x="indicator occurence")+
  facet_grid(algo~dataset)+
  theme(
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(color = "black"))
plot_varimp_indic
ggsave("paper_boxplot_varimp_indic_B.pdf",plot_varimp_indic,height=5,width=7)

varimp_set2 <- varimp_set %>% filter(set!="both")
kruskal_test <- kruskal.test(MeanDecreaseAccuracy ~ set, data = varimp_set2)
print(kruskal_test)

varimp_set2 %>% select(MeanDecreaseAccuracy, set) %>% filter(set=="MB only") %>% summary()
varimp_set2 %>% select(MeanDecreaseAccuracy, set) %>% filter(set=="MG only") %>% summary()


## make pred acc plot families best possible
preds <- read.csv("results_family_typed.csv") %>% select(type, indispec, RF_classification, RF_regression) %>% 
  gather("algo", "value", -type) %>% mutate(MBMG = ifelse(grepl("MB", type), "MB", ifelse(grepl("MG", type), "MG", NA))) %>%
  mutate(dataset = ifelse(grepl("SCO", type), "SCO", ifelse(grepl("NOR", type), "NOR", "both")))

plot_pred<-ggplot(preds,aes(MBMG,value)) +
  geom_boxplot() +
  theme_light()+
  theme(panel.grid = element_blank())+
  labs(y="Impact prediction accuracy", x="Sequencing strategy")+
  facet_grid(.~algo)+
  theme(
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(color = "black"))
plot_pred
ggsave("paper_boxplot_pred.pdf",plot_pred,height=3,width=5)

preds %>% filter(MBMG=="MB") %>% select(value) %>% summarise_each(funs(mean))
preds %>% filter(MBMG=="MG") %>% select(value) %>% summarise_each(funs(mean))


kruskal_test <- kruskal.test(value ~ MBMG, data = preds)
print(kruskal_test)

plot_pred2<-ggplot(preds,aes(MBMG,value)) +
  geom_point() +
  theme_light()+
  theme(panel.grid = element_blank())+
  labs(y="Impact prediction accuracy", x="Sequencing strategy")+
  facet_grid(dataset~algo)+
  theme(
    strip.background = element_rect(fill = "white", color = "white"),
    strip.text = element_text(color = "black"))
plot_pred2
ggsave("paper_boxplot2_pred.pdf",plot_pred2,height=3,width=5)

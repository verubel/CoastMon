# rf bac
# rf3: EQ good and EQ bad only 
library(ggplot2)
library(tibble)
library(randomForest)
library(tidyr)
library(dplyr)

# load metadata
env_eq <- read.csv("metadata_batch123_dna.csv", sep=",") %>% select(-sample) %>% rename(sample=origin) %>%
  select(sample, Rf_3) %>% na.omit() %>% droplevels()


# load ASV table
dna_tab <- read.csv("../stats_sat4/DNA_salmon_rare90.csv") %>% column_to_rownames("X")

# load taxonomic assignment
taxo_key <- read.csv("../stats_sat4/taxo_key.csv") %>% select(-X) %>% filter(kingdom=="Bacteria")

# keep only bacterial reads
taxojoin_bac <- dna_tab[,which(colnames(dna_tab) %in% taxo_key$TAX)]

# keep only ASVs with at least 100 reads
taxojoin_bac <- taxojoin_bac[,which(colSums(taxojoin_bac) >99)]
taxojoin_bac <- taxojoin_bac[which(rownames(taxojoin_bac) %in% env_eq$sample),]
rf_in <- taxojoin_bac %>% rownames_to_column("sample") %>% left_join(env_eq) %>% column_to_rownames("sample")

###### RF pilot based on TAX per EQ
# specify regression or classification
# rf_3: classification! factors
rf_in$Rf_3 <- as.factor(rf_in$Rf_3)


# pilot model
set.seed(1)
pilot_model1 <- randomForest(Rf_3 ~ .,
                             data = rf_in,
                             ntree = 5000,
                             na.action = "na.omit",
                             importance=TRUE)
pilot_model1
varImpPlot(pilot_model1, type=1, n=100)
imp_df <- importance(pilot_model1) %>% as.data.frame() %>% select(MeanDecreaseAccuracy) %>% filter(MeanDecreaseAccuracy>0)
val1<- (summary(imp_df$MeanDecreaseAccuracy))[5]
imp_candidates <- imp_df %>% filter(MeanDecreaseAccuracy>val1)

# make 10 models
# average MDS per model
# them make calculation and indicator assignment

set.seed(12345)
random_numbers <- sample(1:10000, 10)
runs <- paste0("run", 1:10)

for (i in 1:length(random_numbers)) {
  set.seed(random_numbers[i])
  pilot_model <- randomForest(Rf_3 ~ ., 
                              data = rf_in,
                              ntree = 5000, 
                              na.action = "na.omit", 
                              importance=TRUE)
  print(pilot_model)
  imp_df <- importance(pilot_model) %>% as.data.frame() 
  saveRDS(pilot_model, paste0("rf_3/model_run", i, ".rds"))
  write.csv(imp_df, paste0("rf_3/varimp_run", i, ".csv"))
}


# get varimp for all models for each ASV
path <- "/work/TUK-CoastMon/stats_analysis_test4/DNA/rf_tests_sat4/rf_3/"
temp = list.files(path, pattern="*.csv")
temp2 <- paste0(path, temp)
sams <- sapply(strsplit(basename(temp), "*.csv"), `[`, 1)
for (i in 1:length(temp2)) assign(sams[i], read.csv(temp2[i]))

for (i in 1:length(sams)) {
  tab <- get(sams[i]) %>% select(X, MeanDecreaseAccuracy)
  colnames(tab) <- c("TAX", sams[i])
  assign(paste0(sams[i], "tab"), tab) 
}


# join all varimp_run1tab to varimp_run10tab
varimp_list <- list()
for (i in 1:10) {
  run_tab <- paste0("varimp_run", i, "tab")
  varimp_list[[i]] <- get(run_tab)
}
# join the data frames using left_join
varimp_df <- varimp_list[[1]]
for (i in 2:10) {
  varimp_df <- left_join(varimp_df, varimp_list[[i]])
}
rownames(varimp_df) <- NULL
varimp_df <- varimp_df %>% column_to_rownames("TAX")
#write.csv(varimp_df, "rf_3/varimp_pro_model_rf3.csv")


varimp_avg <- rowMeans(varimp_df) %>% 
  as.data.frame() %>% rownames_to_column("TAX")

varimp_plot_in <- varimp_df %>% rownames_to_column("TAX") %>%
  as.data.frame() %>% gather("model", "varimp", -TAX)

varimps_top50 <- as.data.frame(rowMeans(varimp_df)) %>%
  arrange(desc(rowMeans(varimp_df))) %>%
  top_n(50) %>%
  rownames_to_column("TAX")

varimp_df_filt50 <- varimp_plot_in %>% filter(TAX %in% varimps_top50$TAX) %>% droplevels()

plot_varimps2 <- ggplot(varimp_df_filt50, aes(x=reorder(TAX,varimp), y=varimp)) + 
  geom_boxplot() +
  theme_light() +
  theme(axis.text.x = element_text(angle=0), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+
  labs(y="Mean Variable Importance ", x="taxon_id")+
  coord_flip()
plot_varimps2
#ggsave("rf_3/plot_varimp_boxplot_rf3.pdf", plot_varimps2, width=5, height =7)

varimp_df_filt_avg <- varimp_avg %>% filter(TAX %in% varimps_top50$TAX) %>% droplevels()

plot_varimps <- ggplot(varimp_df_filt_avg, aes(x=reorder(TAX,.), y=.)) + 
  geom_bar(stat="identity") +
  theme_light() +
  theme(axis.text.x = element_text(angle=0), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+
  labs(y="Mean Variable Importance ", x="taxon_id")+
  coord_flip()
plot_varimps
#ggsave("rf_3/plot_varimp_avg_rf3.pdf", plot_varimps, width=5, height =7)


# define indicators
varimp_no_neg <- varimp_avg %>% filter(.>0)
val1<- (summary(varimp_no_neg$.))[5]
imp_candidates <- varimp_no_neg %>% filter(.>val1)
dim(imp_candidates)
dim(imp_candidates)[1]/dim(varimp_no_neg)[1]

#write.csv(imp_candidates, "rf_3/imp_candidates_rf3.csv")


##all inds
varimps_top50 <- as.data.frame(rowMeans(varimp_df)) %>%
  arrange(desc(rowMeans(varimp_df))) %>%
  rownames_to_column("TAX")

varimp_df_filt50 <- varimp_plot_in %>% filter(TAX %in% varimps_top50$TAX) %>% droplevels()
varimp_df_filt_avg <- varimp_avg %>% filter(TAX %in% varimps_top50$TAX) %>% droplevels()
plot_varimps <- ggplot(varimp_df_filt_avg, aes(x=reorder(TAX,.), y=./2)) + 
  geom_bar(stat="identity") +
  theme_light() +
  theme(axis.text.x = element_text(angle=0), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+
  labs(y="Mean Variable Importance ", x="taxon_id")+
  coord_flip()
plot_varimps
#ggsave("rf_3/plot_varimp_avg_rf3_allvars.pdf", plot_varimps, width=5, height =7)


varimp_df_filt_avg <- varimp_df_filt_avg %>% filter(.>0)
# analyze the distribution
hist(varimp_df_filt_avg$., breaks = "FD", main = "Permutation Importance")

# estimate the probability density function using kernel density estimation
density <- density(varimp_df_filt_avg$.)

# find local minima in the density estimate
minima <- density$x[which(diff(sign(diff(density$y))) > 0)]

# select the threshold based on the point where the density drops significantly
# choose the first local minimum as the threshold
threshold <- minima[1]  

# plot the density estimate and mark the threshold
plot(density, main = "Density Estimation")
abline(v = threshold, col = "red")

varimp_df_filt_avg <- varimp_df_filt_avg %>% filter(.>threshold)
write.csv(varimp_df_filt_avg, "local_min_indidc_rf3.csv")

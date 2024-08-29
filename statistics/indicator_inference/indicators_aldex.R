# indicators_aldex.R
# Verena Rubel 
# RPTU Kaiserslautern Landau
# 30.10.2023 revised 28.05.2024

library(ALDEx2)
library(plyr)
library(vegan)
library(dplyr)
'%!in%' <- function(x,y)!('%in%'(x,y))

# read table and order columns alphabetically!
tab_MB <- read.csv("../data_created/ASV_Table_MB_10_targets_NONRARE.csv", row.names = 1)[1:65]
tab_MB<- tab_MB %>% dplyr::select(order(colnames(tab_MB)))

# load env with info containing grouping (good/bad)
out_sams <- read.csv("../samples_out_NMDS_MB_MG.csv")
env0 <- read.csv("../data/metadata_MB_MG.csv")
env0 <- env0[order(env0$Sample),]
env <- env0[which(env0$Sample %!in% out_sams$S3), ]
env_cty <- env %>% dplyr::select(Sample, Country)
env_ambi <- env %>% dplyr::select(Sample, AMBI)
env_eqform <- env %>% mutate(EQ_class = ifelse(EQ %in% c(1, 2), "good", "bad")) %>% dplyr::select(Sample, EQ_class)


# MB NOR
env_class_NOR <- env %>% filter(Country=="NOR") %>% 
  mutate(EQ_class = ifelse(EQ %in% c(1, 2), "good", "bad")) %>% dplyr::select(Sample, EQ_class)
tab_MB_NOR <- tab_MB[,which(colnames(tab_MB) %in% env_class_NOR$Sample)]
# load ASV table -> here: read count data per Family level
groups <- env_class_NOR$EQ_class
x.all <- aldex(tab_MB_NOR, groups, mc.samples=128, test="t", effect=TRUE,
               include.sample.summary=FALSE, denom="all", verbose=FALSE) #if more than 2 groups use kw instead of t
png("../results/ASV_indicators_aldex/ALDex2_res_MB_NOR.png")
aldex.plot(x.all, type="MA", test="wilcox", xlab="abundance", ylab="Difference")
dev.off()
boxplot(x.all$wi.eBH)
x.all %>% filter(wi.eBH<0.05)
dim(x.all %>% filter(wi.eBH<0.05))[1]
x.all %>% filter(wi.ep<0.05)
dim(x.all %>% filter(wi.ep<0.05))[1]
write.csv(x.all,"../results/ASV_indicators_aldex/ALDex2_res_MB_NOR.csv", row.names = TRUE)

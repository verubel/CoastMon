# indicators_indispec.R
# Verena Rubel 
# RPTU Kaiserslautern Landau
# 30.10.2023 revised 28.05.2024

library(indispecies)
library(tibble)
library(plyr)
library(dplyr)
library(funrar)
library(vegan)
library(data.table)
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
indval = indicspecies::multipatt(t(tab_MB_NOR), groups, control = how(nperm = 9999))
summary(indval, indvalcomp=TRUE)
#extract table of stats
library(data.table)
indisp.sign<-data.table::as.data.table(indval$sign, keep.rownames=TRUE)
#add adjusted p-value
indisp.sign[ ,p.value.bh:=p.adjust(p.value, method="BH")]
#now you can select only the indicators with adjusted significant p-values
indicspec2 <- indisp.sign[p.value.bh<=0.05, ]
dim(indicspec2)
indicspec3 <- indisp.sign[p.value<=0.05, ]
dim(indicspec3)
write.csv(indisp.sign,"../results/ASV_indicators_indispec/IndicSpec_adj_pval_MB_NOR2.csv", row.names = TRUE)

# indicators_QRS
# Verena Rubel 
# RPTU Kaiserslautern Landau
# 30.10.2023 revised 10.12.2023

library(tidyverse)
library(splines)
library(quantreg)
library(knitr)
library(pracma)
library(ggpubr)
library(phyloseq)
library(broom)
library(ggfortify)
library(tidyr)
library(tibble)
library(ggplot2)
library(quantreg)
library(dplyr)
'%!in%' <- function(x,y)!('%in%'(x,y))

# load env with info containing grouping (good/bad)
out_sams <- read.csv("samples_out_NMDS_MB_MG.csv")
env0 <- read.csv("data/metadata_MB_MG.csv")
env0 <- env0[order(env0$Sample),]
env <- env0[which(env0$Sample %!in% out_sams$S3), ]
env_cty <- env %>% dplyr::select(Sample, Country)
env_ambi <- env %>% dplyr::select(Sample, AMBI)
env_eqform <- env %>% mutate(EQ_class = ifelse(EQ %in% c(1, 2), "good", "bad")) %>% dplyr::select(Sample, EQ_class)
nor_sams<- env %>% filter(Country=="NOR") %>% droplevels()
sco_sams<- env %>% filter(Country=="SCO") %>% droplevels()

# MB both ctys
mb_tab <- read.csv("data_created/reads_fam_MB_clean.csv", row.names = 1)

input_plot2 <- mb_tab %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>% gather("Family", "relab", -Sample) %>%
  left_join(env_ambi)
p2 <-
  ggplot(input_plot2, aes(x = AMBI, y = relab)) +
  geom_point() +
  stat_quantile(method = rq,
                formula = y ~ bs(x, df = 3),
                quantiles = 0.95) +
  facet_wrap(~ Family, scales = 'free')
ggsave("results/indicators_QRS/MB/all_families_MB_splines.pdf",
       p2, width = 15, height = 10)

# get splines data from plot
gb_MB <- ggplot_build(p2) #takes 10 minutes


gb1 <-
  gb_MB$data[[2]] %>% 
  data.frame() %>%
  group_by(group, PANEL) %>%
  nest() %>%
  mutate(
    peak = map(data, ~ .$x[which.max(.$y)]),
    find_peaks = map(
      data,
      ~ findpeaks(
        .$y,
        nups = 0,
        ndowns = 0 #,
        #minpeakheight = max(.$y) / 2.4
      )
    ),
    n_peaks = map_dbl(find_peaks, nrow)
  )

plot_id <-
  gb_MB$plot[[1]] %>%
  distinct(Family) %>%
  arrange(Family)

# get the number of peaks per ASVS and farm 
npeaks <-
  gb1 %>%
  select(n_peaks) %>%
  unnest(n_peaks) %>%
  bind_cols(plot_id)

abund_peaks <-
  npeaks %>%
  group_by(Family) %>%
  summarise_at(vars(n_peaks), list(sum), na.rm = T)

peak_dat <-
  gb1 %>%
  ungroup() %>%
  select(peak) %>%
  unnest(peak) %>%
  bind_cols(plot_id)
write.csv2(peak_dat, "results/indicators_QRS/MB/peak_dat_MB.csv")

# make qrs plots for all families

input_plot1 <- mb_tab %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>% 
  left_join(env_ambi)

for (i in 1:length(rownames(mb_tab))) {
  print(rownames(mb_tab)[i])
  input<- input_plot1 %>% select(Sample, rownames(mb_tab)[i], AMBI)
  colnames(input)[2] <- "relab"
  peak_dat2 <- peak_dat %>% filter(Family==rownames(mb_tab)[i])
  p1 <-
    ggplot(input, aes(x = AMBI, y = sqrt(relab))) +
    geom_point()+
    stat_quantile(method = rq,
                  formula = y ~ bs(x, df = 3),
                  quantiles = 0.95) +
    geom_vline(aes(xintercept = peak_dat2$peak)) +
    #facet_wrap(~ ASV, scales = 'free') +
    labs(y = 'Relative read abundance (sqrt)',
         x = 'AMBI' ,
         parse = T) +
    theme_bw(base_size = 12, base_family ="") %+replace%
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      strip.background = element_blank()
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0,6)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0,NA))
  p1 <- p1+ ggtitle(rownames(mb_tab)[i])
  p1
  ggsave(paste0("results/indicators_QRS/MB/plot_", rownames(mb_tab)[i], ".pdf"),
         p1, width = 6, height = 4.5)
}

#
#
# NOR MB

mb_tab <- read.csv("data_created/reads_fam_MB_clean.csv", row.names = 1)
mb_tab <- mb_tab[,which(colnames(mb_tab) %in% nor_sams$Sample)]


input_plot2 <- mb_tab %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>% gather("Family", "relab", -Sample) %>%
  left_join(env_ambi)
p2 <-
  ggplot(input_plot2, aes(x = AMBI, y = relab)) +
  geom_point() +
  stat_quantile(method = rq,
                formula = y ~ bs(x, df = 3),
                quantiles = 0.95) +
  facet_wrap(~ Family, scales = 'free')
ggsave("results/indicators_QRS/NOR_MB/all_families_MB_NOR_splines.pdf",
       p2, width = 15, height = 10)

# get splines data from plot
gb_NOR_MB <- ggplot_build(p2) #takes 10 minutes


gb1 <-
  gb_NOR_MB$data[[2]] %>% 
  data.frame() %>%
  group_by(group, PANEL) %>%
  nest() %>%
  mutate(
    peak = map(data, ~ .$x[which.max(.$y)]),
    find_peaks = map(
      data,
      ~ findpeaks(
        .$y,
        nups = 0,
        ndowns = 0 #,
        #minpeakheight = max(.$y) / 2.4
      )
    ),
    n_peaks = map_dbl(find_peaks, nrow)
  )

plot_id <-
  gb_NOR_MB$plot[[1]] %>%
  distinct(Family) %>%
  arrange(Family)

# get the number of peaks per ASVS and farm 
npeaks <-
  gb1 %>%
  select(n_peaks) %>%
  unnest(n_peaks) %>%
  bind_cols(plot_id)

abund_peaks <-
  npeaks %>%
  group_by(Family) %>%
  summarise_at(vars(n_peaks), list(sum), na.rm = T)

peak_dat <-
  gb1 %>%
  ungroup() %>%
  select(peak) %>%
  unnest(peak) %>%
  bind_cols(plot_id)
write.csv2(peak_dat, "results/indicators_QRS/NOR_MB/peak_dat_NOR_MB.csv")

# make qrs plots for all families

input_plot1 <- mb_tab %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>% 
  left_join(env_ambi)

for (i in 1:length(rownames(mb_tab))) {
  print(rownames(mb_tab)[i])
  input<- input_plot1 %>% select(Sample, rownames(mb_tab)[i], AMBI)
  colnames(input)[2] <- "relab"
  peak_dat2 <- peak_dat %>% filter(Family==rownames(mb_tab)[i])
  p1 <-
    ggplot(input, aes(x = AMBI, y = sqrt(relab))) +
    geom_point()+
    stat_quantile(method = rq,
                  formula = y ~ bs(x, df = 3),
                  quantiles = 0.95) +
    geom_vline(aes(xintercept = peak_dat2$peak)) +
    #facet_wrap(~ ASV, scales = 'free') +
    labs(y = 'Relative read abundance (sqrt)',
         x = 'AMBI' ,
         parse = T) +
    theme_bw(base_size = 12, base_family ="") %+replace%
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      strip.background = element_blank()
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0,6)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0,NA))
  p1 <- p1+ ggtitle(rownames(mb_tab)[i])
  p1
  ggsave(paste0("results/indicators_QRS/NOR_MB/plot_", rownames(mb_tab)[i], ".pdf"),
         p1, width = 6, height = 4.5)
}


# SCO MB

# MB Regression SCO
# read relab table filtered for SCO
mb_tab <- read.csv("data_created/reads_fam_MB_clean.csv", row.names = 1)
mb_tab <- mb_tab[,which(colnames(mb_tab) %in% sco_sams$Sample)]



input_plot2 <- mb_tab %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>% gather("Family", "relab", -Sample) %>%
  left_join(env_ambi)
p2 <-
  ggplot(input_plot2, aes(x = AMBI, y = relab)) +
  geom_point() +
  stat_quantile(method = rq,
                formula = y ~ bs(x, df = 3),
                quantiles = 0.95) +
  facet_wrap(~ Family, scales = 'free')
ggsave("results/indicators_QRS/SCO_MB/all_families_MB_SCO_splines.pdf",
       p2, width = 15, height = 10)

# get splines data from plot
gb_SCO_MB <- ggplot_build(p2) #takes 10 minutes


gb1 <-
  gb_SCO_MB$data[[2]] %>% 
  data.frame() %>%
  group_by(group, PANEL) %>%
  nest() %>%
  mutate(
    peak = map(data, ~ .$x[which.max(.$y)]),
    find_peaks = map(
      data,
      ~ findpeaks(
        .$y,
        nups = 0,
        ndowns = 0 #,
        #minpeakheight = max(.$y) / 2.4
      )
    ),
    n_peaks = map_dbl(find_peaks, nrow)
  )

plot_id <-
  gb_SCO_MB$plot[[1]] %>%
  distinct(Family) %>%
  arrange(Family)

# get the number of peaks per ASVS and farm 
npeaks <-
  gb1 %>%
  select(n_peaks) %>%
  unnest(n_peaks) %>%
  bind_cols(plot_id)

abund_peaks <-
  npeaks %>%
  group_by(Family) %>%
  summarise_at(vars(n_peaks), list(sum), na.rm = T)

peak_dat <-
  gb1 %>%
  ungroup() %>%
  select(peak) %>%
  unnest(peak) %>%
  bind_cols(plot_id)
write.csv2(peak_dat, "results/indicators_QRS/SCO_MB/peak_dat_SCO_MB.csv")

# make qrs plots for all families

input_plot1 <- mb_tab %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>% 
  left_join(env_ambi)

for (i in 1:length(rownames(mb_tab))) {
  print(rownames(mb_tab)[i])
  input<- input_plot1 %>% select(Sample, rownames(mb_tab)[i], AMBI)
  colnames(input)[2] <- "relab"
  peak_dat2 <- peak_dat %>% filter(Family==rownames(mb_tab)[i])
  p1 <-
    ggplot(input, aes(x = AMBI, y = sqrt(relab))) +
    geom_point()+
    stat_quantile(method = rq,
                  formula = y ~ bs(x, df = 3),
                  quantiles = 0.95) +
    geom_vline(aes(xintercept = peak_dat2$peak)) +
    #facet_wrap(~ ASV, scales = 'free') +
    labs(y = 'Relative read abundance (sqrt)',
         x = 'AMBI' ,
         parse = T) +
    theme_bw(base_size = 12, base_family ="") %+replace%
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      strip.background = element_blank()
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0,6)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0,NA))
  p1 <- p1+ ggtitle(rownames(mb_tab)[i])
  p1
  ggsave(paste0("results/indicators_QRS/SCO_MB/plot_", rownames(mb_tab)[i], ".pdf"),
         p1, width = 6, height = 4.5)
}


##MG

MG_tab <- read.csv("data_created/reads_fam_MG_clean.csv", row.names = 1)
input_plot2 <- MG_tab %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>% gather("Family", "relab", -Sample) %>%
  left_join(env_ambi)
p2 <-
  ggplot(input_plot2, aes(x = AMBI, y = relab)) +
  geom_point() +
  stat_quantile(method = rq,
                formula = y ~ bs(x, df = 3),
                quantiles = 0.95) +
  facet_wrap(~ Family, scales = 'free')
ggsave("results/indicators_QRS/MG/all_families_MG_splines.pdf",
       p2, width = 15, height = 10)

# get splines data from plot
gb_NOR_MG <- ggplot_build(p2) #takes 10 minutes


gb1 <-
  gb_NOR_MG$data[[2]] %>% 
  data.frame() %>%
  group_by(group, PANEL) %>%
  nest() %>%
  mutate(
    peak = map(data, ~ .$x[which.max(.$y)]),
    find_peaks = map(
      data,
      ~ findpeaks(
        .$y,
        nups = 0,
        ndowns = 0 #,
        #minpeakheight = max(.$y) / 2.4
      )
    ),
    n_peaks = map_dbl(find_peaks, nrow)
  )

plot_id <-
  gb_NOR_MG$plot[[1]] %>%
  distinct(Family) %>%
  arrange(Family)

# get the nuMGer of peaks per ASVS and farm 
npeaks <-
  gb1 %>%
  select(n_peaks) %>%
  unnest(n_peaks) %>%
  bind_cols(plot_id)

abund_peaks <-
  npeaks %>%
  group_by(Family) %>%
  summarise_at(vars(n_peaks), list(sum), na.rm = T)

peak_dat <-
  gb1 %>%
  ungroup() %>%
  select(peak) %>%
  unnest(peak) %>%
  bind_cols(plot_id)
write.csv2(peak_dat, "results/indicators_QRS/MG/peak_dat_MG.csv")

# make qrs plots for all families

input_plot1 <- MG_tab %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>% 
  left_join(env_ambi)

for (i in 1:length(rownames(MG_tab))) {
  print(rownames(MG_tab)[i])
  input<- input_plot1 %>% select(Sample, rownames(MG_tab)[i], AMBI)
  colnames(input)[2] <- "relab"
  peak_dat2 <- peak_dat %>% filter(Family==rownames(MG_tab)[i])
  p1 <-
    ggplot(input, aes(x = AMBI, y = sqrt(relab))) +
    geom_point()+
    stat_quantile(method = rq,
                  formula = y ~ bs(x, df = 3),
                  quantiles = 0.95) +
    geom_vline(aes(xintercept = peak_dat2$peak)) +
    #facet_wrap(~ ASV, scales = 'free') +
    labs(y = 'Relative read abundance (sqrt)',
         x = 'AMBI' ,
         parse = T) +
    theme_bw(base_size = 12, base_family ="") %+replace%
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      strip.background = element_blank()
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0,6)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0,NA))
  p1 <- p1+ ggtitle(rownames(MG_tab)[i])
  p1
  ggsave(paste0("results/indicators_QRS/MG/plot_", rownames(MG_tab)[i], ".pdf"),
         p1, width = 6, height = 4.5)
}

# NOR MG

# MG Regression NOR
# read relab table filtered for NOR
MG_tab <- read.csv("data_created/reads_fam_MG_clean.csv", row.names = 1)
MG_tab <- MG_tab[,which(colnames(MG_tab) %in% nor_sams$Sample)]

input_plot2 <- MG_tab %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>% gather("Family", "relab", -Sample) %>%
  left_join(env_ambi)
p2 <-
  ggplot(input_plot2, aes(x = AMBI, y = relab)) +
  geom_point() +
  stat_quantile(method = rq,
                formula = y ~ bs(x, df = 3),
                quantiles = 0.95) +
  facet_wrap(~ Family, scales = 'free')
ggsave("results/indicators_QRS/NOR_MG/all_families_MG_NOR_splines.pdf",
       p2, width = 15, height = 10)

# get splines data from plot
gb_NOR_MG <- ggplot_build(p2) #takes 10 minutes


gb1 <-
  gb_NOR_MG$data[[2]] %>% 
  data.frame() %>%
  group_by(group, PANEL) %>%
  nest() %>%
  mutate(
    peak = map(data, ~ .$x[which.max(.$y)]),
    find_peaks = map(
      data,
      ~ findpeaks(
        .$y,
        nups = 0,
        ndowns = 0 #,
        #minpeakheight = max(.$y) / 2.4
      )
    ),
    n_peaks = map_dbl(find_peaks, nrow)
  )

plot_id <-
  gb_NOR_MG$plot[[1]] %>%
  distinct(Family) %>%
  arrange(Family)

# get the nuMGer of peaks per ASVS and farm 
npeaks <-
  gb1 %>%
  select(n_peaks) %>%
  unnest(n_peaks) %>%
  bind_cols(plot_id)

abund_peaks <-
  npeaks %>%
  group_by(Family) %>%
  summarise_at(vars(n_peaks), list(sum), na.rm = T)

peak_dat <-
  gb1 %>%
  ungroup() %>%
  select(peak) %>%
  unnest(peak) %>%
  bind_cols(plot_id)
write.csv2(peak_dat, "results/indicators_QRS/NOR_MG/peak_dat_NOR_MG.csv")

# make qrs plots for all families

input_plot1 <- MG_tab %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>% 
  left_join(env_ambi)

for (i in 1:length(rownames(MG_tab))) {
  print(rownames(MG_tab)[i])
  input<- input_plot1 %>% select(Sample, rownames(MG_tab)[i], AMBI)
  colnames(input)[2] <- "relab"
  peak_dat2 <- peak_dat %>% filter(Family==rownames(MG_tab)[i])
  p1 <-
    ggplot(input, aes(x = AMBI, y = sqrt(relab))) +
    geom_point()+
    stat_quantile(method = rq,
                  formula = y ~ bs(x, df = 3),
                  quantiles = 0.95) +
    geom_vline(aes(xintercept = peak_dat2$peak)) +
    #facet_wrap(~ ASV, scales = 'free') +
    labs(y = 'Relative read abundance (sqrt)',
         x = 'AMBI' ,
         parse = T) +
    theme_bw(base_size = 12, base_family ="") %+replace%
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      strip.background = element_blank()
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0,6)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0,NA))
  p1 <- p1+ ggtitle(rownames(MG_tab)[i])
  p1
  ggsave(paste0("results/indicators_QRS/NOR_MG/plot_", rownames(MG_tab)[i], ".pdf"),
         p1, width = 6, height = 4.5)
}


# SCO MG

# MG Regression SCO
# read relab table filtered for SCO
MG_tab <- read.csv("data_created/reads_fam_MG_clean.csv", row.names = 1)
MG_tab <- MG_tab[,which(colnames(MG_tab) %in% sco_sams$Sample)]

input_plot2 <- MG_tab %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>% gather("Family", "relab", -Sample) %>%
  left_join(env_ambi)
p2 <-
  ggplot(input_plot2, aes(x = AMBI, y = relab)) +
  geom_point() +
  stat_quantile(method = rq,
                formula = y ~ bs(x, df = 3),
                quantiles = 0.95) +
  facet_wrap(~ Family, scales = 'free')
ggsave("results/indicators_QRS/SCO_MG/all_families_MG_SCO_splines.pdf",
       p2, width = 15, height = 10)

# get splines data from plot
gb_SCO_MG <- ggplot_build(p2) #takes 10 minutes


gb1 <-
  gb_SCO_MG$data[[2]] %>% 
  data.frame() %>%
  group_by(group, PANEL) %>%
  nest() %>%
  mutate(
    peak = map(data, ~ .$x[which.max(.$y)]),
    find_peaks = map(
      data,
      ~ findpeaks(
        .$y,
        nups = 0,
        ndowns = 0 #,
        #minpeakheight = max(.$y) / 2.4
      )
    ),
    n_peaks = map_dbl(find_peaks, nrow)
  )

plot_id <-
  gb_SCO_MG$plot[[1]] %>%
  distinct(Family) %>%
  arrange(Family)

# get the nuMGer of peaks per ASVS and farm 
npeaks <-
  gb1 %>%
  select(n_peaks) %>%
  unnest(n_peaks) %>%
  bind_cols(plot_id)

abund_peaks <-
  npeaks %>%
  group_by(Family) %>%
  summarise_at(vars(n_peaks), list(sum), na.rm = T)

peak_dat <-
  gb1 %>%
  ungroup() %>%
  select(peak) %>%
  unnest(peak) %>%
  bind_cols(plot_id)
write.csv2(peak_dat, "results/indicators_QRS/SCO_MG/peak_dat_SCO_MG.csv")

# make qrs plots for all families

input_plot1 <- MG_tab %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>% 
  left_join(env_ambi)

for (i in 1:length(rownames(MG_tab))) {
  print(rownames(MG_tab)[i])
  input<- input_plot1 %>% select(Sample, rownames(MG_tab)[i], AMBI)
  colnames(input)[2] <- "relab"
  peak_dat2 <- peak_dat %>% filter(Family==rownames(MG_tab)[i])
  p1 <-
    ggplot(input, aes(x = AMBI, y = sqrt(relab))) +
    geom_point()+
    stat_quantile(method = rq,
                  formula = y ~ bs(x, df = 3),
                  quantiles = 0.95) +
    geom_vline(aes(xintercept = peak_dat2$peak)) +
    #facet_wrap(~ ASV, scales = 'free') +
    labs(y = 'Relative read abundance (sqrt)',
         x = 'AMBI' ,
         parse = T) +
    theme_bw(base_size = 12, base_family ="") %+replace%
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      strip.background = element_blank()
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(0,6)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0,NA))
  p1 <- p1+ ggtitle(rownames(MG_tab)[i])
  p1
  ggsave(paste0("results/indicators_QRS/SCO_MG/plot_", rownames(MG_tab)[i], ".pdf"),
         p1, width = 6, height = 4.5)
}



## Plot potential indicators for all 4 options
## include only indicator candidates and peak data

# MB NOR
mb_tab <- read.csv("data_created/reads_fam_MB_clean.csv", row.names = 1)
mb_tab <- mb_tab[,which(colnames(mb_tab) %in% nor_sams$Sample)]
indic_NOR_MB <- read.csv2("results/indicators_QRS/NOR_MB/peak_indic_NOR_MB.csv", sep=";") %>% filter(indic!="n")
table(indic_NOR_MB$indic)
mb_tab_filt <- mb_tab[which(rownames(mb_tab) %in% indic_NOR_MB$Family),]
input_plot2 <- mb_tab_filt %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>% gather("Family", "relab", -Sample) %>%
  left_join(env_ambi)
p2 <-
  ggplot(input_plot2, aes(x = AMBI, y = relab)) +
  geom_point() +
  stat_quantile(method = rq,
                formula = y ~ bs(x, df = 3),
                quantiles = 0.95) +
  facet_wrap(~ Family, scales = 'free') +
  geom_vline(data=indic_NOR_MB, aes(xintercept = peak))
p2
#ggsave("results/indicators_QRS/NOR_MB/INDICATORS_MB_NOR_splines.pdf", p2, width = 15, height = 10)
ggsave("results/indicators_QRS/NOR_MB/INDICATORS_MB_NOR_splines.pdf", p2, width = 15/2, height = 10/2)

# MB SCO
mb_tab <- read.csv("data_created/reads_fam_MB_clean.csv", row.names = 1)
mb_tab <- mb_tab[,which(colnames(mb_tab) %in% sco_sams$Sample)]
indic_SCO_MB <- read.csv2("results/indicators_QRS/SCO_MB/peak_indic_SCO_MB.csv", sep=";") %>% filter(indic!="n")
table(indic_SCO_MB$indic)
mb_tab_filt <- mb_tab[which(rownames(mb_tab) %in% indic_SCO_MB$Family),]
input_plot2 <- mb_tab_filt %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>% gather("Family", "relab", -Sample) %>%
  left_join(env_ambi)
p2 <-
  ggplot(input_plot2, aes(x = AMBI, y = relab)) +
  geom_point() +
  stat_quantile(method = rq,
                formula = y ~ bs(x, df = 3),
                quantiles = 0.95) +
  facet_wrap(~ Family, scales = 'free') +
  geom_vline(data=indic_SCO_MB, aes(xintercept = peak))
p2
ggsave("results/indicators_QRS/SCO_MB/INDICATORS_MB_SCO_splines.pdf", p2, width = 15, height = 10)


# MG NOR
MG_tab <- read.csv("data_created/reads_fam_MG_clean.csv", row.names = 1)
MG_tab <- MG_tab[,which(colnames(MG_tab) %in% nor_sams$Sample)]
indic_NOR_MG <- read.csv2("results/indicators_QRS/NOR_MG/peak_indic_NOR_MG.csv", sep=";") %>% filter(indic!="n")
table(indic_NOR_MG$indic)
MG_tab_filt <- MG_tab[which(rownames(MG_tab) %in% indic_NOR_MG$Family),]
input_plot2 <- MG_tab_filt %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>% gather("Family", "relab", -Sample) %>%
  left_join(env_ambi)
p2 <-
  ggplot(input_plot2, aes(x = AMBI, y = relab)) +
  geom_point() +
  stat_quantile(method = rq,
                formula = y ~ bs(x, df = 3),
                quantiles = 0.95) +
  facet_wrap(~ Family, scales = 'free') +
  geom_vline(data=indic_NOR_MG, aes(xintercept = peak))
p2
ggsave("results/indicators_QRS/NOR_MG/INDICATORS_MG_NOR_splines.pdf", p2, width = 15, height = 10)

# MG SCO
MG_tab <- read.csv("data_created/reads_fam_MG_clean.csv", row.names = 1)
MG_tab <- MG_tab[,which(colnames(MG_tab) %in% sco_sams$Sample)]
indic_SCO_MG <- read.csv2("results/indicators_QRS/SCO_MG/peak_indic_SCO_MG.csv", sep=";") %>% filter(indic!="n")
table(indic_SCO_MG$indic)
MG_tab_filt <- MG_tab[which(rownames(MG_tab) %in% indic_SCO_MG$Family),]
input_plot2 <- MG_tab_filt %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>% gather("Family", "relab", -Sample) %>%
  left_join(env_ambi)
p2 <-
  ggplot(input_plot2, aes(x = AMBI, y = relab)) +
  geom_point() +
  stat_quantile(method = rq,
                formula = y ~ bs(x, df = 3),
                quantiles = 0.95) +
  facet_wrap(~ Family, scales = 'free') +
  geom_vline(data=indic_SCO_MG, aes(xintercept = peak))
p2
ggsave("results/indicators_QRS/SCO_MG/INDICATORS_MG_SCO_splines.pdf", p2, width = 15, height = 10)



### MB both cty
mb_tab <- read.csv("data_created/reads_fam_MB_clean.csv", row.names = 1)
indic_MB <- read.csv2("results/indicators_QRS/MB/peak_indic_MB.csv", sep=";") %>% filter(indic!="n")
table(indic_MB$indic)
mb_tab_filt <- mb_tab[which(rownames(mb_tab) %in% indic_MB$Family),]
input_plot2 <- mb_tab_filt %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>% gather("Family", "relab", -Sample) %>%
  left_join(env_ambi)
p2 <-
  ggplot(input_plot2, aes(x = AMBI, y = relab)) +
  geom_point() +
  stat_quantile(method = rq,
                formula = y ~ bs(x, df = 3),
                quantiles = 0.95) +
  facet_wrap(~ Family, scales = 'free') +
  geom_vline(data=indic_MB, aes(xintercept = peak))
p2
ggsave("results/indicators_QRS/MB/INDICATORS_MB_splines.pdf", p2, width = 15, height = 10)

### MG both cty
MG_tab <- read.csv("data_created/reads_fam_MG_clean.csv", row.names = 1)
indic_MG <- read.csv2("results/indicators_QRS/MG/peak_indic_MG.csv", sep=";") %>% filter(indic!="n")
table(indic_MG$indic)
MG_tab_filt <- MG_tab[which(rownames(MG_tab) %in% indic_MG$Family),]
input_plot2 <- MG_tab_filt %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>% gather("Family", "relab", -Sample) %>%
  left_join(env_ambi)
p2 <-
  ggplot(input_plot2, aes(x = AMBI, y = relab)) +
  geom_point() +
  stat_quantile(method = rq,
                formula = y ~ bs(x, df = 3),
                quantiles = 0.95) +
  facet_wrap(~ Family, scales = 'free') +
  geom_vline(data=indic_MG, aes(xintercept = peak))
p2
ggsave("results/indicators_QRS/MG/INDICATORS_MG_splines.pdf", p2, width = 15, height = 10)

library(tibble)
library(tidyr)
library(ggplot2)
library(dplyr)
library(vegan)
library(tidyverse)
library(reshape2)
library(janitor)
library(dendextend)
library(broom)


###################Calculate Diversity##############################
#read data
asv_tab_all<- read.csv("LIS_DUN_Minion_otu.csv", sep=",")%>% 
  column_to_rownames("OTU") #%>% select(-Taxon)
asv_tab_t <- asv_tab_all %>% t() %>% 
  as.data.frame() 

#Helper Function for alpha diversity
alpha_index<-function (x) {
  cbind.data.frame(
    shannon=diversity(x, index = "shannon")
  ) %>% rownames_to_column("sample")  # this line add a first column containing the rownames (i.e. the sample names)
}
#Calculate Shannon Index
alpha_div <-alpha_index(data)


###########Dendrogram###################
#Read data
asv_tab_all<- read.csv("ASV_table.csv", sep=",") %>% 
  column_to_rownames("asv") %>% select(-Taxon)

asv_tab_t <- asv_tab_all %>% t() %>% 
  as.data.frame() 
#Hellingertransformation
asv_rel<-asv_tab_t %>%  decostand(method="hellinger", MARGIN=1)
rowSums(asv_rel)
#Distancematrix
set.seed(666)
bray_cirpo<-vegdist(asv_rel,method="bray")
#Clusteranalysis
hc_ill<- hclust(bray_cirpo, method = "ward.D2")
#Plotdendrograms
plot(rotate(hc_ill, c(1, 2,3,5,6,4)),ylab = "Bray-Kurtis Dissimilarity",
     cex=1.25, cex.axis=1.25, cex.lab=1.25,hang = -1,
     las=1,xlab="", sub="", main = "Illumina-ASV")
symbols(c(1,4,5),rep(-0.45,3),circles=rep(1,3),add=TRUE,inches = .2,
        bg=c("#44AA99","#FC9272","#882255"),xpd=TRUE)
symbols(c(2,3,6),rep(-0.45,3),squares=rep(1,3),add=TRUE,inches = .4,
        bg=c("#FC9272","#44AA99","#882255"),xpd=TRUE)


###############Significant Differences/ Box-Violin-Plots####################
#Read data
asv_tab_all<- read.csv("Nanopore_Minion_genus_reads_05.csv", sep=",") %>% 
  column_to_rownames("genus")

asv_tab_t <- asv_tab_all %>% t() %>% 
  as.data.frame()

asv_tab_t <- asv_tab_all %>% t() %>% as.data.frame() %>% rownames_to_column("sample")
env <- read.csv("Stations.csv", sep=",")#Sample Code Data

asv_cipro <-as.data.frame(asv_tab_t) 

asv_4_kruskal_cipro <- inner_join(asv_cipro,env,by="sample") %>% column_to_rownames("sample") %>%
  select( -Farm,-Farm_SEQ,-Site)

#Significance testing
df_kruskal_cipro <- asv_4_kruskal_cipro %>% gather(key, value, -Sequencing) %>% 
  group_by(key) %>% 
  do(tidy(wilcox.test(x= .$value, g = .$Sequencing)))
df_kruskal_cipro
with(df_kruskal_cipro, key[p.value<0.05])

df_kruskal_significant_cipro<- df_kruskal_cipro %>% filter(p.value<0.05)

#Reshape data for plotting
asv_4_join_cipro <- asv_tab_all %>% rownames_to_column("key") 

df_4_boxplot_cipro<- inner_join(asv_4_join_cipro,df_kruskal_significant_cipro,by="key") %>% 
  select(-statistic,-p.value,-method,-alternative) %>% as.data.frame()
df_4_boxplot_cipro <- melt(df_4_boxplot_cipro, id.vars = "key")

env_df <- env %>% column_to_rownames("sample") %>% t() %>% as.data.frame() %>% 
  t() %>% as.data.frame() %>%rownames_to_column("variable")
df_4_boxplot_cipro<- inner_join(df_4_boxplot_cipro,env_df,by="variable")
df_4_boxplot_cipro<- mutate(df_4_boxplot_cipro, abundance= value*1)

df_4_boxplot_cipro_rep<-data.frame(lapply(df_4_boxplot_cipro, function(x) {gsub("Minion", "Nanopore", x)}))

#Help function for rounding digits
formatter <- function(...){
  function(x) format(round(x, 1), ...)
}

#Plot Box-violin graph
plot_significant_cipro <- ggplot(df_4_boxplot_cipro_rep, aes(x=key, y=as.numeric(abundance), fill=Sequencing)) + 
  geom_violin(position = position_dodge(0.9))+
  geom_boxplot(width=0.1,color="black",position = position_dodge(0.9))+
  theme_bw(base_size =17)+
  ylab("Relative read abundance in %")+
  labs(x=NULL)+
  labs(fill='Sequencing')+
  facet_wrap(~key,scales = "free")+
  scale_fill_manual(values=c("#DDCC77", "#88CCEE"))+
  theme(axis.text.x=element_blank(),axis.ticks=element_blank(),panel.grid.major   = element_blank())+
  theme(strip.background = element_rect(fill="white"))+
  theme(strip.text = element_text(colour = 'black'))+
  scale_y_continuous(labels = formatter(nsmall = 1))

plot_significant_cipro

###################Barcharts##############################
#Read data
phylum_tab <- read.csv("Genus_top_10_min_illu.csv", sep=",") 
phylum_tab <- read.csv("Top_10_minion_illumina_family.csv", sep=",") 

phylum_grouped_4_plot <- melt(phylum_tab, id.vars = "family")

#Prepare colorblind palette
safe_colorblind_palette_genus <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                                   "#44AA99", "#999933", "#882255",  "#6699CC","#74C476","#003C30",
                                   "#888888","#FC9272","#f0f0f0","#661100")

safe_colorblind_palette_family <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                                    "#44AA99", "#999933", "#882255", "#661100", "#6699CC","#f0f0f0", "#888888")
#Plot Barchart

plot_phylum_mean <- ggplot(phylum_grouped_4_plot, aes(variable, y=value, fill=family)) +  # changed: V7 -> phylum
  geom_bar(width=0.7, stat="identity")+
  scale_y_continuous(labels=function(x) x,expand = c(0.005,0.005)) +
  ylab("Relative read abundance in %") +
  xlab("Sequencing")+
  labs(fill="Family")+
  theme_light(base_size = 18)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), 
        panel.grid.major   = element_blank(),
        panel.spacing.x=unit(0.5, "lines") , panel.spacing.y=unit(5,"lines"),aspect.ratio = 1)+
  labs(fill='Family')+
  scale_fill_manual(values=safe_colorblind_palette_family)
facet_grid(.~factor(upwelling.intensity), 
           scales =  "free_x",space = "free_x")

plot_phylum_mean


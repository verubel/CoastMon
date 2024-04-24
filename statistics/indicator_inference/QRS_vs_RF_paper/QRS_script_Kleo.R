#PACKAGES
# This analysis depends on the following R packages (install if needed)
library(plyr)
library(dplyr)
library(tidyverse)
library(quantreg) 
library(splines) 
library(pracma) 
library(ggpmisc)
library(vegan)


#WORKING DIRECTORY
folderPath <- getwd()

#ANALYSIS

#load data

# List all the input files
fileList1 <- list.files(file.path(folderPath, "DATA"), pattern="asvtable", full.names=TRUE)
fileList2 <- list.files(file.path(folderPath, "DATA"), pattern="metadata", full.names=TRUE)
fileList3 <- list.files(file.path(folderPath, "DATA"), pattern="varimp", full.names=TRUE)
fileList4 <- list.files(file.path(folderPath, "DATA"), pattern="pred", full.names=TRUE)


#Get the top ASVs from each dataset using their proportional table and export in long format

for (country in 1:length(fileList1)) {
  
  # Extract the names of the datasets
  names <- sapply(strsplit(basename(fileList1[country]), "_"), `[`, 1)
  
  # Subset the asv table 
  asvtable <- read.csv(fileList1[country], header = TRUE, row.names = 1)
  
  # and get the proportional table
  proptable <- as.data.frame(round (prop.table (as.matrix(asvtable),2), digits=4)*100)
  
  # subset the metadata
  metadata <- read.csv(fileList2[country], header = TRUE)
  
  # create a table with the 250 most abundant ASVs and their relative abundance
  abund_ASV <-
    proptable %>%
    t(.) %>%
    as.data.frame() %>%
    gather(ASV, prop, 1:length(.)) %>%
    group_by(ASV) %>%
    summarise_at(vars(prop), list(mean), na.rm = T) %>%
    top_n(n = 250, wt = prop) 
  
  # create long format table
  abund_long <-
    proptable %>%
    t(.) %>%
    as.data.frame() %>%
    select (abund_ASV$ASV) %>% 
    rownames_to_column('Sample') %>%
    left_join(., metadata, by='Sample') %>%
    gather(ASV, abund, abund_ASV$ASV) 

# Quantile Regression Splines
  
 # extract the farm names 
 farm <- unique(abund_long$Farm)
  
  # create an empty data frame
  peak_all <- data.frame(rep(NA,250))
  
  # repeat for each farm
  for (i in farm) {
    
    # filter data for each farm
    bact_top <- filter(abund_long, Farm==i)
    
    # create plot with quantile regression and splines with 3 df
    p2 <-
      ggplot(bact_top, aes(x = IQI, y = abund)) +
      geom_point() +
      stat_quantile(method = rq,
                    formula = y ~ bs(x, df = 3),
                    quantiles = 0.95) +
      facet_wrap(~ ASV, scales = 'free')
    
    # get splines data from the plot
    gb <- ggplot_build(p2) 
    
    # generate peak data
    gb1 <-
      gb$data[[2]] %>% 
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
            ndowns = 0,
            minpeakheight = max(.$y) / 2.4
          )
        ),
        n_peaks = map_dbl(find_peaks, nrow)
      )
    
    plot_id <-
      gb$plot[[1]] %>%
      distinct(Farm, ASV) %>%
      arrange(ASV)
    
    # get the number of peaks per AVS and farm
    npeaks <-
      gb1 %>%
      select(n_peaks) %>%
      unnest() %>%
      bind_cols(plot_id) 
    
    # get the peak data
    peak_dat <-
      gb1 %>%
      ungroup() %>%
      select(peak) %>%
      unnest(peak) %>%
      bind_cols(plot_id) 
    
    # merge peaks by Farm and the number of peaks
    peak_table <- 
      full_join (npeaks, peak_dat, by = c("Farm", "ASV"))  %>%
      as.data.frame() %>%
      select(n_peaks:peak) 
    
    # substitute multiple peaks with NA for later calculations
    peak_table$peak[peak_table$n_peaks>2] <- NA
    
    # select only the peak column
    peak <- select(peak_table, peak)
    
    colnames(peak) <- paste (colnames(peak), i, sep="_")

    # fill in the table with peak information
    peak_all <- cbind(peak_all, peak)

  }

    peak_all  <- peak_all[,-c(1)]  
    rownames(peak_all) <- peak_table$ASV
    
    # calculations for ASV peak data
    calc <-data.frame(count = apply(peak_all, MARGIN = 1, function(x) sum(is.na(x))), # number of farms with no peak
                      mean=apply(peak_all, MARGIN = 1, function(x) mean(x, na.rm = T)),# mean peak value across farms
                      sd=apply(peak_all, MARGIN = 1, function(x) sd(x, na.rm = T))) # standard deviation across farms
    
    # combine calculations to peak table
    peak_all <- cbind(peak_all, calc)
    
     # add column with quality of asvs based on peak data
     asv_quality <-
      peak_all %>%
      mutate(
      Quality = as.numeric(as.character(cut(
        sd,
        breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 1),
        # difference where ASV reads peak in relation to IQI
        labels = c(0, 1, 2, 3, 4, 5)
      )))) #assign quality score (ALL 250)
  
    # pick indicators and add column with the Eco-Group
    indicators_pick <-
    asv_quality %>%
    filter(., count==0 | count==1 | count==2) %>% # select ASVs with peaks in at least five out of SEVEN farms
    filter(., Quality==0 | Quality==1 | Quality==2 | Quality==3) %>% # select ASVs with QS 0-3
    mutate(# add a new column that assigns an EG to the ASVs based on the mean peak value
      EG = cut(
        mean,
        # difference where ASV reads peak in relation to IQI
        breaks = c(0, 0.24, 0.43, 0.63, 0.74, 1),
        labels = c('V', 'IV', 'III', 'II', 'I')
      )) %>%
    rownames_to_column("ASV")
  
    # create long format table
    topqual_ASV <- 
      indicators_pick  %>%
      gather(Farm, peak, colnames(indicators_pick)[grepl("peak", colnames(indicators_pick))]) 
  
    topqual_ASV$Farm <- gsub("peak_*", "", topqual_ASV$Farm)
    
    # join with abundance table
    long_top_ASV <- 
      topqual_ASV %>%
      left_join(., abund_long,  by = c('ASV', 'Farm')) 
  
    # calculate eco groups proportions for each sample
    sample  <- unique (long_top_ASV$Sample) 
  
    z <- list()
    
  for (i in sample) {
    sample_abund <- filter (long_top_ASV, Sample==i) 
    SP <- ddply(sample_abund , .(EG), numcolwise(sum))  %>% 
      mutate(prop=abund/sum(abund)*100)  %>%
      select(EG, prop) %>% 
      spread(EG, prop)
    z <- c(z, list(SP))
  }
  
   names(z) <- sample
  
   # table with eco groups
   eg_prop_all <- do.call(rbind, z)
  
   # add molecular AMBI 
   bMBI <-
    eg_prop_all %>%
      mutate(bMBI_borja = (1.5 * II + 3 * III + 4.5 * IV + 6 * V) / 100)
  
# Molecular IQI calculation
    simp_even <- diversity(proptable, "simpson", MARGIN=2)/specnumber(proptable, MARGIN=2)
    molIQI <- ((0.38 * ((1-(bMBI$bMBI_borja/7))/(1-(min(bMBI$bMBI_borja)/7)))  + 
                0.08 * (simp_even / max(simp_even)) +  
                0.54 * (log10(specnumber(proptable, MARGIN=2))/max(log10(specnumber(proptable, MARGIN=2))))) - 0.4)                 /0.6

    molIQI <- as.data.frame(molIQI)  
  
    # add Molecular IQI column 
    indices <- cbind(bMBI, molIQI) %>% rownames_to_column('Sample')
       
    # join with metadata
    indices_meta <- left_join(indices, metadata, by="Sample")
  
    # add pass-fail categories
    indices_cat <-
      indices_meta %>%
          mutate(IQI_group = cut(
               IQI,
               breaks = c(0, 0.63, 1),
               labels = c('fail', 'pass'))) %>%
          mutate(molIQI_group = cut(
               molIQI,
               breaks = c(0, 0.63, 1),
               labels = c('fail', 'pass')))
  
 #Figure 1. Eco-Group distribution across samples in a) Norway (n=138 samples) and b) Scotland (n=92 samples)
  
  eco <- indicators_pick %>%
    group_by(EG) %>% 
    summarise(n = n()) %>%
    as.data.frame()
  eco <-  eco[c(5,4,3,2,1),]
  
  names <- c('Norway','Scotland')
  names2 <- c('a)', 'b)')
  
  b1 <- barplot (eco$n,
                 border="black",
                 xlab="Eco-groups",
                 ylab="Number of indicators",
                 main=paste(names2[country], names[country], sep=" "),
                 names.arg=c("I", "II", "III", "IV", "V"))
  
  
  #Figure 2. Top 20 ASVs assigned with the highest importance value by Random Forest (RF).Indicated with grey are the ASVs which were used as indicators with Quantile Regression Splines (QRS).
  
  varimp <- read.csv(fileList3[country], header=TRUE)
  varimp$inQRS <- ifelse(varimp$X %in% indicators_pick$ASV,"yes","no")
  varimp_sub <- varimp[1:20,]
  varimp_sub$cols <- ifelse(varimp_sub$X %in% indicators_pick$ASV,"grey","white")
  
  names <- c("Norway", "Scotland")
  names.main <- c("a)","b)")

  b1 <- barplot(varimp_sub$varimp, names=varimp_sub$X, las=2, 
                col=varimp_sub$cols,
                ylim=c(0,17), 
                ylab="ASV importance (RF)",
                main = paste(names.main[country], names[country], sep=" "))
  
  legend(x= 10, y=15, bty='n', 
         legend = c('Yes',
                    'No'), 
         fill=c("grey","white"),
         cex=1,
         border = 'black',
         title=expression(bold('Indicator ASV (QRS)')))

  
  #Figure 3 Linear regression plots showing the relationship between the Infaunal Quality Index #(IQI) and the molecular IQI as estimated by Quantile Regression Splines (QRS) and Random Forest #(RF) for a) Norway and b) Scotland salmon farms. The boxes indicate the two environmental #quality categories that IQI assigns the samples (i.e. red for (very good to good environmental #quality samples and green for moderate to poor environmental quality samples). The corresponding #R2 values are given for each regression plot.  
  
  pred <- read.csv(fileList4[country], header=TRUE)

  combined <- left_join(indices_cat, pred, by='Sample') 
  combined <-
    combined %>%
    column_to_rownames("Sample")
  
  names <- c('Norway','Scotland')
  names2 <- c('a)', 'b)')
  
  p1 <-
    ggplot(combined,
           aes(IQI, molIQI)) +
    geom_point(size = 1,alpha = 0) +
    geom_point(aes(shape=Farm))+
    scale_shape_manual(values=c(0, 1, 2, 3, 4, 5, 6)) + 
    geom_smooth(method = "lm", se=FALSE, colour="black")+
    stat_poly_eq(
      formula = y ~ x,
      aes(label =  paste(stat(eq.label),
                         stat(rr.label), sep = "*\", \"*")),
      parse = TRUE,
      size = 4
    ) +
    labs(y = 'Molecular IQI (QRS)',
         x = 'IQI' ,
         parse = T) +
    ylim(0,1) +
    xlim(0,1) +
    ggtitle(paste (names2[country], names[country], sep=" ")) +
    theme_bw(base_size = 15, base_family ="") %+replace%
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      strip.background = element_blank() 
    ) + 
    geom_rect(aes(xmin = 0, xmax = 0.64, ymin = 0, ymax = 0.64), fill=alpha("grey",0),
              colour= "red3") +
    geom_rect(aes(xmin = 0.64, xmax = 1, ymin = 0.64, ymax = 1), fill=alpha("grey",0),
              colour= "darkgreen")
  print(p1)
  
  p2 <-
    ggplot(combined,
           aes(IQI, pred)) +
    geom_point(size = 1,alpha = 0) +
    geom_point(aes(shape=Farm))+
    scale_shape_manual(values=c(0, 1, 2, 3, 4, 5, 6)) + 
    geom_smooth(method = "lm", se=FALSE, colour="black")+
    stat_poly_eq(
      formula = y ~ x,
      aes(label =  paste(stat(eq.label),
                         stat(rr.label), sep = "*\", \"*")),
      parse = TRUE,
      size = 4) +
    labs(y = 'Molecular IQI (RF)',
         x = 'IQI' ,
         parse = T) +
    ylim(0,1) +
    xlim(0,1) +
    ggtitle(" ") +
    theme_bw(base_size = 15, base_family ="") %+replace%
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      strip.background = element_blank() 
    ) + 
    geom_rect(aes(xmin = 0, xmax = 0.64, ymin = 0, ymax = 0.64), fill=alpha("grey",0),
              colour= "red3") +
    geom_rect(aes(xmin = 0.64, xmax = 1, ymin = 0.64, ymax = 1), fill=alpha("grey",0),
              colour= "darkgreen")
  print(p2)
  
  
  #Figure 4 Erroneously predicted samples by Quantile Regression Splines (QRS), Random Forest (RF) and both methods (RF+QRS) for a) Norway and b) Scotland salmon farms.
  
  combined$predQRS <- ifelse(combined$IQI_group==combined$molIQI_group,"Yes","No")
  combined$predRF <- ifelse(combined$IQI_group==combined$class_pred,"Yes","No")
  
  tot <- combined[combined[,"predQRS"] == 'No' | combined[,"predRF"]=="No",]
  common <- combined[combined[,"predQRS"] == 'No' & combined[,"predRF"]=="No",] 
  qrs <- combined[combined[,"predQRS"] == 'No' & combined[,"predRF"]=="Yes",]
  rf <- combined[combined[,"predQRS"] == 'Yes' & combined[,"predRF"]=="No",]
  
  #prepare transposed tables for barplots
  
  common_bar <- 
    common %>%
    select (IQI, molIQI, pred) %>%
    t(.) %>%
    as.data.frame
  colnames(common_bar) <- common$Sample
  
  qrs_bar <- 
    qrs %>%
    select (IQI, molIQI) %>%
    t(.) %>%
    as.data.frame
    colnames(qrs_bar) <- qrs$Sample
    
  rf_bar <- 
    rf %>%
    select (IQI, pred) %>%
    t(.) %>%
    as.data.frame
    colnames(rf_bar) <- rf$Sample
  
  
  #BArplots

  par(mar=c(5.1, 7.5, 3.1, 2.1), cex=1.2, las=2) # set the margins for plotting the data The default is c(5.1, 4.1, 4.1, 2.1)
  par(mfrow=c(2,2)) 
  
  #Barplot1
  cols <- c("white","steelblue","grey48")
  #space= c(c(0,0,0,1,0,0,1,0,0,1,0,0,1,0,0))
  
    common_bar <- as.matrix (common_bar)
    b1 <- barplot (common_bar, beside=TRUE, col=cols, horiz = T,
                   xlab="Biotic Index",
                   #ylab="Sample",
                   main="RF+QRS",
                   cex.axis = 1,
                   cex.lab=1.2)
   
   cols <- c("white","steelblue") 
   qrs_bar <- as.matrix (qrs_bar)
   b2 <- barplot (qrs_bar, beside=TRUE, col=cols, horiz = T,
                   xlab="Biotic Index",
                   #ylab="Sample",
                   main="QRS",
                   cex.axis = 1,
                   cex.lab=1.2)
   
  cols <- c("white","grey48") 
   rf_bar <- as.matrix (rf_bar)
   b3 <- barplot (rf_bar, beside=TRUE, col=cols, horiz = T,
                   xlab="Biotic Index",
                   #ylab="Sample",
                   main="RF",
                   cex.axis = 1,
                   cex.lab=1.2)
  
}      

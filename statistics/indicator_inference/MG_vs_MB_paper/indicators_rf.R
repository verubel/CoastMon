# indicators_rf
# Verena Rubel 
# RPTU Kaiserslautern Landau
# 30.10.2023 revised 10.01.2024

library(tidyr)
library(tibble)
library(ggplot2)
library(randomForest)
library(ggplot2)
library(dplyr)
library(caret) 
'%!in%' <- function(x,y)!('%in%'(x,y))

# load env with info containing grouping (good/bad)
out_sams <- read.csv("samples_out_NMDS_MB_MG.csv")
env0 <- read.csv("data/metadata_MB_MG.csv")
env <- env0[which(env0$Sample %!in% out_sams$S3), ]
env_cty <- env %>% dplyr::select(Sample, Country)
env_ambi <- env %>% dplyr::select(Sample, AMBI)
env_eqform <- env %>% mutate(EQ_class = ifelse(EQ %in% c(1, 2), "good", "bad")) %>% dplyr::select(Sample, EQ_class)
nor_sams<- env %>% filter(Country=="NOR") %>% droplevels()
sco_sams<- env %>% filter(Country=="SCO") %>% droplevels()

# expected: two results per MBMG each because we can do regression using AMBI and classification using EQ classes
# also: seperate analysis for NOR/SCO

# MB Regression NOR
# read relab table filtered for NOR
mb_tab <- read.csv("data_created/reads_fam_MB_clean.csv", row.names = 1)
mb_tab <- mb_tab[,which(colnames(mb_tab) %in% nor_sams$Sample)]

# add label for rf prediction
rf_input <- mb_tab %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>%
  left_join(env_ambi) %>% column_to_rownames("Sample")

### Random Forest pilot

set.seed(42)
model_acc <- randomForest(AMBI ~ ., data=rf_input, ntree=5000, importance = TRUE, proximity = TRUE)
print(model_acc)
importance(model_acc)
plot(model_acc)

# Random forest leave one out approach
### Random Forest LOO
train <- rf_input
uniqueIDs <- rownames(train)          
nruns <- length(uniqueIDs)    # number of cross validation runs: one for each unique ID           
crossclass <- match(rownames(train) , uniqueIDs)  
nobs <- nrow(na.omit(train))
crossPredict <- rep(NA, nobs)
choose_wdh <- 10 #number of repeats
result_rf_matrix <- matrix(nrow = choose_wdh, ncol = nruns)
colnames(result_rf_matrix) <- uniqueIDs
pb <- txtProgressBar(min = 0, max = choose_wdh, style = 3)

for (j in 1:choose_wdh) {
  crossPredict <- rep(NA, nobs)
  for (i in 1:nruns) {
    indtrain <- which(crossclass != i)
    indvalidate <- setdiff(1:nobs, indtrain)
    cat("Run", i, ": training only on observations with ID not", uniqueIDs[i], "\n")
    #IMPORTANT: set seed only for LOO without repetition!!!
    #set.seed(666)
    rf_df_CV <- randomForest(AMBI ~ ., data = train[indtrain,],
                             ntree = 5000, na.action = "na.omit", importance=TRUE)
    crossPredict[indvalidate] <- predict(rf_df_CV, train[indvalidate,])
    result_rf_matrix[j,i] <- crossPredict[indvalidate]
    saveRDS(rf_df_CV, 
            paste0("results/indicators_RF/MB_regression_NOR/models/model_without_", uniqueIDs[i], "_modelrepeat_", j, "X", ".rds"))
  }
  Sys.sleep(0.1)
  setTxtProgressBar(pb, j)
}

#end
close(pb)
print(t(result_rf_matrix))
write.csv(t(result_rf_matrix), "results/indicators_RF/MB_regression_NOR/predictions_MB_regressions_NOR.csv")
vgl_table <- t(result_rf_matrix) %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% gather("run", "value", -Sample) %>%
  group_by(Sample) %>%
  summarize(pred=mean(value), pred_sd=sd(value)) %>%
  left_join(env_ambi) %>%
  droplevels()
colnames(vgl_table) <- c("Sample", "pred", "pred_sd", "truth")

vgl_table_classes <- vgl_table %>% 
  mutate(class_truth=ifelse(truth<=1.2,1,ifelse(truth<=3.3,2,ifelse(truth<=4.3,3,ifelse(truth<=5.5,4,5))))) %>%
  mutate(class_pred=ifelse(pred<=1.2,1,ifelse(pred<=3.3,2,ifelse(pred<=4.3,3,ifelse(pred<=5.5,4,5))))) %>%
  as.data.frame()
vgl_table_classes[,5] <- as.factor(vgl_table_classes[,5])
vgl_table_classes[,6] <- as.factor(vgl_table_classes[,6])

levels_in_t_not_in_p <- setdiff(levels(vgl_table_classes$class_truth), levels(vgl_table_classes$class_pred))
levels(vgl_table_classes$class_pred) <- c(levels(vgl_table_classes$class_pred), levels_in_t_not_in_p)
conf_mat_final <- confusionMatrix(as.factor(vgl_table_classes$class_truth), as.factor(vgl_table_classes$class_pred))

a_tmp <- summary(lm(truth ~ pred,data=vgl_table_classes))
vgl_table_input <- vgl_table_classes %>% left_join(env)
plot_rf_nor <- ggplot(vgl_table_input, aes(x=truth, y=pred, shape=Farm))+
  annotate("rect", xmin=0, xmax=1.2, ymin=0, ymax=1.2, alpha=0, col="blue", fill="white") +
  annotate("rect", xmin=1.2, xmax=3.3, ymin=1.2, ymax=3.3, alpha=0, col="green", fill="white") +
  annotate("rect", xmin=3.3, xmax=4.3, ymin=3.3, ymax=4.3, alpha=0, col="yellow", fill="white") +
  annotate("rect", xmin=4.3, xmax=5.5, ymin=4.3, ymax=5.5, alpha=0, col="darkorange",fill="white") +
  annotate("rect", xmin=5.5, xmax=6.0, ymin=5.5, ymax=6.0, alpha=0, col="red", fill="white")+
  geom_smooth(method = "lm", se=FALSE, color="black", aes(group=1))+
  #geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed")+
  annotate("text",Inf,-Inf,hjust=1,vjust=-0.2,
           label=paste0("R² = ",round(a_tmp$adj.r.squared,digits=2), 
                        paste0("\nAccuracy = ", round(conf_mat_final$overall[1]*100, digits=2),
                               paste0("\nkappa = ",round(conf_mat_final$overall[2], digits=2)))))+
  theme_light()+
  theme(panel.grid = element_blank())+
  labs(x="AMBI", y="predicted AMBI")+
  scale_shape_manual(values=c(10, 15, 16, 17, 18, 9, 8))+
  geom_point()
plot_rf_nor
ggsave("results/indicators_RF/MB_regression_NOR/plot_MB_regression_NOR.pdf", plot_rf_nor, width = 4, height = 4)

# MB regression NOR GM

result_rf_matrix <- read.csv("results/indicators_RF/MB_regression_NOR/predictions_MB_regressions_NOR.csv", row.names = 1)
vgl_table <- result_rf_matrix %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% gather("run", "value", -Sample) %>%
  group_by(Sample) %>%
  summarize(pred=mean(value), pred_sd=sd(value)) %>%
  left_join(env_ambi) %>%
  droplevels()
colnames(vgl_table) <- c("Sample", "pred", "pred_sd", "truth")
vgl_table_classes <- vgl_table %>% 
  mutate(class_truth=ifelse(truth<=3.3, "good", "bad")) %>%
  mutate(class_pred=ifelse(pred<=3.3, "good", "bad")) %>%
  as.data.frame()

levels_in_t_not_in_p <- setdiff(levels(vgl_table_classes$class_truth), levels(vgl_table_classes$class_pred))
levels(vgl_table_classes$class_pred) <- c(levels(vgl_table_classes$class_pred), levels_in_t_not_in_p)
conf_mat_final <- confusionMatrix(as.factor(vgl_table_classes$class_truth), as.factor(vgl_table_classes$class_pred))

a_tmp <- summary(lm(truth ~ pred,data=vgl_table_classes))
vgl_table_input <- vgl_table_classes %>% left_join(env)
plot_rf_NOR <- ggplot(vgl_table_input, aes(x=truth, y=pred, shape=Farm))+
  annotate("rect", xmin=0, xmax=3.3, ymin=0, ymax=3.3, alpha=0, col="blue", fill="white") +
  annotate("rect", xmin=3.3, xmax=6.0, ymin=3.3, ymax=6.0, alpha=0, col="darkorange",fill="white") +
  geom_smooth(method = "lm", se=FALSE, color="black", aes(group=1))+
  #geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed")+
  annotate("text",Inf,-Inf,hjust=1,vjust=-0.2,
           label=paste0("R² = ",round(a_tmp$adj.r.squared,digits=2), 
                        paste0("\nAccuracy = ", round(conf_mat_final$overall[1]*100, digits=2),
                               paste0("\nkappa = ",round(conf_mat_final$overall[2], digits=2)))))+
  theme_light()+
  theme(panel.grid = element_blank())+
  labs(x="AMBI", y="predicted AMBI")+
  scale_shape_manual(values=c(10, 15, 16, 17, 18, 9, 8))+
  geom_point()
plot_rf_NOR
ggsave("results/indicators_RF/MB_regression_NOR/plot_MB_regression_NOR_GM.pdf", plot_rf_NOR, width = 4, height = 4)

# MB Regression SCO
# read relab table filtered for SCO
mb_tab <- read.csv("data_created/reads_fam_MB_clean.csv", row.names = 1)
mb_tab <- mb_tab[,which(colnames(mb_tab) %in% sco_sams$Sample)]
# add label for rf prediction
rf_input <- mb_tab %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>%
  left_join(env_ambi) %>% column_to_rownames("Sample")

### Random Forest LOO
train <- rf_input
uniqueIDs <- rownames(train)          
nruns <- length(uniqueIDs)    # number of cross validation runs: one for each unique ID           
crossclass <- match(rownames(train) , uniqueIDs)  
nobs <- nrow(na.omit(train))
crossPredict <- rep(NA, nobs)
choose_wdh <- 10 #number of repeats
result_rf_matrix <- matrix(nrow = choose_wdh, ncol = nruns)
colnames(result_rf_matrix) <- uniqueIDs
pb <- txtProgressBar(min = 0, max = choose_wdh, style = 3)

for (j in 1:choose_wdh) {
  crossPredict <- rep(NA, nobs)
  for (i in 1:nruns) {
    indtrain <- which(crossclass != i)
    indvalidate <- setdiff(1:nobs, indtrain)
    cat("Run", i, ": training only on observations with ID not", uniqueIDs[i], "\n")
    #IMPORTANT: set seed only for LOO without repetition!!!
    #set.seed(666)
    rf_df_CV <- randomForest(AMBI ~ ., data = train[indtrain,],
                             ntree = 5000, na.action = "na.omit", importance=TRUE)
    crossPredict[indvalidate] <- predict(rf_df_CV, train[indvalidate,])
    result_rf_matrix[j,i] <- crossPredict[indvalidate]
    saveRDS(rf_df_CV, 
            paste0("results/indicators_RF/MB_regression_SCO/models/model_without_", uniqueIDs[i], "_modelrepeat_", j, "X", ".rds"))
  }
  Sys.sleep(0.1)
  setTxtProgressBar(pb, j)
}

#end
close(pb)
print(t(result_rf_matrix))
write.csv(t(result_rf_matrix), "results/indicators_RF/MB_regression_SCO/predictions_MB_regressions_SCO.csv")
vgl_table <- t(result_rf_matrix) %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% gather("run", "value", -Sample) %>%
  group_by(Sample) %>%
  summarize(pred=mean(value), pred_sd=sd(value)) %>%
  left_join(env_ambi) %>%
  droplevels()
colnames(vgl_table) <- c("Sample", "pred", "pred_sd", "truth")

vgl_table_classes <- vgl_table %>% 
  mutate(class_truth=ifelse(truth<=1.2,1,ifelse(truth<=3.3,2,ifelse(truth<=4.3,3,ifelse(truth<=5.5,4,5))))) %>%
  mutate(class_pred=ifelse(pred<=1.2,1,ifelse(pred<=3.3,2,ifelse(pred<=4.3,3,ifelse(pred<=5.5,4,5))))) %>%
  as.data.frame()
vgl_table_classes[,5] <- as.factor(vgl_table_classes[,5])
vgl_table_classes[,6] <- as.factor(vgl_table_classes[,6])

levels_in_t_not_in_p <- setdiff(levels(vgl_table_classes$class_truth), levels(vgl_table_classes$class_pred))
levels(vgl_table_classes$class_pred) <- c(levels(vgl_table_classes$class_pred), levels_in_t_not_in_p)
conf_mat_final <- confusionMatrix(as.factor(vgl_table_classes$class_truth), as.factor(vgl_table_classes$class_pred))

a_tmp <- summary(lm(truth ~ pred,data=vgl_table_classes))
vgl_table_input <- vgl_table_classes %>% left_join(env)
plot_rf_SCO <- ggplot(vgl_table_input, aes(x=truth, y=pred, shape=Farm))+
  annotate("rect", xmin=0, xmax=1.2, ymin=0, ymax=1.2, alpha=0, col="blue", fill="white") +
  annotate("rect", xmin=1.2, xmax=3.3, ymin=1.2, ymax=3.3, alpha=0, col="green", fill="white") +
  annotate("rect", xmin=3.3, xmax=4.3, ymin=3.3, ymax=4.3, alpha=0, col="yellow", fill="white") +
  annotate("rect", xmin=4.3, xmax=5.5, ymin=4.3, ymax=5.5, alpha=0, col="darkorange",fill="white") +
  annotate("rect", xmin=5.5, xmax=6.0, ymin=5.5, ymax=6.0, alpha=0, col="red", fill="white")+
  geom_smooth(method = "lm", se=FALSE, color="black", aes(group=1))+
  #geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed")+
  annotate("text",Inf,-Inf,hjust=1,vjust=-0.2,
           label=paste0("R² = ",round(a_tmp$adj.r.squared,digits=2), 
                        paste0("\nAccuracy = ", round(conf_mat_final$overall[1]*100, digits=2),
                               paste0("\nkappa = ",round(conf_mat_final$overall[2], digits=2)))))+
  theme_light()+
  theme(panel.grid = element_blank())+
  labs(x="AMBI", y="predicted AMBI")+
  scale_shape_manual(values=c(10, 15, 16, 17, 18, 9, 8))+
  geom_point()
plot_rf_SCO
ggsave("results/indicators_RF/MB_regression_SCO/plot_MB_regression_SCO.pdf", plot_rf_SCO, width = 4, height = 4)



# MB regression SCO GM

result_rf_matrix <- read.csv("results/indicators_RF/MB_regression_SCO/predictions_MB_regressions_SCO.csv", row.names = 1)
vgl_table <- result_rf_matrix %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% gather("run", "value", -Sample) %>%
  group_by(Sample) %>%
  summarize(pred=mean(value), pred_sd=sd(value)) %>%
  left_join(env_ambi) %>%
  droplevels()
colnames(vgl_table) <- c("Sample", "pred", "pred_sd", "truth")
vgl_table_classes <- vgl_table %>% 
  mutate(class_truth=ifelse(truth<=3.3, "good", "bad")) %>%
  mutate(class_pred=ifelse(pred<=3.3, "good", "bad")) %>%
  as.data.frame()

levels_in_t_not_in_p <- setdiff(levels(vgl_table_classes$class_truth), levels(vgl_table_classes$class_pred))
levels(vgl_table_classes$class_pred) <- c(levels(vgl_table_classes$class_pred), levels_in_t_not_in_p)
conf_mat_final <- confusionMatrix(as.factor(vgl_table_classes$class_truth), as.factor(vgl_table_classes$class_pred))

a_tmp <- summary(lm(truth ~ pred,data=vgl_table_classes))
vgl_table_input <- vgl_table_classes %>% left_join(env)
plot_rf_SCO <- ggplot(vgl_table_input, aes(x=truth, y=pred, shape=Farm))+
  annotate("rect", xmin=0, xmax=3.3, ymin=0, ymax=3.3, alpha=0, col="blue", fill="white") +
  annotate("rect", xmin=3.3, xmax=6.0, ymin=3.3, ymax=6.0, alpha=0, col="darkorange",fill="white") +
  geom_smooth(method = "lm", se=FALSE, color="black", aes(group=1))+
  #geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed")+
  annotate("text",Inf,-Inf,hjust=1,vjust=-0.2,
           label=paste0("R² = ",round(a_tmp$adj.r.squared,digits=2), 
                        paste0("\nAccuracy = ", round(conf_mat_final$overall[1]*100, digits=2),
                               paste0("\nkappa = ",round(conf_mat_final$overall[2], digits=2)))))+
  theme_light()+
  theme(panel.grid = element_blank())+
  labs(x="AMBI", y="predicted AMBI")+
  scale_shape_manual(values=c(10, 15, 16, 17, 18, 9, 8))+
  geom_point()
plot_rf_SCO
ggsave("results/indicators_RF/MB_regression_SCO/plot_MB_regression_SCO_GM.pdf", plot_rf_SCO, width = 4, height = 4)



# Random forest leave one out approach
# expected: two results per MGMG each because we can do regression using AMBI and classification using EQ classes
# also: seperate analysis for NOR/SCO

# MG Regression NOR
# read relab table filtered for NOR
MG_tab <- read.csv("data_created/reads_fam_MG_clean.csv", row.names = 1)
MG_tab <- MG_tab[,which(colnames(MG_tab) %in% nor_sams$Sample)]

# add label for rf prediction
rf_input <- MG_tab %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>%
  left_join(env_ambi) %>% column_to_rownames("Sample")

### Random Forest LOO
train <- rf_input
uniqueIDs <- rownames(train)          
nruns <- length(uniqueIDs)    # nuMGer of cross validation runs: one for each unique ID           
crossclass <- match(rownames(train) , uniqueIDs)  
nobs <- nrow(na.omit(train))
crossPredict <- rep(NA, nobs)
choose_wdh <- 10 #nuMGer of repeats
result_rf_matrix <- matrix(nrow = choose_wdh, ncol = nruns)
colnames(result_rf_matrix) <- uniqueIDs
pb <- txtProgressBar(min = 0, max = choose_wdh, style = 3)

for (j in 1:choose_wdh) {
  crossPredict <- rep(NA, nobs)
  for (i in 1:nruns) {
    indtrain <- which(crossclass != i)
    indvalidate <- setdiff(1:nobs, indtrain)
    cat("Run", i, ": training only on observations with ID not", uniqueIDs[i], "\n")
    #IMPORTANT: set seed only for LOO without repetition!!!
    #set.seed(666)
    rf_df_CV <- randomForest(AMBI ~ ., data = train[indtrain,],
                             ntree = 5000, na.action = "na.omit", importance=TRUE)
    crossPredict[indvalidate] <- predict(rf_df_CV, train[indvalidate,])
    result_rf_matrix[j,i] <- crossPredict[indvalidate]
    saveRDS(rf_df_CV, 
            paste0("results/indicators_RF/MG_regression_NOR/models/model_without_", uniqueIDs[i], "_modelrepeat_", j, "X", ".rds"))
  }
  Sys.sleep(0.1)
  setTxtProgressBar(pb, j)
}

#end
close(pb)
print(t(result_rf_matrix))
write.csv(t(result_rf_matrix), "results/indicators_RF/MG_regression_NOR/predictions_MG_regressions_NOR.csv")
vgl_table <- t(result_rf_matrix) %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% gather("run", "value", -Sample) %>%
  group_by(Sample) %>%
  summarize(pred=mean(value), pred_sd=sd(value)) %>%
  left_join(env_ambi) %>%
  droplevels()
colnames(vgl_table) <- c("Sample", "pred", "pred_sd", "truth")

vgl_table_classes <- vgl_table %>% 
  mutate(class_truth=ifelse(truth<=1.2,1,ifelse(truth<=3.3,2,ifelse(truth<=4.3,3,ifelse(truth<=5.5,4,5))))) %>%
  mutate(class_pred=ifelse(pred<=1.2,1,ifelse(pred<=3.3,2,ifelse(pred<=4.3,3,ifelse(pred<=5.5,4,5))))) %>%
  as.data.frame()
vgl_table_classes[,5] <- as.factor(vgl_table_classes[,5])
vgl_table_classes[,6] <- as.factor(vgl_table_classes[,6])

levels_in_t_not_in_p <- setdiff(levels(vgl_table_classes$class_truth), levels(vgl_table_classes$class_pred))
levels(vgl_table_classes$class_pred) <- c(levels(vgl_table_classes$class_pred), levels_in_t_not_in_p)
conf_mat_final <- confusionMatrix(as.factor(vgl_table_classes$class_truth), as.factor(vgl_table_classes$class_pred))

a_tmp <- summary(lm(truth ~ pred,data=vgl_table_classes))
vgl_table_input <- vgl_table_classes %>% left_join(env)
plot_rf_nor <- ggplot(vgl_table_input, aes(x=truth, y=pred, shape=Farm))+
  annotate("rect", xmin=0, xmax=1.2, ymin=0, ymax=1.2, alpha=0, col="blue", fill="white") +
  annotate("rect", xmin=1.2, xmax=3.3, ymin=1.2, ymax=3.3, alpha=0, col="green", fill="white") +
  annotate("rect", xmin=3.3, xmax=4.3, ymin=3.3, ymax=4.3, alpha=0, col="yellow", fill="white") +
  annotate("rect", xmin=4.3, xmax=5.5, ymin=4.3, ymax=5.5, alpha=0, col="darkorange",fill="white") +
  annotate("rect", xmin=5.5, xmax=6.0, ymin=5.5, ymax=6.0, alpha=0, col="red", fill="white")+
  geom_smooth(method = "lm", se=FALSE, color="black", aes(group=1))+
  #geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed")+
  annotate("text",Inf,-Inf,hjust=1,vjust=-0.2,
           label=paste0("R² = ",round(a_tmp$adj.r.squared,digits=2), 
                        paste0("\nAccuracy = ", round(conf_mat_final$overall[1]*100, digits=2),
                               paste0("\nkappa = ",round(conf_mat_final$overall[2], digits=2)))))+
  theme_light()+
  theme(panel.grid = element_blank())+
  labs(x="AMBI", y="predicted AMBI")+
  scale_shape_manual(values=c(10, 15, 16, 17, 18, 9, 8))+
  geom_point()
plot_rf_nor
ggsave("results/indicators_RF/MG_regression_NOR/plot_MG_regression_NOR.pdf", plot_rf_nor, width = 4, height = 4)

# MG regression NOR GM

result_rf_matrix <- read.csv("results/indicators_RF/MG_regression_NOR/predictions_MG_regressions_NOR.csv", row.names = 1)
vgl_table <- result_rf_matrix %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% gather("run", "value", -Sample) %>%
  group_by(Sample) %>%
  summarize(pred=mean(value), pred_sd=sd(value)) %>%
  left_join(env_ambi) %>%
  droplevels()
colnames(vgl_table) <- c("Sample", "pred", "pred_sd", "truth")
vgl_table_classes <- vgl_table %>% 
  mutate(class_truth=ifelse(truth<=3.3, "good", "bad")) %>%
  mutate(class_pred=ifelse(pred<=3.3, "good", "bad")) %>%
  as.data.frame()

levels_in_t_not_in_p <- setdiff(levels(vgl_table_classes$class_truth), levels(vgl_table_classes$class_pred))
levels(vgl_table_classes$class_pred) <- c(levels(vgl_table_classes$class_pred), levels_in_t_not_in_p)
conf_mat_final <- confusionMatrix(as.factor(vgl_table_classes$class_truth), as.factor(vgl_table_classes$class_pred))

a_tmp <- summary(lm(truth ~ pred,data=vgl_table_classes))
vgl_table_input <- vgl_table_classes %>% left_join(env)
plot_rf_NOR <- ggplot(vgl_table_input, aes(x=truth, y=pred, shape=Farm))+
  annotate("rect", xmin=0, xmax=3.3, ymin=0, ymax=3.3, alpha=0, col="blue", fill="white") +
  annotate("rect", xmin=3.3, xmax=6.0, ymin=3.3, ymax=6.0, alpha=0, col="darkorange",fill="white") +
  geom_smooth(method = "lm", se=FALSE, color="black", aes(group=1))+
  #geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed")+
  annotate("text",Inf,-Inf,hjust=1,vjust=-0.2,
           label=paste0("R² = ",round(a_tmp$adj.r.squared,digits=2), 
                        paste0("\nAccuracy = ", round(conf_mat_final$overall[1]*100, digits=2),
                               paste0("\nkappa = ",round(conf_mat_final$overall[2], digits=2)))))+
  theme_light()+
  theme(panel.grid = element_blank())+
  labs(x="AMBI", y="predicted AMBI")+
  scale_shape_manual(values=c(10, 15, 16, 17, 18, 9, 8))+
  geom_point()
plot_rf_NOR
ggsave("results/indicators_RF/MG_regression_NOR/plot_MG_regression_NOR_GM.pdf", plot_rf_NOR, width = 4, height = 4)

# MG Regression SCO
# read relab table filtered for SCO
MG_tab <- read.csv("data_created/reads_fam_MG_clean.csv", row.names = 1)
MG_tab <- MG_tab[,which(colnames(MG_tab) %in% sco_sams$Sample)]
# add label for rf prediction
rf_input <- MG_tab %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>%
  left_join(env_ambi) %>% column_to_rownames("Sample")

### Random Forest LOO
train <- rf_input
uniqueIDs <- rownames(train)          
nruns <- length(uniqueIDs)    # number of cross validation runs: one for each unique ID           
crossclass <- match(rownames(train) , uniqueIDs)  
nobs <- nrow(na.omit(train))
crossPredict <- rep(NA, nobs)
choose_wdh <- 10 #nuMGer of repeats
result_rf_matrix <- matrix(nrow = choose_wdh, ncol = nruns)
colnames(result_rf_matrix) <- uniqueIDs
pb <- txtProgressBar(min = 0, max = choose_wdh, style = 3)

for (j in 1:choose_wdh) {
  crossPredict <- rep(NA, nobs)
  for (i in 1:nruns) {
    indtrain <- which(crossclass != i)
    indvalidate <- setdiff(1:nobs, indtrain)
    cat("Run", i, ": training only on observations with ID not", uniqueIDs[i], "\n")
    #IMPORTANT: set seed only for LOO without repetition!!!
    #set.seed(666)
    rf_df_CV <- randomForest(AMBI ~ ., data = train[indtrain,],
                             ntree = 5000, na.action = "na.omit", importance=TRUE)
    crossPredict[indvalidate] <- predict(rf_df_CV, train[indvalidate,])
    result_rf_matrix[j,i] <- crossPredict[indvalidate]
    saveRDS(rf_df_CV, 
            paste0("results/indicators_RF/MG_regression_SCO/models/model_without_", uniqueIDs[i], "_modelrepeat_", j, "X", ".rds"))
  }
  Sys.sleep(0.1)
  setTxtProgressBar(pb, j)
}

#end
close(pb)
print(t(result_rf_matrix))
write.csv(t(result_rf_matrix), "results/indicators_RF/MG_regression_SCO/predictions_MG_regressions_SCO.csv")
vgl_table <- t(result_rf_matrix) %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% gather("run", "value", -Sample) %>%
  group_by(Sample) %>%
  summarize(pred=mean(value), pred_sd=sd(value)) %>%
  left_join(env_ambi) %>%
  droplevels()
colnames(vgl_table) <- c("Sample", "pred", "pred_sd", "truth")

vgl_table_classes <- vgl_table %>% 
  mutate(class_truth=ifelse(truth<=1.2,1,ifelse(truth<=3.3,2,ifelse(truth<=4.3,3,ifelse(truth<=5.5,4,5))))) %>%
  mutate(class_pred=ifelse(pred<=1.2,1,ifelse(pred<=3.3,2,ifelse(pred<=4.3,3,ifelse(pred<=5.5,4,5))))) %>%
  as.data.frame()
vgl_table_classes[,5] <- as.factor(vgl_table_classes[,5])
vgl_table_classes[,6] <- as.factor(vgl_table_classes[,6])

levels_in_t_not_in_p <- setdiff(levels(vgl_table_classes$class_truth), levels(vgl_table_classes$class_pred))
levels(vgl_table_classes$class_pred) <- c(levels(vgl_table_classes$class_pred), levels_in_t_not_in_p)
conf_mat_final <- confusionMatrix(as.factor(vgl_table_classes$class_truth), as.factor(vgl_table_classes$class_pred))

a_tmp <- summary(lm(truth ~ pred,data=vgl_table_classes))
vgl_table_input <- vgl_table_classes %>% left_join(env)
plot_rf_SCO <- ggplot(vgl_table_input, aes(x=truth, y=pred, shape=Farm))+
  annotate("rect", xmin=0, xmax=1.2, ymin=0, ymax=1.2, alpha=0, col="blue", fill="white") +
  annotate("rect", xmin=1.2, xmax=3.3, ymin=1.2, ymax=3.3, alpha=0, col="green", fill="white") +
  annotate("rect", xmin=3.3, xmax=4.3, ymin=3.3, ymax=4.3, alpha=0, col="yellow", fill="white") +
  annotate("rect", xmin=4.3, xmax=5.5, ymin=4.3, ymax=5.5, alpha=0, col="darkorange",fill="white") +
  annotate("rect", xmin=5.5, xmax=6.0, ymin=5.5, ymax=6.0, alpha=0, col="red", fill="white")+
  geom_smooth(method = "lm", se=FALSE, color="black", aes(group=1))+
  #geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed")+
  annotate("text",Inf,-Inf,hjust=1,vjust=-0.2,
           label=paste0("R² = ",round(a_tmp$adj.r.squared,digits=2), 
                        paste0("\nAccuracy = ", round(conf_mat_final$overall[1]*100, digits=2),
                               paste0("\nkappa = ",round(conf_mat_final$overall[2], digits=2)))))+
  theme_light()+
  theme(panel.grid = element_blank())+
  labs(x="AMBI", y="predicted AMBI")+
  scale_shape_manual(values=c(10, 15, 16, 17, 18, 9, 8))+
  geom_point()
plot_rf_SCO
ggsave("results/indicators_RF/MG_regression_SCO/plot_MG_regression_SCO.pdf", plot_rf_SCO, width = 4, height = 4)



# MG regression SCO GM

result_rf_matrix <- read.csv("results/indicators_RF/MG_regression_SCO/predictions_MG_regressions_SCO.csv", row.names = 1)
vgl_table <- result_rf_matrix %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% gather("run", "value", -Sample) %>%
  group_by(Sample) %>%
  summarize(pred=mean(value), pred_sd=sd(value)) %>%
  left_join(env_ambi) %>%
  droplevels()
colnames(vgl_table) <- c("Sample", "pred", "pred_sd", "truth")
vgl_table_classes <- vgl_table %>% 
  mutate(class_truth=ifelse(truth<=3.3, "good", "bad")) %>%
  mutate(class_pred=ifelse(pred<=3.3, "good", "bad")) %>%
  as.data.frame()

levels_in_t_not_in_p <- setdiff(levels(vgl_table_classes$class_truth), levels(vgl_table_classes$class_pred))
levels(vgl_table_classes$class_pred) <- c(levels(vgl_table_classes$class_pred), levels_in_t_not_in_p)
conf_mat_final <- confusionMatrix(as.factor(vgl_table_classes$class_truth), as.factor(vgl_table_classes$class_pred))

a_tmp <- summary(lm(truth ~ pred,data=vgl_table_classes))
vgl_table_input <- vgl_table_classes %>% left_join(env)
plot_rf_SCO <- ggplot(vgl_table_input, aes(x=truth, y=pred, shape=Farm))+
  annotate("rect", xmin=0, xmax=3.3, ymin=0, ymax=3.3, alpha=0, col="blue", fill="white") +
  annotate("rect", xmin=3.3, xmax=6.0, ymin=3.3, ymax=6.0, alpha=0, col="darkorange",fill="white") +
  geom_smooth(method = "lm", se=FALSE, color="black", aes(group=1))+
  #geom_abline(intercept = 0, slope = 1, color="grey", linetype="dashed")+
  annotate("text",Inf,-Inf,hjust=1,vjust=-0.2,
           label=paste0("R² = ",round(a_tmp$adj.r.squared,digits=2), 
                        paste0("\nAccuracy = ", round(conf_mat_final$overall[1]*100, digits=2),
                               paste0("\nkappa = ",round(conf_mat_final$overall[2], digits=2)))))+
  theme_light()+
  theme(panel.grid = element_blank())+
  labs(x="AMBI", y="predicted AMBI")+
  scale_shape_manual(values=c(10, 15, 16, 17, 18, 9, 8))+
  geom_point()
plot_rf_SCO
ggsave("results/indicators_RF/MG_regression_SCO/plot_MG_regression_SCO_GM.pdf", plot_rf_SCO, width = 4, height = 4)




#
#
#
#
#
#
#
#
#
#







# RANDOM FOREST CLASSIFICATION


# MB classification NOR
# read relab table filtered for NOR
mb_tab <- read.csv("data_created/reads_fam_MB_clean.csv", row.names = 1)
mb_tab <- mb_tab[,which(colnames(mb_tab) %in% nor_sams$Sample)]
# add label for rf prediction
rf_input <- mb_tab %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>%
  left_join(env_eqform) %>% column_to_rownames("Sample")
rf_input$EQ_class = factor(rf_input$EQ_class) 
### Random Forest LOO
train <- rf_input
uniqueIDs <- rownames(train)          
nruns <- length(uniqueIDs)    # number of cross validation runs: one for each unique ID           
crossclass <- match(rownames(train) , uniqueIDs)  
nobs <- nrow(na.omit(train))
crossPredict <- rep(NA, nobs)
choose_wdh <- 10 #number of repeats
result_rf_matrix <- matrix(nrow = choose_wdh, ncol = nruns)
colnames(result_rf_matrix) <- uniqueIDs
pb <- txtProgressBar(min = 0, max = choose_wdh, style = 3)

for (j in 1:choose_wdh) {
  crossPredict <- rep(NA, nobs)
  for (i in 1:nruns) {
    indtrain <- which(crossclass != i)
    indvalidate <- setdiff(1:nobs, indtrain)
    cat("Run", i, ": training only on observations with ID not", uniqueIDs[i], "\n")
    #IMPORTANT: set seed only for LOO without repetition!!!
    #set.seed(666)
    rf_df_CV <- randomForest(EQ_class ~ ., data = train[indtrain,],
                             ntree = 5000, na.action = "na.omit", importance=TRUE)
    crossPredict[indvalidate] <- predict(rf_df_CV, train[indvalidate,])
    result_rf_matrix[j,i] <- crossPredict[indvalidate]
    saveRDS(rf_df_CV, 
            paste0("results/indicators_RF/MB_classification_NOR/models/model_without_", uniqueIDs[i], "_modelrepeat_", j, "X", ".rds"))
  }
  Sys.sleep(0.1)
  setTxtProgressBar(pb, j)
}

#end
close(pb)
result_rf_matrix

result_rf_matrix2 <- apply(result_rf_matrix, 2, function(x) ifelse(x == 1, "bad", ifelse(x == 2, "good", x)))
print(t(result_rf_matrix2))
write.csv(t(result_rf_matrix2), "results/indicators_RF/MB_classification_NOR/predictions_MB_classification_NOR.csv")
result_rf_matrix2<- read.csv("results/indicators_RF/MB_classification_NOR/predictions_MB_classification_NOR.csv", row.names = 1)
df <- result_rf_matrix2 %>% as.data.frame()
df2 <- df %>% mutate("majority"=apply(df,1,function(x) names(which.max(table(x)))))

vgl_table <- df2 %>% select(majority) %>%
  rownames_to_column("Sample") %>%
  left_join(env_eqform) %>%
  droplevels()
colnames(vgl_table) <- c("Sample", "pred", "truth")

conf_mat_final <- confusionMatrix(as.factor(vgl_table$truth), as.factor(vgl_table$pred))
print(conf_mat_final)
saveRDS(conf_mat_final, "results/indicators_RF/MB_classification_NOR/confusion_matrix_MB_classification_NOR.csv")




# MB classification SCO
# read relab table filtered for SCO
mb_tab <- read.csv("data_created/reads_fam_MB_clean.csv", row.names = 1)
mb_tab <- mb_tab[,which(colnames(mb_tab) %in% sco_sams$Sample)]
# add label for rf prediction
rf_input <- mb_tab %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>%
  left_join(env_eqform) %>% column_to_rownames("Sample")
rf_input$EQ_class = factor(rf_input$EQ_class) 
### Random Forest LOO
train <- rf_input
uniqueIDs <- rownames(train)          
nruns <- length(uniqueIDs)    # number of cross validation runs: one for each unique ID           
crossclass <- match(rownames(train) , uniqueIDs)  
nobs <- nrow(na.omit(train))
crossPredict <- rep(NA, nobs)
choose_wdh <- 10 #number of repeats
result_rf_matrix <- matrix(nrow = choose_wdh, ncol = nruns)
colnames(result_rf_matrix) <- uniqueIDs
pb <- txtProgressBar(min = 0, max = choose_wdh, style = 3)

for (j in 1:choose_wdh) {
  crossPredict <- rep(NA, nobs)
  for (i in 1:nruns) {
    indtrain <- which(crossclass != i)
    indvalidate <- setdiff(1:nobs, indtrain)
    cat("Run", i, ": training only on observations with ID not", uniqueIDs[i], "\n")
    #IMPORTANT: set seed only for LOO without repetition!!!
    #set.seed(666)
    rf_df_CV <- randomForest(EQ_class ~ ., data = train[indtrain,],
                             ntree = 5000, na.action = "na.omit", importance=TRUE)
    crossPredict[indvalidate] <- predict(rf_df_CV, train[indvalidate,])
    result_rf_matrix[j,i] <- crossPredict[indvalidate]
    saveRDS(rf_df_CV, 
            paste0("results/indicators_RF/MB_classification_SCO/models/model_without_", uniqueIDs[i], "_modelrepeat_", j, "X", ".rds"))
  }
  Sys.sleep(0.1)
  setTxtProgressBar(pb, j)
}

#end
close(pb)
result_rf_matrix

result_rf_matrix2 <- apply(result_rf_matrix, 2, function(x) ifelse(x == 1, "bad", ifelse(x == 2, "good", x)))
print(t(result_rf_matrix2))
write.csv(t(result_rf_matrix2), "results/indicators_RF/MB_classification_SCO/predictions_MB_classification_SCO.csv")

df <- t(result_rf_matrix2) %>% as.data.frame()
df2 <- df %>% mutate("majority"=apply(df,1,function(x) names(which.max(table(x)))))

vgl_table <- df2 %>% select(majority) %>%
  rownames_to_column("Sample") %>%
  left_join(env_eqform) %>%
  droplevels()
colnames(vgl_table) <- c("Sample", "pred", "truth")

conf_mat_final <- confusionMatrix(as.factor(vgl_table$truth), as.factor(vgl_table$pred))
print(conf_mat_final)
saveRDS(conf_mat_final, "results/indicators_RF/MB_classification_SCO/confusion_matrix_MB_classification_SCO.csv")








# RANDOM FOREST CLASSIFICATION


# MG classification NOR
# read relab table filtered for NOR
MG_tab <- read.csv("data_created/reads_fam_MG_clean.csv", row.names = 1)
MG_tab <- MG_tab[,which(colnames(MG_tab) %in% nor_sams$Sample)]
# add label for rf prediction
rf_input <- MG_tab %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>%
  left_join(env_eqform) %>% column_to_rownames("Sample")
rf_input$EQ_class = factor(rf_input$EQ_class) 
### Random Forest LOO
train <- rf_input
uniqueIDs <- rownames(train)          
nruns <- length(uniqueIDs)    # number of cross validation runs: one for each unique ID           
crossclass <- match(rownames(train) , uniqueIDs)  
nobs <- nrow(na.omit(train))
crossPredict <- rep(NA, nobs)
choose_wdh <- 10 #number of repeats
result_rf_matrix <- matrix(nrow = choose_wdh, ncol = nruns)
colnames(result_rf_matrix) <- uniqueIDs
pb <- txtProgressBar(min = 0, max = choose_wdh, style = 3)

for (j in 1:choose_wdh) {
  crossPredict <- rep(NA, nobs)
  for (i in 1:nruns) {
    indtrain <- which(crossclass != i)
    indvalidate <- setdiff(1:nobs, indtrain)
    cat("Run", i, ": training only on observations with ID not", uniqueIDs[i], "\n")
    #IMPORTANT: set seed only for LOO without repetition!!!
    #set.seed(666)
    rf_df_CV <- randomForest(EQ_class ~ ., data = train[indtrain,],
                             ntree = 5000, na.action = "na.omit", importance=TRUE)
    crossPredict[indvalidate] <- predict(rf_df_CV, train[indvalidate,])
    result_rf_matrix[j,i] <- crossPredict[indvalidate]
    saveRDS(rf_df_CV, 
            paste0("results/indicators_RF/MG_classification_NOR/models/model_without_", uniqueIDs[i], "_modelrepeat_", j, "X", ".rds"))
  }
  Sys.sleep(0.1)
  setTxtProgressBar(pb, j)
}

#end
close(pb)
result_rf_matrix

result_rf_matrix2 <- apply(result_rf_matrix, 2, function(x) ifelse(x == 1, "bad", ifelse(x == 2, "good", x)))
print(t(result_rf_matrix2))
write.csv(t(result_rf_matrix2), "results/indicators_RF/MG_classification_NOR/predictions_MG_classification_NOR.csv")
result_rf_matrix2<- read.csv("results/indicators_RF/MG_classification_NOR/predictions_MG_classification_NOR.csv", row.names = 1)
df <- result_rf_matrix2 %>% as.data.frame()
df2 <- df %>% mutate("majority"=apply(df,1,function(x) names(which.max(table(x)))))

vgl_table <- df2 %>% select(majority) %>%
  rownames_to_column("Sample") %>%
  left_join(env_eqform) %>%
  droplevels()
colnames(vgl_table) <- c("Sample", "pred", "truth")

conf_mat_final <- confusionMatrix(as.factor(vgl_table$truth), as.factor(vgl_table$pred))
print(conf_mat_final)
saveRDS(conf_mat_final, "results/indicators_RF/MG_classification_NOR/confusion_matrix_MG_classification_NOR.csv")




# MG classification SCO
# read relab table filtered for SCO
MG_tab <- read.csv("data_created/reads_fam_MG_clean.csv", row.names = 1)
MG_tab <- MG_tab[,which(colnames(MG_tab) %in% sco_sams$Sample)]
# add label for rf prediction
rf_input <- MG_tab %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>%
  left_join(env_eqform) %>% column_to_rownames("Sample")
rf_input$EQ_class = factor(rf_input$EQ_class) 
### Random Forest LOO
train <- rf_input
uniqueIDs <- rownames(train)          
nruns <- length(uniqueIDs)    # number of cross validation runs: one for each unique ID           
crossclass <- match(rownames(train) , uniqueIDs)  
nobs <- nrow(na.omit(train))
crossPredict <- rep(NA, nobs)
choose_wdh <- 10 #number of repeats
result_rf_matrix <- matrix(nrow = choose_wdh, ncol = nruns)
colnames(result_rf_matrix) <- uniqueIDs
pb <- txtProgressBar(min = 0, max = choose_wdh, style = 3)

for (j in 1:choose_wdh) {
  crossPredict <- rep(NA, nobs)
  for (i in 1:nruns) {
    indtrain <- which(crossclass != i)
    indvalidate <- setdiff(1:nobs, indtrain)
    cat("Run", i, ": training only on observations with ID not", uniqueIDs[i], "\n")
    #IMPORTANT: set seed only for LOO without repetition!!!
    #set.seed(666)
    rf_df_CV <- randomForest(EQ_class ~ ., data = train[indtrain,],
                             ntree = 5000, na.action = "na.omit", importance=TRUE)
    crossPredict[indvalidate] <- predict(rf_df_CV, train[indvalidate,])
    result_rf_matrix[j,i] <- crossPredict[indvalidate]
    saveRDS(rf_df_CV, 
            paste0("results/indicators_RF/MG_classification_SCO/models/model_without_", uniqueIDs[i], "_modelrepeat_", j, "X", ".rds"))
  }
  Sys.sleep(0.1)
  setTxtProgressBar(pb, j)
}

#end
close(pb)
result_rf_matrix

result_rf_matrix2 <- apply(result_rf_matrix, 2, function(x) ifelse(x == 1, "bad", ifelse(x == 2, "good", x)))
print(t(result_rf_matrix2))
write.csv(t(result_rf_matrix2), "results/indicators_RF/MG_classification_SCO/predictions_MG_classification_SCO.csv")

df <- t(result_rf_matrix2) %>% as.data.frame()
df2 <- df %>% mutate("majority"=apply(df,1,function(x) names(which.max(table(x)))))

vgl_table <- df2 %>% select(majority) %>%
  rownames_to_column("Sample") %>%
  left_join(env_eqform) %>%
  droplevels()
colnames(vgl_table) <- c("Sample", "pred", "truth")

conf_mat_final <- confusionMatrix(as.factor(vgl_table$truth), as.factor(vgl_table$pred))
print(conf_mat_final)
saveRDS(conf_mat_final, "results/indicators_RF/MG_classification_SCO/confusion_matrix_MG_classification_SCO.csv")




###
### extract indicators from regression models -> get varimp
###

# read all models
# get all varimps
# and combine

# MB NOR 
path <- "results/indicators_RF/MB_regression_NOR/models/"
temp = list.files(path, pattern="*rds")
temp2 <- paste0(path, temp)
sams <- sapply(strsplit(basename(temp), ".rds"), `[`, 1)
sams2 <- sapply(strsplit(basename(sams), "without_"), `[`, 2)
for (i in 1:length(temp2)) assign(sams2[i], readRDS(temp2[i]))

sams2[1]
head(importance(Ak_0_modelrepeat_1X))
data <- get(sams2[1])
head(importance(data))
model_0 <- importance(data) %>% as.data.frame() %>% rownames_to_column("ASV")

data2 <- get(sams2[2])
importance(get(sams2[2]))
model_1 <- importance(data2) %>% as.data.frame() %>% rownames_to_column("ASV")

model_comb <- model_0 %>% full_join(model_1, by="ASV")

sams3 <- paste0(sams2, "_VARIMP")
for (i in 1:length(sams2)) assign(sams3[i], importance(get(sams2[i])))

result_matrix <- matrix(ncol= length(sams3), nrow=55)
colnames(result_matrix) <- sams3

anfang <- get(sams3[1])
result_matrix[,1] <- anfang[,1]
rownames(result_matrix) <- rownames(anfang)

for (i in 1:length(sams3)) {
  vals_tmp <- get(sams3[i])
  result_matrix[,i] <- vals_tmp[,1]
}

write.csv(result_matrix, "results/indicators_RF/MB_regression_NOR/varimp_pro_model_MB_regression_NOR.csv")

imp_mean_df <- as.data.frame(rowMeans(result_matrix))
colnames(imp_mean_df) <- "varimp"
imp_mean_df <- imp_mean_df %>% arrange(-imp_mean_df)

write.csv(imp_mean_df, "results/indicators_RF/MB_regression_NOR/varimp_mean_allmodels_MB_regression_NOR.csv")

imp_mean_df_over0 <- imp_mean_df %>% filter(varimp>0) %>% droplevels() %>% rownames_to_column("Family")

plot_varimps2 <- ggplot(imp_mean_df_over0, aes(x=reorder(Family,varimp), y=varimp)) + 
  geom_bar(stat="identity") +
  theme_light() +
  theme(axis.text.x = element_text(angle=0), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+
  labs(y="Mean Variable Importance", x="Family")+
  coord_flip()
plot_varimps2
ggsave("results/indicators_RF/MB_regression_NOR/varimp_MB_regression_NOR.pdf", plot_varimps2, width=5, height =7)

in_new_plot3 <- as.data.frame(result_matrix) %>% rownames_to_column("family") %>% gather("a", "b", -family)
plot_varimps3 <- ggplot(in_new_plot3, aes(x=reorder(family,b), y=b)) + 
  geom_boxplot()+
  theme_light() +
  theme(axis.text.x = element_text(angle=0), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+
  labs(y="Mean Variable Importance", x="Family")+
  coord_flip()
plot_varimps3
ggsave("results/indicators_RF/MB_regression_NOR/varimp_MB_regression_NOR_BOXPLOT.pdf", plot_varimps3, width=5, height =7)

find_indicators <-  as.data.frame(result_matrix)
find_indicators$min <- do.call(pmin, find_indicators)
find_indicators_clean <- find_indicators %>% filter(min>0)
write.csv(find_indicators_clean, "results/indicators_RF/MB_regression_NOR/varimp_mean_allmodels_MB_regression_NOR_MINIMUM.csv")


#
# MB SCO 
path <- "results/indicators_RF/MB_regression_SCO/models/"
temp = list.files(path, pattern="*rds")
temp2 <- paste0(path, temp)
sams <- sapply(strsplit(basename(temp), ".rds"), `[`, 1)
sams2 <- sapply(strsplit(basename(sams), "without_"), `[`, 2)
for (i in 1:length(temp2)) assign(sams2[i], readRDS(temp2[i]))

sams2[1]
head(importance(Ak_0_modelrepeat_1X))
data <- get(sams2[1])
head(importance(data))
model_0 <- importance(data) %>% as.data.frame() %>% rownames_to_column("ASV")

data2 <- get(sams2[2])
importance(get(sams2[2]))
model_1 <- importance(data2) %>% as.data.frame() %>% rownames_to_column("ASV")

model_comb <- model_0 %>% full_join(model_1, by="ASV")

sams3 <- paste0(sams2, "_VARIMP")
for (i in 1:length(sams2)) assign(sams3[i], importance(get(sams2[i])))

result_matrix <- matrix(ncol= length(sams3), nrow=55)
colnames(result_matrix) <- sams3

anfang <- get(sams3[1])
result_matrix[,1] <- anfang[,1]
rownames(result_matrix) <- rownames(anfang)

for (i in 1:length(sams3)) {
  vals_tmp <- get(sams3[i])
  result_matrix[,i] <- vals_tmp[,1]
}

write.csv(result_matrix, "results/indicators_RF/MB_regression_SCO/varimp_pro_model_MB_regression_SCO.csv")

imp_mean_df <- as.data.frame(rowMeans(result_matrix))
colnames(imp_mean_df) <- "varimp"
imp_mean_df <- imp_mean_df %>% arrange(-imp_mean_df)

write.csv(imp_mean_df, "results/indicators_RF/MB_regression_SCO/varimp_mean_allmodels_MB_regression_SCO.csv")

imp_mean_df_over0 <- imp_mean_df %>% filter(varimp>0) %>% droplevels() %>% rownames_to_column("Family")

plot_varimps2 <- ggplot(imp_mean_df_over0, aes(x=reorder(Family,varimp), y=varimp)) + 
  geom_bar(stat="identity") +
  theme_light() +
  theme(axis.text.x = element_text(angle=0), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+
  labs(y="Mean Variable Importance", x="Family")+
  coord_flip()
plot_varimps2
ggsave("results/indicators_RF/MB_regression_SCO/varimp_MB_regression_SCO.pdf", plot_varimps2, width=5, height =7)

### NEW PLOT?

in_new_plot3 <- as.data.frame(result_matrix) %>% rownames_to_column("family") %>% gather("a", "b", -family)
plot_varimps3 <- ggplot(in_new_plot3, aes(x=reorder(family,b), y=b)) + 
  geom_boxplot()+
  theme_light() +
  theme(axis.text.x = element_text(angle=0), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+
  labs(y="Mean Variable Importance", x="Family")+
  coord_flip()
plot_varimps3
ggsave("results/indicators_RF/MB_regression_SCO/varimp_MB_regression_SCO_BOXPLOT.pdf", plot_varimps3, width=5, height =7)

find_indicators <-  as.data.frame(result_matrix)
find_indicators$min <- do.call(pmin, find_indicators)
find_indicators_clean <- find_indicators %>% filter(min>0)
dim(find_indicators_clean)
write.csv(find_indicators_clean, "results/indicators_RF/MB_regression_SCO/varimp_mean_allmodels_MB_regression_SCO_MINIMUM.csv")


###
### extract indicators from regression models -> get varimp
###

# read all models
# get all varimps
# and combine

# MG NOR 
path <- "results/indicators_RF/MG_regression_NOR/models/"
temp = list.files(path, pattern="*rds")
temp2 <- paste0(path, temp)
sams <- sapply(strsplit(basename(temp), ".rds"), `[`, 1)
sams2 <- sapply(strsplit(basename(sams), "without_"), `[`, 2)
for (i in 1:length(temp2)) assign(sams2[i], readRDS(temp2[i]))

sams2[1]
head(importance(Ak_0_modelrepeat_1X))
data <- get(sams2[1])
head(importance(data))
model_0 <- importance(data) %>% as.data.frame() %>% rownames_to_column("ASV")
asv_num <- dim(model_0)[1]
data2 <- get(sams2[2])
importance(get(sams2[2]))
model_1 <- importance(data2) %>% as.data.frame() %>% rownames_to_column("ASV")

model_comb <- model_0 %>% full_join(model_1, by="ASV")

sams3 <- paste0(sams2, "_VARIMP")
for (i in 1:length(sams2)) assign(sams3[i], importance(get(sams2[i])))

result_matrix <- matrix(ncol= length(sams3), nrow=asv_num)
colnames(result_matrix) <- sams3

anfang <- get(sams3[1])
result_matrix[,1] <- anfang[,1]
rownames(result_matrix) <- rownames(anfang)

for (i in 1:length(sams3)) {
  vals_tmp <- get(sams3[i])
  result_matrix[,i] <- vals_tmp[,1]
}

write.csv(result_matrix, "results/indicators_RF/MG_regression_NOR/varimp_pro_model_MG_regression_NOR.csv")

imp_mean_df <- as.data.frame(rowMeans(result_matrix))
colnames(imp_mean_df) <- "varimp"
imp_mean_df <- imp_mean_df %>% arrange(-imp_mean_df)

write.csv(imp_mean_df, "results/indicators_RF/MG_regression_NOR/varimp_mean_allmodels_MG_regression_NOR.csv")

imp_mean_df_over0 <- imp_mean_df %>% filter(varimp>0) %>% droplevels() %>% rownames_to_column("Family")

plot_varimps2 <- ggplot(imp_mean_df_over0, aes(x=reorder(Family,varimp), y=varimp)) + 
  geom_bar(stat="identity") +
  theme_light() +
  theme(axis.text.x = element_text(angle=0), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+
  labs(y="Mean Variable Importance", x="Family")+
  coord_flip()
plot_varimps2
ggsave("results/indicators_RF/MG_regression_NOR/varimp_MG_regression_NOR.pdf", plot_varimps2, width=5, height =7)

### NEW PLOT?

in_new_plot3 <- as.data.frame(result_matrix) %>% rownames_to_column("family") %>% gather("a", "b", -family)
plot_varimps3 <- ggplot(in_new_plot3, aes(x=reorder(family,b), y=b)) + 
  geom_boxplot()+
  theme_light() +
  theme(axis.text.x = element_text(angle=0), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+
  labs(y="Mean Variable Importance", x="Family")+
  coord_flip()
plot_varimps3
ggsave("results/indicators_RF/MG_regression_NOR/varimp_MG_regression_NOR_BOXPLOT.pdf", plot_varimps3, width=5, height =7)

find_indicators <-  as.data.frame(result_matrix) 
find_indicators$min <- do.call(pmin, find_indicators)
find_indicators_clean <- find_indicators %>% filter(min>0)
write.csv(find_indicators_clean, "results/indicators_RF/MG_regression_NOR/varimp_mean_allmodels_MG_regression_NOR_MINIMUM.csv")


#
# MB SCO 
path <- "results/indicators_RF/MG_regression_SCO/models/"
temp = list.files(path, pattern="*rds")
temp2 <- paste0(path, temp)
sams <- sapply(strsplit(basename(temp), ".rds"), `[`, 1)
sams2 <- sapply(strsplit(basename(sams), "without_"), `[`, 2)
for (i in 1:length(temp2)) assign(sams2[i], readRDS(temp2[i]))

sams2[1]
head(importance(Ak_0_modelrepeat_1X))
data <- get(sams2[1])
head(importance(data))
dim(importance(data))
model_0 <- importance(data) %>% as.data.frame() %>% rownames_to_column("ASV")

data2 <- get(sams2[2])
importance(get(sams2[2]))
model_1 <- importance(data2) %>% as.data.frame() %>% rownames_to_column("ASV")

model_comb <- model_0 %>% full_join(model_1, by="ASV")

sams3 <- paste0(sams2, "_VARIMP")
for (i in 1:length(sams2)) assign(sams3[i], importance(get(sams2[i])))

result_matrix <- matrix(ncol= length(sams3), nrow=71)
colnames(result_matrix) <- sams3

anfang <- get(sams3[1])
result_matrix[,1] <- anfang[,1]
rownames(result_matrix) <- rownames(anfang)

for (i in 1:length(sams3)) {
  vals_tmp <- get(sams3[i])
  result_matrix[,i] <- vals_tmp[,1]
}

write.csv(result_matrix, "results/indicators_RF/MG_regression_SCO/varimp_pro_model_MG_regression_SCO.csv")

imp_mean_df <- as.data.frame(rowMeans(result_matrix))
colnames(imp_mean_df) <- "varimp"
imp_mean_df <- imp_mean_df %>% arrange(-imp_mean_df)

write.csv(imp_mean_df, "results/indicators_RF/MG_regression_SCO/varimp_mean_allmodels_MG_regression_SCO.csv")

imp_mean_df_over0 <- imp_mean_df %>% filter(varimp>0) %>% droplevels() %>% rownames_to_column("Family")

plot_varimps2 <- ggplot(imp_mean_df_over0, aes(x=reorder(Family,varimp), y=varimp)) + 
  geom_bar(stat="identity") +
  theme_light() +
  theme(axis.text.x = element_text(angle=0), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+
  labs(y="Mean Variable Importance", x="Family")+
  coord_flip()
plot_varimps2
ggsave("results/indicators_RF/MG_regression_SCO/varimp_MG_regression_SCO.pdf", plot_varimps2, width=5, height =7)

in_new_plot3 <- as.data.frame(result_matrix) %>% rownames_to_column("family") %>% gather("a", "b", -family)
plot_varimps3 <- ggplot(in_new_plot3, aes(x=reorder(family,b), y=b)) + 
  geom_boxplot()+
  theme_light() +
  theme(axis.text.x = element_text(angle=0), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+
  labs(y="Mean Variable Importance", x="Family")+
  coord_flip()
plot_varimps3
ggsave("results/indicators_RF/MG_regression_SCO/varimp_MG_regression_SCO_BOXPLOT.pdf", plot_varimps3, width=5, height =7)

find_indicators <-  as.data.frame(result_matrix)
find_indicators$min <- do.call(pmin, find_indicators)
find_indicators_clean <- find_indicators %>% filter(min>0)
dim(find_indicators_clean)
write.csv(find_indicators_clean, "results/indicators_RF/MG_regression_SCO/varimp_mean_allmodels_MG_regression_SCO_MINIMUM.csv")


#
#
#
## variable importance classification RF
#
#
#


# MB NOR 
path <- "results/indicators_RF/MB_classification_NOR/models/"
temp = list.files(path, pattern="*rds")
temp2 <- paste0(path, temp)
sams <- sapply(strsplit(basename(temp), ".rds"), `[`, 1)
sams2 <- sapply(strsplit(basename(sams), "without_"), `[`, 2)
for (i in 1:length(temp2)) assign(sams2[i], readRDS(temp2[i]))

sams2[1]
head(importance(Ak_0_modelrepeat_1X))
data <- get(sams2[1])
head(importance(data))
model_0 <- importance(data) %>% as.data.frame() %>% rownames_to_column("ASV")

data2 <- get(sams2[2])
importance(get(sams2[2]))
dim(importance(get(sams2[2])))
model_1 <- importance(data2) %>% as.data.frame() %>% rownames_to_column("ASV")

model_comb <- model_0 %>% full_join(model_1, by="ASV")

sams3 <- paste0(sams2, "_VARIMP")
for (i in 1:length(sams2)) assign(sams3[i], importance(get(sams2[i])))

result_matrix <- matrix(ncol= length(sams3), nrow=55)
colnames(result_matrix) <- sams3

anfang <- get(sams3[1])
result_matrix[,1] <- anfang[,3]
rownames(result_matrix) <- rownames(anfang)

for (i in 1:length(sams3)) {
  vals_tmp <- get(sams3[i])
  result_matrix[,i] <- vals_tmp[,3]
}

write.csv(result_matrix, "results/indicators_RF/MB_classification_NOR/varimp_pro_model_MB_classification_NOR.csv")

imp_mean_df <- as.data.frame(rowMeans(result_matrix))
colnames(imp_mean_df) <- "varimp"
imp_mean_df <- imp_mean_df %>% arrange(-imp_mean_df)

write.csv(imp_mean_df, "results/indicators_RF/MB_classification_NOR/varimp_mean_allmodels_MB_classification_NOR.csv")

imp_mean_df_over0 <- imp_mean_df %>% filter(varimp>0) %>% droplevels() %>% rownames_to_column("Family")

plot_varimps2 <- ggplot(imp_mean_df_over0, aes(x=reorder(Family,varimp), y=varimp)) + 
  geom_bar(stat="identity") +
  theme_light() +
  theme(axis.text.x = element_text(angle=0), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+
  labs(y="Mean Variable Importance", x="Family")+
  coord_flip()
plot_varimps2
ggsave("results/indicators_RF/MB_classification_NOR/varimp_MB_classification_NOR.pdf", plot_varimps2, width=5, height =7)

in_new_plot3 <- as.data.frame(result_matrix) %>% rownames_to_column("family") %>% gather("a", "b", -family)
plot_varimps3 <- ggplot(in_new_plot3, aes(x=reorder(family,b), y=b)) + 
  geom_boxplot()+
  theme_light() +
  theme(axis.text.x = element_text(angle=0), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+
  labs(y="Mean Variable Importance", x="Family")+
  coord_flip()
plot_varimps3
ggsave("results/indicators_RF/MB_classification_NOR/varimp_MB_classification_NOR_BOXPLOT.pdf", plot_varimps3, width=5, height =7)

find_indicators <-  as.data.frame(result_matrix)
find_indicators$min <- do.call(pmin, find_indicators)
find_indicators_clean <- find_indicators %>% filter(min>0)
dim(find_indicators_clean)
write.csv(find_indicators_clean, "results/indicators_RF/MB_classification_NOR/varimp_mean_allmodels_MB_classification_NOR_MINIMUM.csv")


## MB SCO

# MB SCO 
path <- "results/indicators_RF/MB_classification_SCO/models/"
temp = list.files(path, pattern="*rds")
temp2 <- paste0(path, temp)
sams <- sapply(strsplit(basename(temp), ".rds"), `[`, 1)
sams2 <- sapply(strsplit(basename(sams), "without_"), `[`, 2)
for (i in 1:length(temp2)) assign(sams2[i], readRDS(temp2[i]))

sams2[1]
head(importance(Ak_0_modelrepeat_1X))
data <- get(sams2[1])
head(importance(data))
model_0 <- importance(data) %>% as.data.frame() %>% rownames_to_column("ASV")

data2 <- get(sams2[2])
importance(get(sams2[2]))
dim(importance(get(sams2[2])))
model_1 <- importance(data2) %>% as.data.frame() %>% rownames_to_column("ASV")

model_comb <- model_0 %>% full_join(model_1, by="ASV")

sams3 <- paste0(sams2, "_VARIMP")
for (i in 1:length(sams2)) assign(sams3[i], importance(get(sams2[i])))

result_matrix <- matrix(ncol= length(sams3), nrow=55)
colnames(result_matrix) <- sams3

anfang <- get(sams3[1])
result_matrix[,1] <- anfang[,3]
rownames(result_matrix) <- rownames(anfang)

for (i in 1:length(sams3)) {
  vals_tmp <- get(sams3[i])
  result_matrix[,i] <- vals_tmp[,3]
}

write.csv(result_matrix, "results/indicators_RF/MB_classification_SCO/varimp_pro_model_MB_classification_SCO.csv")

imp_mean_df <- as.data.frame(rowMeans(result_matrix))
colnames(imp_mean_df) <- "varimp"
imp_mean_df <- imp_mean_df %>% arrange(-imp_mean_df)

write.csv(imp_mean_df, "results/indicators_RF/MB_classification_SCO/varimp_mean_allmodels_MB_classification_SCO.csv")

imp_mean_df_over0 <- imp_mean_df %>% filter(varimp>0) %>% droplevels() %>% rownames_to_column("Family")

plot_varimps2 <- ggplot(imp_mean_df_over0, aes(x=reorder(Family,varimp), y=varimp)) + 
  geom_bar(stat="identity") +
  theme_light() +
  theme(axis.text.x = element_text(angle=0), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+
  labs(y="Mean Variable Importance", x="Family")+
  coord_flip()
plot_varimps2
ggsave("results/indicators_RF/MB_classification_SCO/varimp_MB_classification_SCO.pdf", plot_varimps2, width=5, height =7)

in_new_plot3 <- as.data.frame(result_matrix) %>% rownames_to_column("family") %>% gather("a", "b", -family)
plot_varimps3 <- ggplot(in_new_plot3, aes(x=reorder(family,b), y=b)) + 
  geom_boxplot()+
  theme_light() +
  theme(axis.text.x = element_text(angle=0), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+
  labs(y="Mean Variable Importance", x="Family")+
  coord_flip()
plot_varimps3
ggsave("results/indicators_RF/MB_classification_SCO/varimp_MB_classification_SCO_BOXPLOT.pdf", plot_varimps3, width=5, height =7)

find_indicators <-  as.data.frame(result_matrix)
find_indicators$min <- do.call(pmin, find_indicators)
find_indicators_clean <- find_indicators %>% filter(min>0)
dim(find_indicators_clean)
write.csv(find_indicators_clean, "results/indicators_RF/MB_classification_SCO/varimp_mean_allmodels_MB_classification_SCO_MINIMUM.csv")



## MG

# MG NOR 
path <- "results/indicators_RF/MG_classification_NOR/models/"
temp = list.files(path, pattern="*rds")
temp2 <- paste0(path, temp)
sams <- sapply(strsplit(basename(temp), ".rds"), `[`, 1)
sams2 <- sapply(strsplit(basename(sams), "without_"), `[`, 2)
for (i in 1:length(temp2)) assign(sams2[i], readRDS(temp2[i]))

sams2[1]
head(importance(Ak_0_modelrepeat_1X))
data <- get(sams2[1])
head(importance(data))
model_0 <- importance(data) %>% as.data.frame() %>% rownames_to_column("ASV")

data2 <- get(sams2[2])
importance(get(sams2[2]))
dim(importance(get(sams2[2])))
model_1 <- importance(data2) %>% as.data.frame() %>% rownames_to_column("ASV")

model_coMG <- model_0 %>% full_join(model_1, by="ASV")

sams3 <- paste0(sams2, "_VARIMP")
for (i in 1:length(sams2)) assign(sams3[i], importance(get(sams2[i])))

result_matrix <- matrix(ncol= length(sams3), nrow=71)
colnames(result_matrix) <- sams3

anfang <- get(sams3[1])
result_matrix[,1] <- anfang[,3]
rownames(result_matrix) <- rownames(anfang)

for (i in 1:length(sams3)) {
  vals_tmp <- get(sams3[i])
  result_matrix[,i] <- vals_tmp[,3]
}

write.csv(result_matrix, "results/indicators_RF/MG_classification_NOR/varimp_pro_model_MG_classification_NOR.csv")

imp_mean_df <- as.data.frame(rowMeans(result_matrix))
colnames(imp_mean_df) <- "varimp"
imp_mean_df <- imp_mean_df %>% arrange(-imp_mean_df)

write.csv(imp_mean_df, "results/indicators_RF/MG_classification_NOR/varimp_mean_allmodels_MG_classification_NOR.csv")

imp_mean_df_over0 <- imp_mean_df %>% filter(varimp>0) %>% droplevels() %>% rownames_to_column("Family")

plot_varimps2 <- ggplot(imp_mean_df_over0, aes(x=reorder(Family,varimp), y=varimp)) + 
  geom_bar(stat="identity") +
  theme_light() +
  theme(axis.text.x = element_text(angle=0), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+
  labs(y="Mean Variable Importance", x="Family")+
  coord_flip()
plot_varimps2
ggsave("results/indicators_RF/MG_classification_NOR/varimp_MG_classification_NOR.pdf", plot_varimps2, width=5, height =7)

in_new_plot3 <- as.data.frame(result_matrix) %>% rownames_to_column("family") %>% gather("a", "b", -family)
plot_varimps3 <- ggplot(in_new_plot3, aes(x=reorder(family,b), y=b)) + 
  geom_boxplot()+
  theme_light() +
  theme(axis.text.x = element_text(angle=0), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+
  labs(y="Mean Variable Importance", x="Family")+
  coord_flip()
plot_varimps3
ggsave("results/indicators_RF/MG_classification_NOR/varimp_MG_classification_NOR_BOXPLOT.pdf", plot_varimps3, width=5, height =7)

find_indicators <-  as.data.frame(result_matrix)
find_indicators$min <- do.call(pmin, find_indicators)
find_indicators_clean <- find_indicators %>% filter(min>0)
dim(find_indicators_clean)
write.csv(find_indicators_clean, "results/indicators_RF/MG_classification_NOR/varimp_mean_allmodels_MG_classification_NOR_MINIMUM.csv")


## MG SCO

# MG SCO 
path <- "results/indicators_RF/MG_classification_SCO/models/"
temp = list.files(path, pattern="*rds")
temp2 <- paste0(path, temp)
sams <- sapply(strsplit(basename(temp), ".rds"), `[`, 1)
sams2 <- sapply(strsplit(basename(sams), "without_"), `[`, 2)
for (i in 1:length(temp2)) assign(sams2[i], readRDS(temp2[i]))

sams2[1]
head(importance(Ak_0_modelrepeat_1X))
data <- get(sams2[1])
head(importance(data))
model_0 <- importance(data) %>% as.data.frame() %>% rownames_to_column("ASV")

data2 <- get(sams2[2])
importance(get(sams2[2]))
dim(importance(get(sams2[2])))
model_1 <- importance(data2) %>% as.data.frame() %>% rownames_to_column("ASV")

model_coMG <- model_0 %>% full_join(model_1, by="ASV")

sams3 <- paste0(sams2, "_VARIMP")
for (i in 1:length(sams2)) assign(sams3[i], importance(get(sams2[i])))

result_matrix <- matrix(ncol= length(sams3), nrow=71)
colnames(result_matrix) <- sams3

anfang <- get(sams3[1])
result_matrix[,1] <- anfang[,3]
rownames(result_matrix) <- rownames(anfang)

for (i in 1:length(sams3)) {
  vals_tmp <- get(sams3[i])
  result_matrix[,i] <- vals_tmp[,3]
}

write.csv(result_matrix, "results/indicators_RF/MG_classification_SCO/varimp_pro_model_MG_classification_SCO.csv")

imp_mean_df <- as.data.frame(rowMeans(result_matrix))
colnames(imp_mean_df) <- "varimp"
imp_mean_df <- imp_mean_df %>% arrange(-imp_mean_df)

write.csv(imp_mean_df, "results/indicators_RF/MG_classification_SCO/varimp_mean_allmodels_MG_classification_SCO.csv")

imp_mean_df_over0 <- imp_mean_df %>% filter(varimp>0) %>% droplevels() %>% rownames_to_column("Family")

plot_varimps2 <- ggplot(imp_mean_df_over0, aes(x=reorder(Family,varimp), y=varimp)) + 
  geom_bar(stat="identity") +
  theme_light() +
  theme(axis.text.x = element_text(angle=0), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+
  labs(y="Mean Variable Importance", x="Family")+
  coord_flip()
plot_varimps2
ggsave("results/indicators_RF/MG_classification_SCO/varimp_MG_classification_SCO.pdf", plot_varimps2, width=5, height =7)

in_new_plot3 <- as.data.frame(result_matrix) %>% rownames_to_column("family") %>% gather("a", "b", -family)
plot_varimps3 <- ggplot(in_new_plot3, aes(x=reorder(family,b), y=b)) + 
  geom_boxplot()+
  theme_light() +
  theme(axis.text.x = element_text(angle=0), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())+
  labs(y="Mean Variable Importance", x="Family")+
  coord_flip()
plot_varimps3
ggsave("results/indicators_RF/MG_classification_SCO/varimp_MG_classification_SCO_BOXPLOT.pdf", plot_varimps3, width=5, height =7)

find_indicators <-  as.data.frame(result_matrix)
find_indicators$min <- do.call(pmin, find_indicators)
find_indicators_clean <- find_indicators %>% filter(min>0)
dim(find_indicators_clean)
write.csv(find_indicators_clean, "results/indicators_RF/MG_classification_SCO/varimp_mean_allmodels_MG_classification_SCO_MINIMUM.csv")




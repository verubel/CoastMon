# IQI prediction based on prop-table

# load required libraries 
library(tibble)
library(randomForest)
library(caret)
library(tidyr)
library(dplyr)


# get rarefied/transformed table
nor_tab_rare <- read.csv("proptable_nor250.csv", row.names = 1) %>%
  column_to_rownames("ASV")

# get metadata table
env <- read.csv2("metadata.csv") %>% filter(Country=="Norway")
env_iqi <- env %>% select(Sample, IQI)
env_farm <- env %>% select(Sample, Farm)

# if no 4th root tansformation only transform initial table
rf_input_df <- t(nor_tab_rare) %>% as.data.frame() %>% 
  rownames_to_column("Sample") %>% left_join(env_iqi) %>%
  column_to_rownames("Sample")

train <- rf_input_df
uniqueIDs <- rownames(train)          
nruns <- length(uniqueIDs)    # number of crossvalidation runs = one for each unique ID for leave one out approach          
crossclass <- match(rownames(train) , uniqueIDs)  
nobs <- nrow(na.omit(train))
crossPredict <- rep(NA, nobs)
choose_wdh <- 10 # number of reps

# construct empty matrix to save the results
result_rf_matrix <- matrix(nrow = choose_wdh, ncol = nruns)
colnames(result_rf_matrix) <- uniqueIDs

# include progress bar
pb <- txtProgressBar(min = 0, max = choose_wdh, style = 3)

for (j in 1:choose_wdh) {
  crossPredict <- rep(NA, nobs)
  for (i in 1:nruns) {
    indtrain <- which(crossclass != i)
    indvalidate <- setdiff(1:nobs, indtrain)
    cat("Run", i, ": training only on observations with ID not", uniqueIDs[i], "\n")
    #IMPORTANT: do not set seed here because of 10-fold repetiton!
    #set.seed(666)
    rf_df_CV <- randomForest(IQI ~ ., data = train[indtrain,],
                             ntree = 500, na.action = "na.omit", importance=TRUE)
    crossPredict[indvalidate] <- predict(rf_df_CV, train[indvalidate,])
    result_rf_matrix[j,i] <- crossPredict[indvalidate]
    saveRDS(rf_df_CV, 
            paste0("models/model_NOR_relab_without_", uniqueIDs[i], "_modelrepeat_", j, ".rds"))
  }
  Sys.sleep(0.1)
  setTxtProgressBar(pb, j)
}

#end
close(pb)
print(t(result_rf_matrix))
write.csv(t(result_rf_matrix), "result_pred_nor.csv")
# 

# read matrix of results
result_tab <- read.csv("result_pred_nor.csv", row.names = 1)
vgl <- result_tab %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample") %>% gather("run", "value", -Sample) %>%
  group_by(Sample) %>%
  summarize(pred=mean(value), pred_sd=sd(value))
colnames(vgl) <- c("Sample", "pred", "truth")

# join metadata to results and rename columns
vgl_table <- vgl %>%
  left_join(env_iqi) %>%
  droplevels()
colnames(vgl_table) <- c("Sample", "pred", "pred_sd", "truth")

# categorize the predictions
vgl_table_classes <- vgl_table %>% 
  mutate(class_truth=ifelse(truth>=0.64,"pass", "fail")) %>%
  mutate(class_pred=ifelse(pred>=0.64,"pass", "fail")) %>%
  left_join(env)

# write table of predictions and create confusion matrix
write.csv(vgl_table_classes[,1:6], "table_of_predictions_relab_norway.csv")

CM_tmp <- confusionMatrix(as.factor(vgl_table_classes$class_truth),
                          as.factor(vgl_table_classes$class_pred))
						  
						  
						  
						  

# GET VARAMPLE IMPORTANCE MEASUES:
# read all models
# get all varimps
# and combine


# load all models
path <- "/models/"
temp = list.files(path, pattern="*rds")
temp2 <- paste0("/path/to/models/", temp)
sams <- sapply(strsplit(basename(temp), ".rds"), `[`, 1)
sams2 <- sapply(strsplit(basename(sams), "without_"), `[`, 2)
for (i in 1:length(temp2)) assign(sams2[i], readRDS(temp2[i]))

# get importance measure of first model
sams2[1]
head(importance(AK_0_G1_RepA_modelrepeat_1))
data <- get(sams2[1])
head(importance(data))
model_0 <- importance(data) %>% as.data.frame() %>% rownames_to_column("ASV")

# start construction of data frame containing variable importances
data2 <- get(sams2[2])
importance(get(sams2[2]))
model_1 <- importance(data2) %>% as.data.frame() %>% rownames_to_column("ASV")

model_comb <- model_0 %>% full_join(model_1, by="ASV")

sams3 <- paste0(sams2, "_VARIMP")
for (i in 1:length(sams2)) assign(sams3[i], importance(get(sams2[i])))

result_matrix <- matrix(ncol= length(sams3), nrow=250)
colnames(result_matrix) <- sams3

anfang <- get(sams3[1])
result_matrix[,1] <- anfang[,1]
rownames(result_matrix) <- rownames(anfang)

# get variable importances for each ASV and write in output file
for (i in 1:length(sams3)) {
  vals_tmp <- get(sams3[i])
  result_matrix[,i] <- vals_tmp[,1]
}

write.csv(result_matrix, "varimp_pro_model_nor.csv")

# get mean variable importances for each ASV and write in output file in descending order
imp_mean_df <- as.data.frame(rowMeans(result_matrix))
colnames(imp_mean_df) <- "varimp"
imp_mean_df <- imp_mean_df %>% arrange(-imp_mean_df)
write.csv(imp_mean_df, "varimp_mean_nor.csv")

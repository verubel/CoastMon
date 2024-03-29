# random_forest_regression_LOO.R
# Verena Rubel 
# RPTU Kaiserslautern Landau
# created 30.10.2023 revised 10.01.2024
# regression
# Leave one out LOO
# AMBI prediction

# load required modules
library(tidyr)
library(tibble)
library(randomForest)
library(dplyr)

# load mapping file with info containing grouping for prediction
env <- read.csv("data/mapping.csv")
env_ambi <- env %>% dplyr::select(Sample, AMBI)

# random forest leave one out approach
# input read table
mb_tab <- read.csv("data_created/reads.csv", row.names = 1)

# add label for rf prediction
rf_input <- mb_tab %>% t() %>% as.data.frame() %>% rownames_to_column("Sample") %>%
  left_join(env_ambi) %>% column_to_rownames("Sample")

# random forest LOO with 10 repeats 
train <- rf_input
uniqueIDs <- rownames(train)          
nruns <- length(uniqueIDs)   
crossclass <- match(rownames(train) , uniqueIDs)  
nobs <- nrow(na.omit(train))
crossPredict <- rep(NA, nobs)
choose_wdh <- 10 
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
            paste0("results/models/model_without_", uniqueIDs[i], "_modelrepeat_", j, "X", ".rds"))
  }
  Sys.sleep(0.1)
  setTxtProgressBar(pb, j)
}

#end
close(pb)
print(t(result_rf_matrix))
write.csv(t(result_rf_matrix), "results/predictions_regressions.csv")

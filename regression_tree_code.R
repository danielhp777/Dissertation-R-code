library(rpart)
library(rpart.plot)
library(dplyr)
library(ggplot2)

data <- read.table("Data2.txt",header=TRUE)

# separating data (6 measurements per subject for most subjects, m=481)
time <- c(0,2,4,8,16,24)
calwk <- data$calwk
logrna <- as.numeric(data$logrna)
nnrti <- factor(data$nnrti)
txday <- data$txday
cens <- factor(data$cens)
trtarm <- factor(data$trtarm)
patid <- data$patid

# correctly defined data frame
hiv <- data.frame(patid, nnrti, calwk, txday, logrna, cens, trtarm)

#Split the data into training and testing data, 80:20 split, ie 385:96 subjects
set.seed(158)
patients <- seq(1:481)
test_ids <- sample(patients, size = 96, replace = FALSE) # choosing the patient ID's for test data
test_ids <- sort(test_ids)
train_ids <- setdiff(patients, test_ids)

hiv_train <- hiv %>% filter(patid %in% train_ids)
hiv_test <- hiv %>% filter(patid %in% test_ids)

# remove patient id
hiv_train <- hiv_train[,-1]
hiv_test <- hiv_test[,-1]

# fit the regression tree to training data
model1 <- rpart(logrna ~ ., data = hiv_train, method = "anova") # seems to be the best

rpart.plot(model1, type = 3, digits = 3)

# check how tree performs on the training set
predictions_train <- predict(model1, hiv_train)

# plot predictions vs actual responses for training data
D1 <- data.frame(real_values = hiv_train$logrna, Predictions = predictions_train)

ggplot(D1, aes(x = Predictions, y = real_values)) + 
  geom_point(alpha = 0.6, color = "red") +
  geom_smooth(method = "lm", se = FALSE, color = "black") + # best fit line
  #geom_abline(slope = 1, intercept = 0, color = 'green', linewidth = 1) +
  labs(
    x = "Predicted viral load", y = "Actual viral load",
    title = 
    "Predicted vs Actual viral loads for the training data (regression tree)") +
  theme_minimal()

# mse for training data
mse_train <- mean((predictions_train - hiv_train$logrna)^2)

# predict viral loads for testing data
predictions_test <- predict(model1, hiv_test)

# plot predictions vs actual responses for testing data
D2 <- data.frame(real_values = hiv_test$logrna, Predictions = predictions_test)

ggplot(D2, aes(x = Predictions, y = real_values)) + 
  geom_point(alpha = 0.6, aes(color = 'Data points')) +
  geom_smooth(method = "lm", se = FALSE, aes(color = 'Regression line')) + # best fit line
  geom_line(linewidth = 1, data = data.frame(x = range(D2$Predictions), y = range(D2$Predictions)), aes(x=x, y=y, color = 'Identity line')) +
  scale_color_manual(name = 'Legend', values = c('Data points'= 'purple', 'Regression line' = 'black', 'Identity line' = 'green')) +         
  labs(
    x = "Predicted viral load", y = "Actual viral load",
    title = 
      "Predicted vs Actual viral loads for the testing data (regression tree)") +
  theme_minimal()

# mse for testing data
mse_test <- mean((predictions_test - hiv_test$logrna)^2)

# mae for testing data
mae_test <- mean(abs(predictions_test - hiv_test$logrna))

# R-squared for actual vs predicted responses (located in summary of linear model)
reg_tree_r_squared <- lm(hiv_test$logrna ~ predictions_test)
summary(reg_tree_r_squared)

# correlation between actual and predicted responses
cor(hiv_test$logrna, predictions_test)


# XGBoost
library(xgboost)
library(Matrix)
# sort data (XGBoost requires data as DMatrices)
hiv_train_x <- hiv_train[,-4]
hiv_train_y <- hiv_train$logrna
hiv_test_x <- hiv_test[,-4]
hiv_test_y <- hiv_test$logrna

# xgboost requires numeric entries, so we convert factors to integers
hiv_train_x[] <- lapply(hiv_train_x, function(col){ if (is.factor(col)) as.integer(col) else col})
hiv_test_x[] <- lapply(hiv_test_x, function(col){ if (is.factor(col)) as.integer(col) else col})

# ensure training and testing data are matrices, not data frames
hiv_train_x <- as.matrix(hiv_train_x)
hiv_test_x <- as.matrix(hiv_test_x)

xgb_training <- xgb.DMatrix(data = hiv_train_x, label = hiv_train_y)
xgb_testing <- xgb.DMatrix(data = hiv_test_x, label = hiv_test_y)

# tune eta, nrounds and max_depth parameters using cross-validation
eta <- c(0.1, 0.3, 0.5)
max_depth <- c(3,4,5)
best_parameters <- list() # tracks which eta and max_depth works best
best_nrounds <- 0
best_rmse <- 10000 # used to track which eta and max_depth produces lowest root mean square error
for(i in 1:3){
  for(j in 1:3){
    parameters <- list(eta = eta[i], gamma = 0, eval_metric = 'rmse', max_depth = max_depth[j], 
                       tree_method = 'hist', enable_categorical = TRUE)
    cv <- xgb.cv(params = parameters, data = xgb_training, nrounds = 100, nfold = 5, early_stopping_rounds = 15,
                 verbose = 0) # prevents output from showing
    
    if(cv$evaluation_log$test_rmse_mean[cv$best_iteration] < best_rmse){
      best_parameters <- parameters # chooses eta and max_depth such that rmse is reduced in training
      best_nrounds <- cv$best_iteration # chooses best nrounds
    }
  }
} # we used best_parameters as the parameters in XGBoost


# train boosted tree
xgb_tree <- xgb.train(params = best_parameters, data = xgb_training, nrounds = best_nrounds,
                      watchlist = list(train = xgb_training, test = xgb_testing),
                      print_every_n = 25)

# check how xgboost performs on training data
xgb_predictions_train <- predict(xgb_tree, xgb_training)

# plot predictions vs actual responses for training data
D99 <- data.frame(real_values = hiv_train$logrna, Predictions = xgb_predictions_train)

ggplot(D99, aes(x = Predictions, y = real_values)) + 
  geom_point(alpha = 0.6, aes(color = 'Data points')) +
  geom_smooth(method = "lm", se = FALSE, aes(color = 'Regression line')) + # best fit line
  geom_line(linewidth = 1, data = data.frame(x = range(D99$Predictions), y = range(D99$Predictions)), aes(x=x, y=y, color = 'Identity line')) +
  scale_color_manual(name = 'Legend', values = c('Data points'= 'red', 'Regression line' = 'black', 'Identity line' = 'green')) +         
  labs(
    x = "Predicted viral load", y = "Actual viral load",
    title = 
      "Predicted vs Actual viral loads for the training data (XGBoost)") +
  theme_minimal()

# mse for training data
xgb_mse_train <- mean((xgb_predictions_train  - hiv_train$logrna)^2)



# test performance for testing data
xgb_predictions_test <- predict(xgb_tree, xgb_testing)

# plot predictions vs actual responses for training data
D98 <- data.frame(real_values = hiv_test$logrna, Predictions = xgb_predictions_test)

ggplot(D98, aes(x = Predictions, y = real_values)) + 
  geom_point(alpha = 0.6, aes(color = 'Data points')) +
  geom_smooth(method = "lm", se = FALSE, aes(color = 'Regression line')) + # best fit line
  geom_line(linewidth = 1, data = data.frame(x = range(D98$Predictions), y = range(D98$Predictions)), aes(x=x, y=y, color = 'Identity line')) +
  scale_color_manual(name = 'Legend', values = c('Data points'= 'purple', 'Regression line' = 'black', 'Identity line' = 'green')) +         
  labs(
    x = "Predicted viral load", y = "Actual viral load",
    title = 
      "Predicted vs Actual viral loads for the testing data (XGBoost)") +
  theme_minimal()

# mse for testing data
xgb_mse_test <- mean((xgb_predictions_test - hiv_test$logrna)^2)

# mae for testing data
xgb_mae_test <- mean(abs(xgb_predictions_test - hiv_test$logrna))

# correlation between predictions vs actual responses
cor(hiv_test$logrna, xgb_predictions_test)


# other plots:
# plot the actual boosted trees
xgb.plot.tree(model = xgb_tree, trees = 0) # 1st iteration
xgb.plot.tree(model = xgb_tree, trees = 5) # final iteration

# plot showing the importance of each feature with respect to the tree architecture
importance <- xgb.importance(model = xgb_tree)
xgb.plot.importance(importance)

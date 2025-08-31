library(dplyr)
library(tidyr)
library(keras3)
library(tensorflow)
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

# defining the input and output data for the RNN
hiv_train_x <- hiv_train[ , -5]
hiv_train_y <- hiv_train[ , 5, drop = FALSE]
hiv_test_x <- hiv_test[ , -5]
hiv_test_y <- hiv_test[ , 5, drop = FALSE]

# we need to split the data into dataframes for each patient and rearrange into the 3D form of patient x time x covariates (ie 385 x 6 x 5 for training data, 96 x 6 x 5 for testing data)
# we need matrices instead of dataframes for keras to work properly
hiv_train_x <- hiv_train_x %>% select(-patid, -calwk) %>% data.matrix()
hiv_test_x <- hiv_test_x %>% select(-patid, -calwk) %>% data.matrix()

hiv_train_x <- array(hiv_train_x, dim = c(385, 6, 4))
hiv_test_x <- array(hiv_test_x, dim = c(96, 6, 4))

# same for output data
hiv_train_y <- hiv_train_y %>% data.matrix()
hiv_test_y <- hiv_test_y %>% data.matrix()

hiv_train_y <- array(hiv_train_y, dim = c(385, 6, 1))
hiv_test_y <- array(hiv_test_y, dim = c(96, 6, 1))

# RNN
# we are looking to predict the logrna's across 6 time points, so output layer should have 6 neurons (or 1 output layer with 6 cycles?)
# input layer: we input the covariates nnrti, calwk, txday, cens and trtarm recurrently through the 6 time steps, for each patient
# form RNN
RNN1 <- keras_model_sequential()
RNN1 |> 
  layer_lstm(units = 128, activation = 'relu', recurrent_activation = 'sigmoid',
             input_shape = c(6,4), return_sequences = TRUE) |>
    layer_dense(units = 32, activation = 'relu') |>
  layer_lstm(units = 128, activation = 'relu', recurrent_activation = 'sigmoid',
             return_sequences = TRUE) |>
    layer_dense(units = 32, activation = 'relu')|>
  layer_lstm(units = 128, activation = 'relu', recurrent_activation = 'sigmoid',
             return_sequences = TRUE) |>
    layer_dense(units = 1, activation = 'relu')

# compile RNN (RMSprop commonly used for time-series data, mae useful since we use mse in loss)
RNN1 |> compile(optimizer = 'RMSprop', loss = 'mse', metrics = 'mae')

# train the RNN using training data
output <- RNN1 |> fit(hiv_train_x, hiv_train_y, epochs = 20, batch_size = 32)

# assess training performance by predicting training data
prediction <- RNN1 |> predict(hiv_train_x)

# plot predictions vs actual responses
D3 <- data.frame(real_values = c(hiv_train_y), Predictions = c(prediction))

ggplot(D3, aes(x = Predictions, y = real_values)) + 
  geom_point(alpha = 0.6, aes(color = 'Data points')) +
  geom_smooth(method = "lm", se = FALSE, aes(color = 'Regression line')) + # best fit line
  geom_line(linewidth = 1, data = data.frame(x = range(D3$Predictions), y = range(D3$Predictions)), aes(x=x, y=y, color = 'Identity line')) +
  scale_color_manual(name = 'Legend', values = c('Data points'= 'red', 'Regression line' = 'black', 'Identity line' = 'green')) +         
  labs(
    x = "Predicted viral load", y = "Actual viral load",
    title = 
      "Predicted vs Actual viral loads for the training data (RNN)") +
  theme_minimal()

# compare with real training responses
mse <- mean((prediction - hiv_train_y)^2)



# try other models:
RNN2 <- keras_model_sequential()
RNN2 |> 
  layer_lstm(units = 128, activation = 'relu', recurrent_activation = 'sigmoid',
             input_shape = c(6,4), return_sequences = TRUE) |>
  layer_dense(units = 32, activation = 'relu') |>
  layer_lstm(units = 128, activation = 'relu', recurrent_activation = 'sigmoid',
             return_sequences = TRUE) |>
  layer_dense(units = 32, activation = 'relu')|>
  layer_lstm(units = 128, activation = 'relu', recurrent_activation = 'sigmoid',
             return_sequences = TRUE) |>
  layer_dense(units = 32, activation = 'relu')|>
  layer_lstm(units = 128, activation = 'relu', recurrent_activation = 'sigmoid',
             return_sequences = TRUE) |>
  layer_dense(units = 1, activation = 'relu')


RNN2 |> compile(optimizer = 'RMSprop', loss = 'mse', metrics = 'mae')

output <- RNN2 |> fit(hiv_train_x, hiv_train_y, epochs = 20, batch_size = 32)

prediction2 <- RNN2 |> predict(hiv_train_x)

# plot predictions vs actual responses
D6 <- data.frame(real_values = c(hiv_train_y), Predictions = c(prediction2))

ggplot(D6, aes(x = Predictions, y = real_values)) + 
  geom_point(alpha = 0.6, aes(color = 'Data points')) +
  geom_smooth(method = "lm", se = FALSE, aes(color = 'Regression line')) + # best fit line
  geom_line(linewidth = 1, data = data.frame(x = range(D6$Predictions), y = range(D6$Predictions)), aes(x=x, y=y, color = 'Identity line')) +
  scale_color_manual(name = 'Legend', values = c('Data points'= 'red', 'Regression line' = 'black', 'Identity line' = 'green')) +         
  labs(
    x = "Predicted viral load", y = "Actual viral load",
    title = 
      "Predicted vs Actual viral loads for the training data (RNN)") +
  theme_minimal()

# compare with real training responses
mse2 <- mean((prediction2 - hiv_train_y)^2) 


# other attempted models:

RNN3 <- keras_model_sequential()
RNN3 |> 
  layer_lstm(units = 128, activation = 'relu', recurrent_activation = 'sigmoid',
             input_shape = c(6,4), return_sequences = TRUE) |>
  layer_dense(units = 32, activation = 'relu') |>
  layer_lstm(units = 128, activation = 'relu', recurrent_activation = 'sigmoid',
             return_sequences = TRUE) |>
  layer_dense(units = 32, activation = 'relu')|>
  layer_lstm(units = 128, activation = 'relu', recurrent_activation = 'sigmoid',
             return_sequences = TRUE) |>
  layer_dense(units = 32, activation = 'relu')|>
  layer_lstm(units = 128, activation = 'relu', recurrent_activation = 'sigmoid',
             return_sequences = TRUE) |>
  layer_dense(units = 32, activation = 'relu')|>
  layer_lstm(units = 128, activation = 'relu', recurrent_activation = 'sigmoid',
             return_sequences = TRUE) |>
  layer_dense(units = 1, activation = 'relu')


RNN3 |> compile(optimizer = 'RMSprop', loss = 'mse', metrics = 'mae')

output <- RNN3 |> fit(hiv_train_x, hiv_train_y, epochs = 20, batch_size = 32)
prediction3 <- RNN3 |> predict(hiv_train_x)

# plot predictions vs actual responses
D7 <- data.frame(real_values = c(hiv_train_y), Predictions = c(prediction3))

ggplot(D7, aes(x = Predictions, y = real_values)) + 
  geom_point(alpha = 0.6, aes(color = 'Data points')) +
  geom_smooth(method = "lm", se = FALSE, aes(color = 'Regression line')) + # best fit line
  geom_line(linewidth = 1, data = data.frame(x = range(D7$Predictions), y = range(D7$Predictions)), aes(x=x, y=y, color = 'Identity line')) +
  scale_color_manual(name = 'Legend', values = c('Data points'= 'red', 'Regression line' = 'black', 'Identity line' = 'green')) +         
  labs(
    x = "Predicted viral load", y = "Actual viral load",
    title = 
      "Predicted vs Actual viral loads for the training data (RNN)") +
  theme_minimal()

# compare with real training responses
mse3 <- mean((prediction3 - hiv_train_y)^2) 


# change LSTM units
RNN4 <- keras_model_sequential()
RNN4 |> 
  layer_lstm(units = 256, activation = 'relu', recurrent_activation = 'sigmoid',
             input_shape = c(6,4), return_sequences = TRUE) |>
  layer_dense(units = 32, activation = 'relu') |>
  layer_lstm(units = 256, activation = 'relu', recurrent_activation = 'sigmoid',
             return_sequences = TRUE) |>
  layer_dense(units = 32, activation = 'relu')|>
  layer_lstm(units = 256, activation = 'relu', recurrent_activation = 'sigmoid',
             return_sequences = TRUE) |>
  layer_dense(units = 32, activation = 'relu')|>
  layer_lstm(units = 256, activation = 'relu', recurrent_activation = 'sigmoid',
             return_sequences = TRUE) |>
  layer_dense(units = 1, activation = 'relu')

RNN4 |> compile(optimizer = 'RMSprop', loss = 'mse', metrics = 'mae')

output <- RNN4 |> fit(hiv_train_x, hiv_train_y, epochs = 20, batch_size = 32)
prediction4 <- RNN4 |> predict(hiv_train_x)

# plot predictions vs actual responses
D8 <- data.frame(real_values = c(hiv_train_y), Predictions = c(prediction4))

ggplot(D8, aes(x = Predictions, y = real_values)) + 
  geom_point(alpha = 0.6, aes(color = 'Data points')) +
  geom_smooth(method = "lm", se = FALSE, aes(color = 'Regression line')) + # best fit line
  geom_line(linewidth = 1, data = data.frame(x = range(D8$Predictions), y = range(D8$Predictions)), aes(x=x, y=y, color = 'Identity line')) +
  scale_color_manual(name = 'Legend', values = c('Data points'= 'red', 'Regression line' = 'black', 'Identity line' = 'green')) +         
  labs(
    x = "Predicted viral load", y = "Actual viral load",
    title = 
      "Predicted vs Actual viral loads for the training data (RNN)") +
  theme_minimal()

# compare with real training responses
mse4 <- mean((prediction4 - hiv_train_y)^2)


# increase units in dense layers:
RNN5 <- keras_model_sequential()
RNN5 |> 
  layer_lstm(units = 128, activation = 'relu', recurrent_activation = 'sigmoid',
             input_shape = c(6,4), return_sequences = TRUE) |>
  layer_dense(units = 64, activation = 'relu') |>
  layer_lstm(units = 128, activation = 'relu', recurrent_activation = 'sigmoid',
             return_sequences = TRUE) |>
  layer_dense(units = 64, activation = 'relu')|>
  layer_lstm(units = 128, activation = 'relu', recurrent_activation = 'sigmoid',
             return_sequences = TRUE) |>
  layer_dense(units = 64, activation = 'relu')|>
  layer_lstm(units = 128, activation = 'relu', recurrent_activation = 'sigmoid',
             return_sequences = TRUE) |>
  layer_dense(units = 1, activation = 'relu')


RNN5 |> compile(optimizer = 'RMSprop', loss = 'mse', metrics = 'mae')

output <- RNN5 |> fit(hiv_train_x, hiv_train_y, epochs = 20, batch_size = 32)
prediction5 <- RNN5 |> predict(hiv_train_x)

# plot predictions vs actual responses
D9 <- data.frame(real_values = c(hiv_train_y), Predictions = c(prediction5))

ggplot(D9, aes(x = Predictions, y = real_values)) + 
  geom_point(alpha = 0.6, aes(color = 'Data points')) +
  geom_smooth(method = "lm", se = FALSE, aes(color = 'Regression line')) + # best fit line
  geom_line(linewidth = 1, data = data.frame(x = range(D9$Predictions), y = range(D9$Predictions)), aes(x=x, y=y, color = 'Identity line')) +
  scale_color_manual(name = 'Legend', values = c('Data points'= 'red', 'Regression line' = 'black', 'Identity line' = 'green')) +         
  labs(
    x = "Predicted viral load", y = "Actual viral load",
    title = 
      "Predicted vs Actual viral loads for the training data (RNN)") +
  theme_minimal()

# compare with real training responses
mse5 <- mean((prediction5 - hiv_train_y)^2)

# remove assertion of sigmoid function on LSTM units
RNN6 <- keras_model_sequential()
RNN6 |> 
  layer_lstm(units = 128, activation = 'relu',
             input_shape = c(6,4), return_sequences = TRUE) |>
  layer_dense(units = 32, activation = 'relu') |>
  layer_lstm(units = 128, activation = 'relu',
             return_sequences = TRUE) |>
  layer_dense(units = 32, activation = 'relu')|>
  layer_lstm(units = 128, activation = 'relu',
             return_sequences = TRUE) |>
  layer_dense(units = 32, activation = 'relu')|>
  layer_lstm(units = 128, activation = 'relu',
             return_sequences = TRUE) |>
  layer_dense(units = 1, activation = 'relu')


RNN6 |> compile(optimizer = 'RMSprop', loss = 'mse', metrics = 'mae')

output <- RNN6 |> fit(hiv_train_x, hiv_train_y, epochs = 20, batch_size = 32)
prediction6 <- RNN6 |> predict(hiv_train_x)

# plot predictions vs actual responses
D10 <- data.frame(real_values = c(hiv_train_y), Predictions = c(prediction6))

ggplot(D10, aes(x = Predictions, y = real_values)) + 
  geom_point(alpha = 0.6, aes(color = 'Data points')) +
  geom_smooth(method = "lm", se = FALSE, aes(color = 'Regression line')) + # best fit line
  geom_line(linewidth = 1, data = data.frame(x = range(D10$Predictions), y = range(D10$Predictions)), aes(x=x, y=y, color = 'Identity line')) +
  scale_color_manual(name = 'Legend', values = c('Data points'= 'red', 'Regression line' = 'black', 'Identity line' = 'green')) +         
  labs(
    x = "Predicted viral load", y = "Actual viral load",
    title = 
      "Predicted vs Actual viral loads for the training data (RNN)") +
  theme_minimal()

# compare with real training responses
mse6 <- mean((prediction6 - hiv_train_y)^2)

# conclude that RNN2 is the best model, based on mean square error

# use RNN2 to predict viral loads of test data and compare with true values:
test_prediction <- RNN2 |> predict(hiv_test_x)

D4 <- data.frame(real_values = c(hiv_test_y),
                 Predictions = c(test_prediction))

ggplot(D4, aes(x = Predictions, y = real_values)) + 
  geom_point(alpha = 0.6, aes(color = 'Data points')) +
  geom_smooth(method = "lm", se = FALSE, aes(color = 'Regression line')) + # best fit line
  geom_line(linewidth = 1, data = data.frame(x = range(D4$Predictions), y = range(D4$Predictions)), aes(x=x, y=y, color = 'Identity line')) +
  scale_color_manual(name = 'Legend', values = c('Data points'= 'purple', 'Regression line' = 'black', 'Identity line' = 'green')) +         
  labs(
    x = "Predicted viral load", y = "Actual viral load",
    title = 
      "Predicted vs Actual viral loads for the testing data (RNN)") +
  theme_minimal()

# mse for test data
test_mse <- mean((test_prediction - hiv_test_y)^2)

# mae
test_mae <- mean(abs(hiv_test_y - test_prediction))

# correlation between predicted and actual responses
cor(hiv_test_y, test_prediction)

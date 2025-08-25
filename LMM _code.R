library(nlme)
library(dplyr)
library(ggplot2)
# load and process data
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

hiv <- data.frame(patid, nnrti, calwk, txday, logrna, cens, trtarm)
set.seed(158)
patients <- seq(1:481)
test_ids <- sample(patients, size = 96, replace = FALSE) # choosing the patient ID's for test data
test_ids <- sort(test_ids)
train_ids <- setdiff(patients, test_ids)

hiv_train <- hiv %>% filter(patid %in% train_ids)
hiv_test <- hiv %>% filter(patid %in% test_ids)

# Linear Mixed-Effect Models: choosing the best LMM to use (random slope models)
# independent within subject correlation structure
lme1 <- lme(logrna ~ nnrti + calwk + txday + cens + trtarm,
            data = hiv_train, random = ~ calwk | patid, method = 'ML')
summary(lme1)

# compound symmetry structure for within subject correlation
lme2 <- update(lme1, correlation = corCompSymm(form = ~ 1 | patid))
summary(lme2) # ar and arma not appropriate here

# fewer fixed effects, random slope with independent within subject correlation
lme3 <- lme(logrna ~ trtarm * calwk, data = hiv_train,
            random = ~calwk | patid, method = 'ML')
summary(lme3)

lme4 <- update(lme3, correlation = corCompSymm(form = ~ 1 | patid))
summary(lme4)

anova(lme1, lme2, lme3, lme4) # we choose lme1

# improve parameter estimates by using REML instead of ML
lme5 <- update(lme1, method = 'REML')

# predictions
lme5_pred <- predict(lme5, newdata = hiv_test, level = 0) # level asserts that we only use fixed effects in prediction, since expectation of random components is 0

# plot predictions vs actual responses
D5 <- data.frame(real_values = hiv_test$logrna, Predictions = lme5_pred)

ggplot(D5, aes(x = Predictions, y = real_values)) + 
  geom_point(alpha = 0.6, aes(color = 'Data points')) +
  geom_smooth(method = "lm", se = FALSE, aes(color = 'Regression line')) + # best fit line
  geom_line(linewidth = 1, data = data.frame(x = range(D5$Predictions), y = range(D5$Predictions)), aes(x=x, y=y, color = 'Identity line')) +
  scale_color_manual(name = 'Legend', values = c('Data points'= 'deeppink3', 'Regression line' = 'black', 'Identity line' = 'green')) +         
  labs(
    x = "Predicted viral load", y = "Actual viral load",
    title = 
      "Predicted vs Actual viral loads for the testing data (LMM)") +
  theme_minimal()

# mse
mse_lmm <- mean((lme5_pred - hiv_test$logrna)^2)

# mae
mae_lmm <- mean(abs(lme5_pred - hiv_test$logrna))

# correlation between predicted and actual responses
cor(hiv_test$logrna, lme5_pred)

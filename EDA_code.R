library(ggplot2)
library(dplyr)
data=read.table("Data2.txt",header=TRUE)
head(data)

# plot of logrna as a whole
ggplot(data, aes(x = calwk, y = logrna, group = patid, color = factor(trtarm))) +
  geom_line(size = 1) +
  geom_point(size = 1) +
  facet_wrap(~ trtarm, ncol = 4, labeller = label_both) +  # 3 columns: one for each treatment arm
  labs(title = "Longitudinal Log RNA Viral Load by Treatment Arm",
       x = "Calendar Week", y = "Log10 RNA Viral Load", 
  ) +
  theme_minimal() + theme(legend.position = "none")

# separating data (6 measurements per subject for most subjects, m=481)
time <- c(0,2,4,8,16,24)
calwk <- data$calwk
logrna <- as.numeric(data$logrna)
nnrti <- factor(data$nnrti)
txday <- data$txday
cens <- factor(data$cens)
trtarm <- factor(data$trtarm)

# averages for each trtarm
mean_logrna <- aggregate(logrna ~ trtarm + calwk, data = data, FUN = mean, na.rm = TRUE)
ggplot(mean_logrna, aes(x = calwk, y = logrna, color = factor(trtarm), group = trtarm)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(title = "Mean log10 Viral Load Over Time by Treatment", x = "Week",
       y = "Mean log10 Viral Load", 
       color = "Treatment Arm") + 
  scale_color_discrete(labels = 
                         c("Saquinavir", "Indinavir", "Nelfinavir", "Placebo")) +
  theme_minimal()

# creating new data frames that only contain one treatment arm
trt1 <- NULL 
trt2 <- NULL
trt3 <- NULL
trt4 <- NULL
for(i in 1:length(data$patid)){
  if(data$trtarm[i]==1){
    trt1 <- rbind(trt1, data[i,])
  }
  if(data$trtarm[i]==2){
    trt2 <- rbind(trt2, data[i,])
  }
  if(data$trtarm[i]==3){
    trt3 <- rbind(trt3, data[i,])
  }
  if(data$trtarm[i]==4){
    trt4 <- rbind(trt4, data[i,])
  }
}

# calculate the mean and sd logrna values at each time point for each treatment arm
trt1_means <- aggregate(logrna ~ calwk, data = trt1, FUN = mean)
trt2_means <- aggregate(logrna ~ calwk, data = trt2, FUN = mean)
trt3_means <- aggregate(logrna ~ calwk, data = trt3, FUN = mean)
trt4_means <- aggregate(logrna ~ calwk, data = trt4, FUN = mean)
trt1_sd <- aggregate(logrna ~ calwk, data = trt1, FUN= sd)
trt2_sd <- aggregate(logrna ~ calwk, data = trt2, FUN= sd)
trt3_sd <- aggregate(logrna ~ calwk, data = trt3, FUN= sd)
trt4_sd <- aggregate(logrna ~ calwk, data = trt4, FUN= sd)
trt1_mean_sd <- cbind(trt1_means, trt1_sd[,2])
trt2_mean_sd <- cbind(trt2_means, trt2_sd[,2])
trt3_mean_sd <- cbind(trt3_means, trt3_sd[,2])
trt4_mean_sd <- cbind(trt4_means, trt4_sd[,2])

# plotting each treatment arm data (not included in dissertation)
ggplot(trt1_mean_sd, aes(x=calwk, y=logrna)) + geom_line(size=1) + geom_point(size=1.5) +
  geom_ribbon(aes(ymin = trt1_means[,2] - trt1_sd[,2],
                  ymax = trt1_means[,2] + trt1_sd[,2]),
              alpha = 0.4, fill = "blue") + 
  labs(title = 'Mean log-10 HIV RNA values for patients receiving saquinavir (treatment 1)') + theme_minimal()

ggplot(trt2_mean_sd, aes(x=calwk, y=logrna)) + geom_line(size=1) + geom_point(size=1.5) +
  geom_ribbon(aes(ymin = trt2_means[,2] - trt2_sd[,2],
                  ymax = trt2_means[,2] + trt2_sd[,2]),
              alpha = 0.4, fill = "red") + 
  labs(title = 'Mean log-10 HIV RNA values for patients receiving indinavir (treatment 2)') + theme_minimal()

ggplot(trt3_mean_sd, aes(x=calwk, y=logrna)) + geom_line(size=1) + geom_point(size=1.5) +
  geom_ribbon(aes(ymin = trt3_means[,2] - trt3_sd[,2],
                  ymax = trt3_means[,2] + trt3_sd[,2]),
              alpha = 0.4, fill = "green") + 
  labs(title = 'Mean log-10 HIV RNA values for patients receiving nelfinavir (treatment 3)') + theme_minimal()

ggplot(trt4_mean_sd, aes(x=calwk, y=logrna)) + geom_line(size=1) + geom_point(size=1.5) +
  geom_ribbon(aes(ymin = trt4_means[,2] - trt4_sd[,2],
                  ymax = trt4_means[,2] + trt4_sd[,2]),
              alpha = 0.4, fill = "purple") + 
  labs(title = 'Mean log-10 HIV RNA values for patients receiving the placebo (treatment 4)') + theme_minimal()

# Proportions of groups where HIV cured after week 24
# saquinavir:
cured1 <- sum(trt1$logrna < 3 & trt1$calwk == 24)
infected1 <- sum(trt1$logrna >= 3 & trt1$calwk == 24)
cured_percentage1 <- 100 * cured1/(cured1 + infected1)

# indinavir:
cured2 <- sum(trt2$logrna < 3 & trt2$calwk == 24)
infected2 <- sum(trt2$logrna >= 3 & trt2$calwk == 24) 
cured_percentage2 <- 100 * cured2/(cured2 + infected2)

# nelfanavir:
cured3 <- sum(trt3$logrna < 3 & trt3$calwk == 24)
infected3 <- sum(trt3$logrna >= 3 & trt3$calwk == 24) 
cured_percentage3 <- 100 * cured3/(cured3 + infected3)

# placebo:
cured4 <- sum(trt4$logrna < 3 & trt4$calwk == 24)
infected4 <- sum(trt4$logrna >= 3 & trt4$calwk == 24)
cured_percentage4 <- 100 * cured4/(cured4 + infected4)

# bar chart of patients cured after 24 weeks
df_barchart <- data.frame(Treatment = c("Saquinavir", "Saquinavir", 
                                        "Indinavir", "Indinavir", 
                                        "Nelfinavir", "Nelfinavir",
                                        "Placebo", "Placebo"),
                          State = c("Cured", "Infected",
                                     "Cured", "Infected",
                                     "Cured", "Infected",
                                     "Cured", "Infected"),
                          count = c(cured1, infected1,
                                    cured2, infected2,
                                    cured3, infected3,
                                    cured4, infected4),
                          percentage_cured = c(cured_percentage1, NA,
                                               cured_percentage2, NA,
                                               cured_percentage3, NA,
                                               cured_percentage4, NA))
ggplot(df_barchart, aes(x= Treatment, y = count, fill = State)) + 
  geom_bar(stat = "identity") + 
  labs(title = "Proportion of participants cured from HIV after 24 weeks", 
       x = "Treatment", y = "Number of participants", 
       fill = "Cured or Infected?") + 
  geom_text(aes(label = ifelse(!is.na(percentage_cured), paste0(round(percentage_cured,1), "%"), "")),
            position = position_stack(vjust = 0.5)) +
  theme_minimal()

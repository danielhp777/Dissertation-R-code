students <- paste("Student", 1:5)
ages <- seq(6, 9, by = 0.5)  # 6, 6.5, ..., 9

# Create a data frame of heights
data <- expand.grid(Age = ages, Student = students)
data$Height <- round(rnorm(nrow(data), mean = 120 + 5 * data$Age, sd = 2), 1)

# Base R plot
plot(data$Age, data$Height, type = "n", 
     xlab = "Age (years)", ylab = "Height (cm)", 
     main = "Height Measurements")
colors <- rainbow(length(students))

# Add points and lines per student
for (i in seq_along(students)) {
  subset_data <- subset(data, Student == students[i])
  lines(subset_data$Age, subset_data$Height, col = colors[i], lwd = 2)
  points(subset_data$Age, subset_data$Height, col = colors[i], pch = 16)
}
legend("topleft", legend = students, col = colors, pch = 16, lty = 1, lwd = 2)

set.seed(456)  # For reproducibility
students <- paste("Student", 1:5)
ages <- seq(6, 9, by = 0.5)  # 6, 6.5, ..., 9

# Create a data frame of weights
data_w <- expand.grid(Age = ages, Student = students)
data_w$Weight <- round(rnorm(nrow(data_w), mean = 20 + 2.5 * data_w$Age, sd = 1.5), 1)

# Base R plot
plot(data_w$Age, data_w$Weight, type = "n", 
     xlab = "Age (years)", ylab = "Weight (kg)", 
     main = "Longitudinal Weight Measurements")
colors <- rainbow(length(students))

# Add points and lines per student
for (i in seq_along(students)) {
  subset_data <- subset(data_w, Student == students[i])
  lines(subset_data$Age, subset_data$Weight, col = colors[i], lwd = 2)
  points(subset_data$Age, subset_data$Weight, col = colors[i], pch = 16)
}

legend("topleft", legend = students, col = colors, pch = 16, lty = 1, lwd = 2)

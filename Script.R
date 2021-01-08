##### Install required packages

if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org")
if(!require(dplyr)) install.packages("dplyr", repos = "http://cran.us.r-project.org")
if(!require(readr)) install.packages("readr", repos = "http://cran.us.r-project.org")
if(!require(caret)) install.packages("caret", repos = "http://cran.us.r-project.org")
if(!require(zoo)) install.packages("zoo", repos = "http://cran.us.r-project.org")
if(!require(DMwR)) install.packages("DMwR", repos = "http://cran.us.r-project.org")
if(!require(stats)) install.packages("stats", repos = "http://cran.us.r-project.org")
if(!require(e1071)) install.packages("e1071", repos = "http://cran.us.r-project.org")
if(!require(smotefamily)) install.packages("smotefamily", repos = "http://cran.us.r-project.org")
if(!require(rms)) install.packages("rms", repos = "http://cran.us.r-project.org")
if(!require(gridExtra)) install.packages("gridExtra", repos = "http://cran.us.r-project.org")
if(!require(smoother)) install.packages("smoother", repos = "http://cran.us.r-project.org")
if(!require(GeneCycle)) install.packages("GeneCycle", repos = "http://cran.us.r-project.org")


library(tidyverse)
library(dplyr)
library(readr)
library(caret)
library(zoo)
library(DMwR)
library(stats)
library(e1071)
library(smotefamily)
library(rms)
library(gridExtra)
library(smoother)
library(GeneCycle)

### Load the data
# Load training set and transform into data frame
dl_1 <- "~\\Downloads\\exoTrain.csv"
train_set <- read_csv(file = dl_1)
train_set <- as.data.frame(train_set)

# Load testing set and transform into data frame
dl_2 <- "~\\Downloads\\exoTest.csv"
test_set <- read_csv(file = dl_2)
test_set <- as.data.frame(test_set)

# Turn LABEL column into factor
train_set[,1] <- as.factor(train_set[,1])
test_set[,1] <- as.factor(test_set[,1])

#### Plot two examples of stars with exoplanets and two of stars without any

## Plot the flux over time for a star with an exoplanet (index 2)

# Add a row that accounts for 'time'
with_exo_2 <- rbind(seq(1, 3197), train_set[2,-1])

# Transpose the data frame (turn columns into rows and viceversa)
with_exo_2 <- as.data.frame(t(as.matrix(with_exo_2)))

# Provide names to the columns
names(with_exo_2) <- c("time", "fluxes")

# Create the plot: flux over time for the first thousand moments in time
plot_with_2 <- with_exo_2 %>% ggplot(aes(x=time, y=fluxes)) + 
  geom_line() +
  scale_x_continuous(limits=c(1,1000))+
  scale_y_continuous(limits=c(-400, 300))+
  ggtitle("Star 2 (with exoplanet)")

## Flux over time for a star without an exoplanet (index 151)

# We repeat the same procedure

without_exo_151 <- rbind(seq(1,3197), train_set[151, -1])
without_exo_151 <- as.data.frame(t(as.matrix(without_exo_151)))
names(without_exo_151) <- c(x="time", y="fluxes")
plot_without_151 <- without_exo_151 %>% ggplot(aes(x=time, y=fluxes)) + 
  geom_line() +
  scale_x_continuous(limits=c(1,1000))+
  scale_y_continuous(limits=c(-400,300))+
  ggtitle("Star 151 (without exoplanet)")

## Flux over time for a star with an exoplanet (index 5)

with_exo_5 <- rbind(seq(1, 3197), train_set[5,-1])
with_exo_5 <- as.data.frame(t(as.matrix(with_exo_5)))
names(with_exo_5) <- c("time", "fluxes")
plot_with_5 <- with_exo_5 %>% ggplot(aes(x=time, y=fluxes)) + 
  geom_line() +
  scale_x_continuous(limits=c(1,1000))+
  scale_y_continuous(limits=c(-1200, 800))+
  ggtitle("Star 5 (with exoplanet)")

## Flux over time for a star without an exoplanet (index 220)

without_exo_220 <- rbind(seq(1,3197), train_set[220, -1])
without_exo_220 <- as.data.frame(t(as.matrix(without_exo_220)))
names(without_exo_220) <- c(x="time", y="fluxes")
plot_without_220 <- without_exo_220 %>% ggplot(aes(x=time, y=fluxes)) + 
  geom_line() +
  scale_x_continuous(limits=c(1,1000))+
  scale_y_continuous(limits=c(-400,300))+
  ggtitle("Star 220 (without exoplanet)")

## Include all four plots

grid.arrange(plot_with_2, plot_without_151, plot_with_5, plot_without_220, ncol = 2, nrow = 2)
rm(with_exo_2, with_exo_5, without_exo_151, without_exo_220, plot_with_2, plot_with_5, plot_without_151, plot_without_220)

###############################################################################################
### Use SMOTE Data Augment to balance the training set

## Create new training set by adding new examples of the minority class (LABEL 2)

# Since the SMOTE function produces an object of class "gen_data", we only extract the "data" section
set.seed(1, sample.kind = "Rounding")
SMOTE_train <- smotefamily::SMOTE(train_set[,-1], train_set[,1])$data

# We extract the last column which includes the LABEL
LABEL <- as.factor(SMOTE_train[,3198])

# We add the LABEL as the first column
SMOTE_train <- cbind(LABEL, SMOTE_train[,-3198])


##############################################################################################
### Data Normalization 

## Select subset
set.seed(1, sample.kind = "Rounding")
ind <- createDataPartition(SMOTE_train$LABEL, p = 0.1)$Resample1
subset <- SMOTE_train[ind, ]
rm(ind)

## We calculate the z-score for each flux value, for each star

norm_SMOTE_train <- as.data.frame(scale(as.numeric(subset[1,-1]))) # normalize the flux values for the first star
index <- 1 # set initial index

# Normalize the flux values over time for each star, and add them to the existing norm_SMOTE_train
repeat {
  index <- index + 1
  if (index > nrow(subset)){break}
  print(index)
  scaled_index <- as.data.frame(scale(as.numeric(subset[index, -1])))
  norm_SMOTE_train <- cbind(norm_SMOTE_train, scaled_index)
}
rm(scaled_index, SMOTE_train)


# Since we now have the stars as columns, change rows for columns and add back the LABEL column
norm_SMOTE_train<- as.data.frame(t(as.matrix(norm_SMOTE_train)))
LABEL <- subset[,1]
norm_SMOTE_train <- cbind(LABEL, norm_SMOTE_train)


## Plot raw data vs normalized data
raw_1 <- rbind(seq(1,3197), subset[1, -1])
raw_1 <- as.data.frame(t(as.matrix(raw_1)))
names(raw_1) <- c(x="time", y="fluxes")
plot_raw_1 <- raw_1 %>% ggplot(aes(x=time, y=fluxes)) + 
  geom_line() +
  ggtitle("Star 1 (raw fluxes)")
norm_1 <- rbind(seq(1,3197), norm_SMOTE_train[1, -1])
norm_1 <- as.data.frame(t(as.matrix(norm_1)))
names(norm_1) <- c(x="time", y="fluxes")
plot_norm_1 <- norm_1 %>% ggplot(aes(x=time, y=fluxes)) + 
  geom_line() +
  ggtitle("Star 1 (normalized fluxes)")

grid.arrange(plot_raw_1, plot_norm_1, ncol = 2)
rm(raw_1, norm_1, plot_raw_1, plot_norm_1)


#####################################################
### Gaussian Filter to make our data smoother

index <- 1
smooth_norm_SMOTE_train <- smth.gaussian(as.numeric(norm_SMOTE_train[1,-1]), alpha = 15, tails = TRUE)

repeat{
  index <- index + 1
  if (index > nrow(norm_SMOTE_train)) {break}
  print(index)
  smooth_index <- smth.gaussian(as.numeric(norm_SMOTE_train[index, -1]), alpha = 12, tails = TRUE)
  smooth_norm_SMOTE_train <- as.data.frame(rbind(smooth_norm_SMOTE_train, smooth_index))
}
rm(smooth_index, index)
smooth_norm_SMOTE_train <- cbind(LABEL, smooth_norm_SMOTE_train)

# Plot norm data vs and smooth data

norm_100 <- rbind(seq(1,3197), norm_SMOTE_train[100, -1])
norm_100 <- as.data.frame(t(as.matrix(norm_100)))
names(norm_100) <- c(x="time", y="fluxes")
plot_norm_100 <- norm_100 %>% ggplot(aes(x=time, y=fluxes)) + 
  geom_line() +
  scale_x_continuous(limits = c(0, 1000))+
  scale_y_continuous(limits = c) +
  ggtitle("Star 100 (normalized fluxes)")

smooth_100 <- rbind(seq(1,3197), smooth_norm_SMOTE_train[100, -1])
smooth_100 <- as.data.frame(t(as.matrix(smooth_100)))
names(smooth_100) <- c(x="time", y="fluxes")
plot_smooth_100 <- smooth_100 %>% ggplot(aes(x=time, y=fluxes))+
  geom_line()+
  scale_x_continuous(limits = c(0, 1000))+
  ggtitle("Star 100 (smooth norm fluxes)")

grid.arrange(plot_norm_100, plot_smooth_100, ncol = 2)
rm(norm_100, smooth_100, plot_norm_100, plot_smooth_100)

#############################################################################################
### Perform Fast Fourier Transform on each entry (get values of intensity as a function of frequency)

fourier_train <- fft(as.numeric(smooth_norm_SMOTE_train[1,-1]))
fourier_train <- as.data.frame(GeneCycle::periodogram(as.numeric(fourier_train))$spec)
fourier_train <- as.data.frame(t(as.matrix(fourier_train)))

# Take the frequency values we will use when plotting
freq <- GeneCycle::periodogram(as.numeric(fft(as.numeric(smooth_norm_SMOTE_train[1,-1]))))$freq


# Fourier transform applied to each row
index <- 1
repeat{
  index <- index +1
  if (index>nrow(smooth_norm_SMOTE_train)){break}
  print(index)
  index_fou <- fft(as.numeric(smooth_norm_SMOTE_train[index,-1]))
  index_fou <- as.data.frame(GeneCycle::periodogram(as.numeric(index_fou))$spec)
  index_fou <- as.data.frame(t(as.matrix(index_fou)))
  fourier_train <- rbind(fourier_train, index_fou)
}

LABEL <- smooth_norm_SMOTE_train[,1]
fourier_train <- cbind(LABEL, fourier_train)
rm(index_fou, LABEL)

# Plot intensity as a function of frequency (star 1, with an exoplanet)

data <- as.data.frame(t(as.matrix(fourier_train[1,-1])))
freq <- as.data.frame(freq)
data <- cbind(freq, data)
names(data) <- c("frequency", "intensity")
data %>% ggplot(aes(x=frequency, y=intensity)) +
  geom_line()+
  ggtitle("Intensity(frequency) for a star with an exoplanet") 

# Plot intensity as a function of frequency (star 505, without an exoplanet)

data_2 <- as.data.frame(t(as.matrix(fourier_train[505,-1])))
data_2 <- cbind(freq, data_2)
names(data_2) <- c("frequency", "intensity")
data_2 %>% ggplot(aes(x=frequency, y=intensity)) +
  geom_line()+
  ggtitle("Intensity(frequency) for a star without an exoplanet")

###########################################################################################
### Normalize, smooth and Fourier-transform our test set

## Normalizing our test set

norm_test <- as.data.frame(scale(as.numeric(test_set[1,-1]))) # normalize the flux values for the first star
index <- 1 # set initial index

# Normalize the flux values over time for each star, and add them to the existing norm_SMOTE_train
repeat {
  index <- index + 1
  if (index > nrow(test_set)){break}
  print(index)
  scaled_index <- as.data.frame(scale(as.numeric(test_set[index, -1])))
  norm_test <- cbind(norm_test, scaled_index)
}
rm(scaled_index)


# Since we now have the stars as columns, change rows for columns and add back the LABEL column

norm_test <- as.data.frame(t(as.matrix(norm_test)))
LABEL <- test_set[,1]
norm_test <- cbind(LABEL, norm_test)


## Smoothing our test set

index <- 1
smooth_test <- smth.gaussian(as.numeric(norm_test[1,-1]), alpha = 15, tails = TRUE)

repeat{
  index <- index + 1
  if (index > nrow(norm_test)) {break}
  print(index)
  smooth_index <- smth.gaussian(as.numeric(norm_test[index, -1]), alpha = 12, tails = TRUE)
  smooth_test <- as.data.frame(rbind(smooth_test, smooth_index))
}
rm(smooth_index, index)
smooth_test <- cbind(LABEL, smooth_test)


## Apply the Fourier transform to our test set

fourier_test <- fft(as.numeric(smooth_test[1,-1]))
fourier_test <- as.data.frame(GeneCycle::periodogram(as.numeric(fourier_test))$spec)
fourier_test <- as.data.frame(t(as.matrix(fourier_test)))


index <- 1
repeat{
  index <- index +1
  if (index>nrow(smooth_test)) {break}
  print(index)
  index_fou <- fft(as.numeric(smooth_test[index,-1]))
  index_fou <- as.data.frame(GeneCycle::periodogram(as.numeric(index_fou))$spec)
  index_fou <- as.data.frame(t(as.matrix(index_fou)))
  fourier_test <- rbind(fourier_test, index_fou)
}

LABEL <- smooth_test[,1]
fourier_test <- cbind(LABEL, fourier_test)

rm(index_fou, LABEL, index)



############################################################################################
### Fitting our models 

fit_svm <- svm(LABEL ~., data = fourier_train)
pred_svm <- predict(fit_svm, fourier_test)
svm_results <- confusionMatrix(pred_svm, test_set$LABEL, positive = "2")$table

confusionMatrix(pred_svm, test_set$LABEL, positive="2")

fit_lm <- lm(LABEL~., data = fourier_train)
pred_lm <- predict(fit_lm, fourier_test)

# we give those predicted values equal or above 1.5 a LABEL value of 2, and a value of 1 to those lower than 1.5
pred_lm[which(pred_lm>=1.5)] <- 2
pred_lm[which(pred_lm<1.5)] <- 1
pred_lm <- as.factor(pred_lm)
lm_results <- confusionMatrix(pred_lm, test_set$LABEL, positive = "2")$table


###########################################################################################
### Results
lm_recall <- lm_results[2,2]/(lm_results[2,2] + lm_results[1,2])
lm_precision <- lm_results[2,2] / (lm_results[2,2]+lm_results[2,1])

svm_recall <- svm_results[2,2]/(svm_results[2,2]+svm_results[1,2])
svm_precision <- svm_results[2,2]/(svm_results[2,2]+svm_results[2,1])

success_metrics <- tibble("Method"=c("lm", "svm"), "Recall"=c(lm_recall, svm_recall), "Precision"=c(lm_precision, svm_precision))

success_metrics



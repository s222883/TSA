library(ggplot2)
library(dplyr)
library(ggfortify)

# Load data
df <- read.csv("C:\\Users\\Diego\\Documents\\DTU\\Spring2023\\TSA\\Project\\A1_co2.txt", sep = " ")

test_date <- 2018
df <- df %>% mutate(partition = case_when(year >= test_date ~ "test", 
                                     TRUE ~ "train"))
# Plotting
ggplot(data = df, aes(x = time, y = co2, color = partition )) +
  geom_line(size = 0.5) + labs(x = "Year", y = "Co2") + 
  scale_color_discrete(name="")

ggsave("DTU\\Spring2023\\TSA\\Project\\ej1.png", width=10, height=5.5)

### EXERCISE 1
## Part 1
df_train <- df[df$partition == "train",] 
Y_train <- df_train$co2
N <- nrow(df_train)
T <- seq(1, N)
p <- 12
col1 <- rep(1, nrow(df_train))
col2 <- T
col3 <- sin((2*pi*T)/p)
col4 <- cos((2*pi*T)/p)
col3 <- sin((2*pi/p) * T)
x <- matrix(c(col1, col2, col3, col4), nrow=N, ncol=4, byrow = F)

# Normal equations
theta <- solve(t(x) %*% x)%*%t(x)%*%Y_train
Y_pred <- x %*% theta
df_train["pred"] = Y_pred 

## Part 2: Present the estimated parameters including a measure of uncertainty
eps = Y_train - Y_pred
var(eps)

## Part 3: Plotting
ggplot(data = df_train, aes(x = time)) +
  geom_line(aes(y = co2, color = "Actual"), size = 0.5) + 
  geom_line(aes(y = pred, color = "Prediction"), size = 0.5) + 
  labs(x = "Year", y = "Co2") + 
  scale_color_discrete(name="")

## Part 4 & 5:
# Relaxation algorithm:
# 1. Start with an approximation (given structure)

ro <- cor(eps[1:N-1], eps[2:N])
sigma_f <- function(ro) {
  S <- matrix(ro, nrow = N, ncol = N)
  for (i in 1:N) {
    for (j in i:N) {
      S[i, j] <- S[i, j]^(abs(i - j))
      S[j, i] <- S[i, j]
    }
  }  
  return(S);
}
S <- sigma_f(ro)
# 2. Solve the normal equations 
for (i in 1:5){
  theta = solve(t(x) %*% solve(S,tol = 1e-30) %*% x, tol = 1e-30) %*% (t(x) %*% solve(S, tol = 1e-30) %*% Y_train)
  Y_pred <- x %*% theta
  eps <- Y_train - Y_pred
  ro <- cor(eps[1:N-1], eps[2:N])
  S <- sigma_f(ro)
}

# Part 6
theta_wls <- solve(t(x) %*% solve(S,tol = 1e-30) %*% x, tol = 1e-30) %*% (t(x) %*% solve(S, tol = 1e-30) %*% Y_train)
Y_pred_wls <- x %*% theta_wls
eps_wls <- Y_train - Y_pred_wls
var_wls = (1/(N-4)) * t(eps_wls) %*% solve(S) %*% (eps_wls)
df_train["pred_wls"] = Y_pred_wls 
alpha = 0.05
uncert <- vector(length = 4)
for (i in 1:N){
  uncert[i] <- qt(alpha/2, N-4) * sqrt(var_wls) * sqrt(1+t(x[i,] %*% solve(t(x) %*% x) %*% x[i,]))
}
df_train["uncert"] <- uncert
ggplot(data = df_train, aes(x = time)) +
  geom_line(aes(y = co2, color = "Actual"), size = 0.5) + 
  geom_line(aes(y = pred, color = "Prediction"), size = 0.5) + 
  geom_line(aes(y = pred_wls, color = "Prediction WLS"), size = 0.5) +
  geom_line(aes(y = pred_wls + uncert), size = 0.2) +
  geom_line(aes(y = pred_wls - uncert), size = 0.2) +
  labs(x = "Year", y = "Co2") + 
  scale_color_discrete(name="")

## EXERCISE 3
# 3.1, 3.2, 3.3
# Part 1
# L <- matrix(c(c(1,0,0,0), c(1,1,0,0), c(0,0, cos(2*pi/p), sin(2*pi/p)), c(0,0, -sin(2*pi/p), cos(2*pi/p))), nrow = 4, ncol=4, byrow=T)
# f_0 <- matrix(c(1,0,0,1)) 
f_j <- function(k) {
  return(matrix(c(1, k, sin((2*pi/p)*k), cos((2*pi/p)*k)), nrow = 4))
}
init <- 10
Y_one_step <- 0*vector(length = N)
theta_one_step <- 0*matrix(nrow = N, ncol = 4)
error_one_step <- 0*vector(length = N)
lambda <- 0.9
step <- 1
params <- 4
sigma_one_step <- 0*vector(length = N)
conf_int <- 0*vector(length = N)
var_e_one_step <- 0*vector(length = N)
for (i in init:(N-1)){
  sq <- seq(0,i-1)
  T <- sum(lambda^sq)
  sum_F <- matrix(0, nrow=4, ncol=4)
  sum_h <- matrix(0, nrow=4, ncol=1)
  for (j in 0:(i-1)){
    sum_F <- sum_F + lambda^j * f_j(-j) %*% t(f_j(-j))
    sum_h <- sum_h + lambda^j * f_j(-j) * Y_train[i - j]
  }
  F_N <- sum_F
  h_N <- sum_h
  theta_one_step[i,] <- solve(F_N, tol = 1e-38) %*% h_N
  Y_one_step[i+1] <- t(f_j(step)) %*% theta_one_step[i,]
  error_one_step[i+1] <- Y_train[i+1] - Y_one_step[i+1]
  # Computing the sigma
  sum_sigma <- 0
  for (j in 0:(i-1)){
    sum_sigma <- sum_sigma + lambda^j * (Y_train[i-j] - t(f_j(-j)) %*% theta_one_step[i,])^2
  }
  sigma_one_step[i] <- sum_sigma/(T - params)
  var_e_one_step[i+1] <- sigma_one_step[i] * (1 + (1 + t(f_j(step)) %*% solve(F_N, tol = 1e-38) %*% f_j(step)))
  conf_int[i+1] <- qt(alpha/2, N-params) * sqrt(sigma_one_step[i] * (1 + t(f_j(step)) %*% solve(F_N, tol = 1e-38) %*% f_j(step)))
}
df_train["Prediction"] <- Y_one_step
df_train["Conf_int"] <- conf_int
df_train["Sd"] <- sqrt(sigma_one_step)
df_train["var_e"] <- var_e_one_step

# Plot Prediction
ggplot(data = df_train[init:N,], aes(x = time)) +
  geom_line(aes(y = co2, color = "Actual"), size = 0.5) + 
  geom_line(aes(y = Prediction, color = "One step prediction"), size = 0.5) +
  labs(x = "Year", y = "Co2") +
  scale_color_discrete(name="")

df_train <- df_train %>% mutate(Upper_Bound = Prediction + Conf_int, 
                    Lower_Bound = Prediction - Conf_int) 
df_train <- df_train %>% mutate(Error = co2 - Prediction)

# Plotting
# Plot errors and sds
ggplot(data = df_train[10:N,], aes(x = time)) +
  geom_line(aes(y = Error, color = "Error"))+
  geom_line(aes(y = var_e, color = "Sd"))+
  labs(x = "Year", y = "Error") +
  scale_color_discrete(name="")

ggplot(data = df_train[10:N,], aes(x = time)) +
  geom_line(aes(y = co2, color = "Actual"), size = 0.5) + 
  geom_line(aes(y = One_Step, color = "One step prediction"), size = 0.5) +
  geom_line(aes(y = Upper_Bound, color = "+ CI"), linetype=4, size = 0.3) +
  geom_line(aes(y = Lower_Bound, color = "- CI"), linetype=4, size = 0.3) +
  labs(x = "Year", y = "Co2") + 
  scale_color_discrete(name="")


## Exercise 3.
df_test <- df[df$partition=="test",]
N_test <- nrow(df_test)
steps <- seq(1, N_test)
df_zoom <- df_train[df_train$year>=2010,]
N_train <- nrow(df_zoom)
sum_F <- matrix(0, nrow=4, ncol=4)
sum_h <- matrix(0, nrow=4, ncol=1)
for (j in 0:(N_train-1)){
  sum_F <- sum_F + lambda^j * f_j(-j) %*% t(f_j(-j))
  sum_h <- sum_h + lambda^j * f_j(-j) * df_zoom$co2[N_train - j]
}
F_N <- sum_F
h_N <- sum_h
theta_test <- solve(F_N, tol = 1e-38) %*% h_N
Y_test <- 0*vector(length = N_test)
error_test <- 0*vector(length = N_test)
sigma_test <- 0*vector(length = N_test)
conf_int_test <- 0*vector(length = N_test)
for (i in steps){
  step <- i
  Y_test[i] <- t(f_j(step)) %*% theta_test
  error_test[i] <- df_zoom$co2[i] - Y_test[i]
  # Computing the sigma
  sum_sigma <- 0
  for (j in 0:(N_train-1)){
    sum_sigma <- sum_sigma + lambda^j * (df_zoom$co2[N_train-j] - t(f_j(-step)) %*% theta_test)^2
  }
  sigma_test[i] <- sum_sigma/((N_train+1)-step)
  conf_int_test[i] <- qt(alpha/2, N_train-params) * 
    sqrt(sigma_test[i] * (1 + t(f_j(step)) %*% solve(F_N) %*% f_j(step)))
}
df_test["Prediction"] <- Y_test
df_test["Conf_int"] <- conf_int_test
df_test["Sd"] <- sqrt(sigma_test)
df_test <- df_test %>% mutate(Upper_Bound = Prediction - Conf_int, 
                                Lower_Bound = Prediction + Conf_int) 
df_test <- df_test %>% mutate(Error = abs(co2 - Prediction))


df_full <- dplyr::bind_rows(df_zoom, df_test)
ggplot(data = df_full, aes(x = time)) +
  geom_line(aes(y = co2, color = "Actual"), size = 0.5) + 
  geom_line(aes(y = One_Step, color = "Prediction"), size = 0.5) +
  geom_line(aes(y = Prediction, color = "Prediction"), size = 0.5) +
  geom_line(aes(y = Upper_Bound, color = "Upper Bound"), size = 0.5, linetype=2) +
  geom_line(aes(y = Lower_Bound, color = "Lower Bound"), size = 0.5, linetype=2) +
  labs(x = "Year", y = "Co2") + 
  scale_color_discrete(name="")


## Exercise 4.
df_lambdas <- df_train[df_train$year >= 2002,]
N_lambdas <- nrow(df_lambdas)
lambdas <- seq(0.01,1,length.out=100)
M <- length(lambdas)
lambda_error <- 0*vector(length = N_lambdas)
total_error <- 0*vector(length = M)
Y_train_lambdas <- df_lambdas$co2
Y_lambdas <- 0*vector(length = N)
theta_lambdas <- 0*matrix(nrow = N, ncol = 4)
for (k in 1:M){
  for (i in 100:N_lambdas-1){
    sum_F <- matrix(0,nrow=4, ncol=4)
    sum_h <- matrix(c(0,0,0,0), nrow=4)
    sum_sigma <- 0
    for (j in 0:(i-1)){
      sum_F <- sum_F + lambdas[k]^j * f_j(-j) %*% t(f_j(-j))
      sum_h <- sum_h + lambdas[k]^j * f_j(-j) * Y_train_lambdas[i - j]
    }
    F_N <- sum_F
    h_N <- sum_h
    theta_lambdas[i,] <- solve(F_N,tol = 1e-38) %*% h_N
    Y_lambdas[i+1] <- t(f_j(step)) %*% theta_lambdas[i,]
    lambda_error[i+1] <- Y_train_lambdas[i+1] - Y_lambdas[i+1]
  }
  total_error[k] <- sum(lambda_error^2)
}
plot(lambdas, total_error)
opt_idx <- which.min(total_error)
lambdas[opt_idx]


# F_N <- F_N + lambda^i * f_j(-i) %*% t(f_j(-i))
# h_N <- lambda * solve(L) %*% h_N + f_j(0) * Y_train[i+1]
# theta <- solve(F_N) %*% h_N
# Y_one_step[i+1] <- t(f_j(1)) %*%  theta
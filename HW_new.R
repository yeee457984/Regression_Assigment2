####################### 1.(b) ######################
library(MASS)
library(SIS)
library(Matrix)

n <- 200       # sample size
p <- 1500      # dimension
rho <- 0.7     # correlation
beta <- rep(1, 5)  # true coefficients
nrep <- 1000   # number of repetitions

# Initialize results matrix
results <- matrix(0, nrep, 4)
colnames(results) <- c("L1_norm", "L2_norm", "SPE", "SEN")

# Function to simulate data according to model (1)
simulate_data <- function(n, p, rho, beta) {
  # Create correlation matrix
  Sigma <- matrix(rho, p, p)
  diag(Sigma) <- 1
  Sigma[,5] <- Sigma[5,] <- sqrt(rho)
  Sigma[5,5] <- 1
  
  # Generate X from multivariate normal
  X <- mvrnorm(n, mu=rep(0,p), Sigma)
  
  # Generate y according to model (1)
  epsilon <- rnorm(n)
  y <- X[,1:4] %*% beta[1:4] - 4*sqrt(rho)*X[,5] + epsilon
  
  return(list(X=X, y=y))
}

# Simulation
set.seed(123)
for(i in 1:nrep) {
  data <- simulate_data(n, p, rho, beta)
  sis_fit <- suppressWarnings(
    SIS(data$X, data$y, 
        family="gaussian", 
        penalty="lasso", 
        nsis=floor(n/log(n)), 
        tune="bic")
  )
  
  selected <- sis_fit$ix
  beta_est <- rep(0, p)
  beta_est[selected] <- sis_fit$coef.est
  
  # L1 norm of difference
  results[i,1] <- sum(abs(beta_est[1:5] - c(beta[1:4], -4*sqrt(rho))))
  
  # L2 norm of difference
  results[i,2] <- sqrt(sum((beta_est[1:5] - c(beta[1:4], -4*sqrt(rho)))^2))
  
  # SPE (Sensitivity)
  true_nonzero <- c(1:5)
  results[i,3] <- sum(selected %in% true_nonzero) / length(true_nonzero)
  
  # SEN (Specificity)
  true_zero <- (6:p)
  results[i,4] <- sum(!(true_zero %in% selected)) / length(true_zero)
}

# Print average results
round(colMeans(results), 4)

### Ans ###
# L1_norm   L2_norm    SPE      SEN 
# 5.3987    3.7785   0.8154    0.9790

####################### 1.(c) ######################
n = 200  # 樣本數
p = 40   # 變數個數
beta = c(rep(1, 5), rep(0, p-5))  # 真實的beta值
nsim = 1000  # 模擬次數
lasso_results = matrix(0, nsim, 4)  # 存放lasso的結果
alasso_results = matrix(0, nsim, 4)  # 存放adaptive lasso的結果
colnames(lasso_results) = colnames(alasso_results) = c("L1", "L2", "SPE", "SEN")

for(i in 1:nsim) {
  Sigma = matrix(0, p, p)
  for(j in 1:p) {
    for(k in 1:p) {
      Sigma[j,k] = 0.5^abs(j-k)
    }
  }
  
  X = mvrnorm(n, rep(0,p), Sigma)
  e = rnorm(n, 0, 1)
  y = X %*% beta + e
  
  # Lasso
  lasso = glmnet(X, y, family="gaussian", alpha=1)
  obj1 = -deviance(lasso)
  k = lasso$df
  BIC_lasso = log(n)*k - obj1
  lambda.lasso = which.min(BIC_lasso)
  lasso.beta = as.vector(lasso$beta[,lambda.lasso])
  
  # Adaptive Lasso估計
  beta1 = solve(t(X) %*% X) %*% (t(X) %*% y)
  w3 = 1/abs(as.vector(beta1))
  
  alasso = glmnet(X, y, family="gaussian", alpha=1, penalty.factor=w3)
  obj2 = -deviance(alasso)
  k = alasso$df
  BIC_alasso = log(n)*k - obj2
  lambda.alasso = which.min(BIC_alasso)
  alasso.beta = as.vector(alasso$beta[,lambda.alasso])
  
  # L1 norm
  lasso_results[i,1] = sum(abs(lasso.beta - beta))
  alasso_results[i,1] = sum(abs(alasso.beta - beta))
  
  # L2 norm
  lasso_results[i,2] = sqrt(sum((lasso.beta - beta)^2))
  alasso_results[i,2] = sqrt(sum((alasso.beta - beta)^2))
  
  # SPE (Specificity)
  lasso_results[i,3] = sum((lasso.beta[6:p] == 0)) / (p-5)
  alasso_results[i,3] = sum((alasso.beta[6:p] == 0)) / (p-5)
  
  # SEN (Sensitivity)
  lasso_results[i,4] = sum((lasso.beta[1:5] != 0)) / 5
  alasso_results[i,4] = sum((alasso.beta[1:5] != 0)) / 5
}

lasso_mean = colMeans(lasso_results)
alasso_mean = colMeans(alasso_results)

result_table = rbind(lasso_mean, alasso_mean)
rownames(result_table) = c("lasso", "adaptive lasso")
print(round(result_table, 4))

#                  L1     L2    SPE    SEN
# lasso          0.5218 0.2562 0.9711   1
# adaptive lasso 0.4100 0.2102 0.9927   1

####################### 2.(a) Logistic Regressio ######################
library(glmnet)
library(ggplot2)
data = read.csv("endometria.csv")

model_glm <- glm(case ~ age + gall + hyper + estrogen + drugs,
                 family = binomial(link = "logit"),
                 data = data)

# 顯示詳細的模型摘要
summary_glm <- summary(model_glm)
print("標準邏輯迴歸結果：")
print(summary_glm)


####################### 2.(b) LASSO & Adaptive LASSO ######################
x <- model.matrix(case ~ age + gall + hyper + estrogen + drugs, data)[,-1]
y <- data$case

# LASSO
set.seed(123)
cv_lasso <- cv.glmnet(x, y, family="binomial", alpha=1)
lasso_coef <- coef(cv_lasso, s="lambda.min")

# Adaptive LASSO
# 使用初始邏輯迴歸係數計算權重
weights <- 1/abs(coef(model_glm))[-1]
weights[is.infinite(weights)] <- max(weights[!is.infinite(weights)]) * 100
cv_alasso <- cv.glmnet(x, y, family="binomial", alpha=1, 
                       penalty.factor=weights)
alasso_coef <- coef(cv_alasso, s="lambda.min")

# 係數比較表
all_variables <- unique(c(rownames(lasso_coef), names(coef(model_glm))))
coef_matrix <- matrix(NA, nrow=length(all_variables), ncol=3)
colnames(coef_matrix) <- c("GLM", "LASSO", "Adaptive_LASSO")
rownames(coef_matrix) <- all_variables
coef_matrix[names(coef(model_glm)), "GLM"] <- coef(model_glm)
coef_matrix[rownames(lasso_coef), "LASSO"] <- as.vector(lasso_coef)
coef_matrix[rownames(alasso_coef), "Adaptive_LASSO"] <- as.vector(alasso_coef)
coef_comparison <- as.data.frame(coef_matrix)
coef_comparison$Variable <- rownames(coef_matrix)
print(coef_comparison)


# 變數選擇結果
get_selected_vars <- function(coef_vector, threshold = 1e-4) {
  vars <- rownames(coef_vector)[-1]
  selected <- abs(coef_vector[-1]) > threshold
  return(vars[selected])
}

print("LASSO選擇的變數：")
print(get_selected_vars(lasso_coef)) # gall, estrogen, drugs   

print("Adaptive LASSO選擇的變數：") # gall, estrogen
print(get_selected_vars(alasso_coef))

####################### 2.(c) ######################


#######################  3.(a) Linear and Polynomial Fit ######################
library(locpol)
library(KernSmooth)
data2 <- read.table("wool.txt")
colnames(data2) <- c("xt", "yt")
linear_model <- lm(yt ~ xt, data = data2)
poly_model <- lm(yt ~ poly(xt, 10, raw=TRUE), data = data2)

par(mfrow=c(1,1))
plot(data2$xt, data2$yt,
     xlab = "Time (weeks since Jan 1, 1976)",
     ylab = "Log Price Ratio (19μm price/floor price)",
     main = "Wool Price Analysis: Linear and Polynomial Fits",
     pch = 16,
     cex.lab = 1.2,
     cex.axis = 1.1)
lines(sort(data2$xt), predict(linear_model)[order(data2$xt)], 
      col = "blue", lwd = 2)
lines(sort(data2$xt), predict(poly_model)[order(data2$xt)], 
      col = "red", lwd = 2)
legend("topright", 
       legend = c("Observed Data", "Linear Fit", "Polynomial Fit (degree 10)"),
       col = c("black", "blue", "red"),
       pch = c(16, NA, NA),
       lty = c(NA, 1, 1),
       lwd = c(NA, 2, 2),
       cex = 0.8)

####################### 3.(b) Local Linear Kernel Estimation #######################
bandwidth <- dpill(data2$xt, data2$yt)
print(paste("Optimal bandwidth:", round(bandwidth, 4)))

# 使用locpol進行局部線性擬合
ll_fit <- locpol(yt ~ xt, data = data2, 
                 bw = bandwidth,
                 kernel = gaussK,
                 deg = 1,
                 xeval = data2$xt)

par(mfrow=c(1,1))
plot(data2$xt, data2$yt,
     xlab = "Time (weeks since Jan 1, 1976)",
     ylab = "Log Price Ratio (19μm price/floor price)",
     main = paste0("Local Linear Regression (h = ", round(bandwidth, 4), ")"),
     pch = 16,
     cex.lab = 1.2,
     cex.axis = 1.1)
ci <- confInterval(ll_fit)
lines(ll_fit$lpFit[,1], ll_fit$lpFit[,2], col = "green", lwd = 2)
lines(ci$xgrid, ci$upper, col = "green", lty = 2)
lines(ci$xgrid, ci$lower, col = "green", lty = 2)
legend("topright", 
       legend = c("Observed Data", "Local Linear Fit", "95% Confidence Band"),
       col = c("black", "green", "green"),
       pch = c(16, NA, NA),
       lty = c(NA, 1, 2),
       cex = 0.8)


####################### 3.(c) Compare Fitted Values #######################
linear_rss <- sum((predict(linear_model) - data2$yt)^2)
poly_rss <- sum((predict(poly_model) - data2$yt)^2)

ll_fitted <- fitted(ll_fit)
local_rss <- sum((ll_fitted - data2$yt)^2)

cat("\n模型比較結果：\n")
cat("1. 線性回歸 RSS:", round(linear_rss, 4), "\n")  # 0.9628
cat("2. 10次多項式回歸 RSS:", round(poly_rss, 4), "\n") # 0.1061
cat("3. 局部線性回歸 RSS:", round(local_rss, 4), "\n") # 0.0259
cat("   選擇的帶寬:", round(bandwidth, 4), "\n") # 3.624

improve_poly <- (linear_rss - poly_rss) / linear_rss * 100
improve_local <- (linear_rss - local_rss) / linear_rss * 100
cat("\n相較於線性回歸的改善：\n")
cat("多項式回歸改善了", round(improve_poly, 2), "%\n")
cat("局部線性回歸改善了", round(improve_local, 2), "%\n")

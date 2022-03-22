library(mvtnorm)
library(tidyverse)
source("functions/BayesKalmJoint.R")
source("functions/BayesKalmUneq2.R")
source("functions/ffbsJoint.R")



# set.seed(123)
numits <- 1000
cs <- 2
n <- 100
sigeps <- c(25, 20)
SIGeps <- diag(sigeps)
Sigma <- matrix(c(10, 2, 2, 10), 2)
TT <- 10
B1 <- c(4, 2, 1)
B2 <- c(-3, 0, 1)

B <- cbind(B1, B2)

etas <- rmvnorm(TT * n, sigma = Sigma)


mydatraw <- data.frame(
  id = sort(rep(1:n, TT)),
  time2 = 1:TT,
  eta1 = etas[,1],
  eta2 = etas[,2]
) %>%
  group_by(id) %>%
  mutate(
    
    alpha1 = cumsum(eta1), 
    alpha2 = cumsum(eta2),
    x1 = time2,
    x2 = rbinom(1, 1, .5) * time2,
    oind = c(sample(c(TRUE, rep(FALSE, TT-2))), FALSE),
    oind = any(oind) & time2 >= time2[oind],
    oind = time2*oind,
    x3 = ifelse(oind > 0, time2 - min(oind[oind> 0]), 0),
    Kern1 = B1[1] * x1 + B1[2] * x2 + B1[3] * x3,
    Kern2 = B2[1] * x1 + B2[2] * x2 + B2[3] * x3,
    # Kern1 = B1[1] * x1,
    # Kern2 = B2[1] * x1,
    y1 = alpha1  + Kern1 + rnorm(TT, sd = sqrt(sigeps[1])),
    y2 = alpha2  + Kern2 + rnorm(TT, sd = sqrt(sigeps[2])),
    # Keep = TRUE
    Keep = c(TRUE, sample(c(TRUE, FALSE), TT - 2, replace = TRUE), TRUE)
  ) 

mydat <- mydatraw %>% 
  filter(Keep) %>%
  ungroup()


outcomes <- c("y1", "y2")
predictors <- c("x1", "x2", "x3")
timevar <- "time2"
initialization <- "Bayes"
id <- "id"
silence <- FALSE



bkjoint <- BayesKalmJoint(mydat, outcomes, predictors, timevar, id)




thin <- 1
# burn <- numits/2 + 1
burn <- 1
burnp <- burn/thin 



CI <- quantile(bkjoint$Eta.Track[1,1,seq(burn, numits, thin)], c(.025, .975))
Covered <- CI[1] < Sigma[1,1] & CI[2] > Sigma[1,1]
plot(bkjoint$Eta.Track[1,1,seq(1, numits, thin)], main = paste("ETA11 Covered: ", Covered))
abline(h = CI[1], col = "green")
abline(h = CI[2], col = "green")
abline(h = Sigma[1,1], col = "red")
abline(v = burnp, lty = 2)


CI <- quantile(bkjoint$Eta.Track[1,2,seq(burn, numits, thin)], c(.025, .975))
Covered <- CI[1] < Sigma[1,2] & CI[2] > Sigma[1,2]
plot(bkjoint$Eta.Track[1,2,seq(1, numits, thin)], main = paste("ETA21 Covered: ", Covered))
abline(h = CI[1], col = "green")
abline(h = CI[2], col = "green")
abline(h = Sigma[1,2], col = "red")
abline(v = burnp, lty = 2)


CI <- quantile(bkjoint$Eta.Track[2,2,seq(burn, numits, thin)], c(.025, .975))
Covered <- CI[1] < Sigma[2,2] & CI[2] > Sigma[2,2]
plot(bkjoint$Eta.Track[2,2,seq(1, numits, thin)], main = paste("ETA22 Covered: ", Covered))
abline(h = CI[1], col = "green")
abline(h = CI[2], col = "green")
abline(h = Sigma[2,2], col = "red")
abline(v = burnp, lty = 2)





CI <- quantile(bkjoint$Eps.Track[seq(burn, numits, thin), 1], c(.025, .975))
Covered <- CI[1] < diag(SIGeps)[1] & CI[2] > diag(SIGeps)[1]
plot(bkjoint$Eps.Track[seq(1, numits, thin), 1], main = paste("EPS1 Covered: ", Covered))
abline(h = CI[1], col = "green")
abline(h = CI[2], col = "green")
abline(h = diag(SIGeps)[1], col = "red")
abline(v = burnp, lty = 2)


CI <- quantile(bkjoint$Eps.Track[seq(burn, numits, thin), 2], c(.025, .975))
Covered <- CI[1] < diag(SIGeps)[2] & CI[2] > diag(SIGeps)[2]
plot(bkjoint$Eps.Track[seq(1, numits, thin), 2], main = paste("EPS2 Covered: ", Covered))
abline(h = CI[1], col = "green")
abline(h = CI[2], col = "green")
abline(h = diag(SIGeps)[2], col = "red")
abline(v = burnp, lty = 2)


betadf <- map(1:numits, function(i){
  bdf <- bkjoint[[1]][,,i] %>%
    as.data.frame()
  colnames(bdf) <- paste("C", 1:ncol(bdf), sep = "")
  row.names(bdf) <- paste("B", 1:nrow(bdf), sep = "")
  
  bdf %>% 
    rownames_to_column("Beta") %>%
    gather("outcome", "value", starts_with("C")) %>%
    mutate(Param = paste(Beta, outcome, sep = ""), iteration = i)
}) %>%
  do.call("rbind", .)


betadf %>%
  group_by(Param) %>%
  nest() %>%
  mutate(plots = map2(data, Param,function(df, title){
    
    df <- filter(df, iteration %in% seq(1, numits, thin))
    
    i <- substr(df$Beta[1], 2, nchar(df$Beta[1])) %>% as.numeric()
    j <- substr(df$outcome[1], 2, nchar(df$outcome[1])) %>% as.numeric()
    
    CI <- quantile(df$value[df$iteration >= burn], c(0.025, 0.975))
    B.mean <- mean(df$value[df$iteration >= burn])
    Covered <- CI[1] < B[i,j] & CI[2] > B[i,j]
    
    df %>%
      ggplot(aes(x = iteration, y = value)) + 
      geom_point() +
      geom_hline(yintercept = B[i,j], color = "red") +
      geom_hline(aes(yintercept = B.mean),color = "blue") +
      geom_hline(aes(yintercept = CI[1]),color = "green") +
      geom_hline(aes(yintercept = CI[2]),color = "green") +
      geom_vline(aes(xintercept = burnp), linetype = 2) +
      
      ggtitle(paste(title, "Covered:", Covered))
  })) %>%
  .[["plots"]]





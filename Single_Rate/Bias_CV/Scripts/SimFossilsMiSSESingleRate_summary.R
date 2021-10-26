

library(hisse)
library(tidyverse)
library(Metrics)

userID  <-  Sys.info()['user']


#set dirs based on userID
switch(userID,
"jeremybeaulieu" = { # assume run from dir selon_work/
    out.dir <- "/Users/jeremybeaulieu/hisse_fossil_sims/Single_Rate/Bias_CV/Bias_CV_xls/";
    single.rate.dir <- "/Users/jeremybeaulieu/hisse_fossil_sims/Single_Rate/Bias_CV/SingleRate/scenario";},
"jmbeauli" = { # assume run from dir selon_work/
    out.dir <- "/Users/jeremybeaulieu/hisse_fossil_sims/Single_Rate/Bias_CV/Bias_CV_xls/";
    single.rate.dir <- "/Users/jmbeauli/hisse_fossil_sims/Single_Rate/Bias_CV/SingleRate/scenario";}
)


GetTable <- function(prefix, psi.true){
    
    files <- system(paste("ls -1 ", prefix, ".*", sep=""), intern=TRUE)
    tmp <- c()
    for(i in 1:length(files)){
        load(files[i])
        tmp <- rbind(tmp, local_results)
    }
    tmp <- cbind(tmp, psi.true=psi.true)
    return(tmp)
}


GetVariance <- function(x){
    var <- mean((mean(x) - x)^2)
    return(var)
}


full <- c()
true.lambda <- 0.4
true.mu.vector <- c(0.3, 0.2, 0.1)
true.psi.vector <- c(0.01, 0.05, 0.1)
for (mu_index in seq_along(true.mu.vector)) {
    for (psi_index in seq_along(true.psi.vector)) {
        
        prefix <- paste(single.rate.dir, true.lambda, true.mu.vector[mu_index], true.psi.vector[psi_index], sep=".")
        tmp <- GetTable(prefix, psi.true=true.psi.vector[psi_index])
        for (i in 2:ncol(tmp)) {
            tmp[,i] <- as.numeric(tmp[,i])
        }
        full <- rbind(full, tmp)
    }
}

results <- full
results$error.eps <- as.numeric(results$eps0A)-as.numeric(results$true.eps)
results$error.turnover <- as.numeric(results$turnover0A)-as.numeric(results$true.turnover)
results$true.netdiv <- as.numeric(results$true.lambda - results$true.mu)
results$lambda <- as.numeric(results$turnover0A) / (1 + as.numeric(results$eps0A))
results$mu <- (as.numeric(results$turnover0A) * as.numeric(results$eps0A)) / (1 + as.numeric(results$eps0A))
results$net_div <- as.numeric(results$lambda)-as.numeric(results$mu)
results$error.netdiv <- as.numeric(results$net_div)-as.numeric(results$true.netdiv)
results$m.prop <- as.numeric(results$m.used) / as.numeric(results$m)
results$k.prop <- as.numeric(results$k.used) / as.numeric(results$k)
results$misse <- as.factor(results$misse)
results$psi.error <- as.numeric(results$psi) - as.numeric(results$psi.true)
results_summary <- results %>% group_by(true.mu, psi.true, misse) %>% summarise(mean_turnover_error=mean(error.turnover), median_turnover_error=median(error.turnover), var_turnover=GetVariance(turnover0A), mean_eps_error=mean(error.eps), median_eps_error=median(error.eps), var_eps=GetVariance(eps0A), mean_netdiv_error=mean(error.netdiv), median_netdiv_error=median(error.netdiv), var_netdiv=GetVariance(net_div), m_prop_median=median(m.prop), m_median=median(m.used), k_prop_median=median(k.prop), k_median=median(k.used), mean.psi.error=mean(psi.error), median.psi.error=median(psi.error), var.psi=GetVariance(psi))

file.name <- paste0(out.dir, "TableSummarySingleRate.csv")
write.csv(results_summary, file=file.name)






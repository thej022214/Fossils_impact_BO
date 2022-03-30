

library(hisse)
library(tidyverse)
library(Metrics)
library(viridis)


userID  <-  Sys.info()['user']

#set dirs based on userID
switch(userID,
"jeremybeaulieu" = { # assume run from dir selon_work/
    out.dir <- "/Users/jeremybeaulieu/Fossils_impact_BO/Tables_Figures/";
    out2.dir <- "/Users/jeremybeaulieu/Fossils_impact_BO/MiSSE_two_rates/Bias_CV_xls/";
    scen1.dir <- "/Users/jeremybeaulieu/Fossils_impact_BO/MiSSE_two_rates/Scenario1_FINAL/scenario1.200.";
    scen2.dir <- "/Users/jeremybeaulieu/Fossils_impact_BO/MiSSE_two_rates/Scenario2_FINAL/scenario2.200.";
    scen3.dir <- "/Users/jeremybeaulieu/Fossils_impact_BO/MiSSE_two_rates/Scenario3_FINAL/scenario3.200.";},
"jmbeauli" = { # assume run from dir selon_work/
    out.dir <- "/Users/jmbeauli/Fossils_impact_BO/Tables_Figures/";
    out2.dir <- "/Users/jmbeauli/Fossils_impact_BO/Two_Rate/Bias_CV/";
    scen1.dir <- "/Users/jmbeauli/Fossils_impact_BO/Two_Rate/Scenario1_FINAL/200_taxa/scenario1.200.";
    scen2.dir <- "/Users/jmbeauli/Fossils_impact_BO/Two_Rate/Scenario2_FINAL/200_taxa/scenario2.200.";
    scen3.dir <- "/Users/jmbeauli/Fossils_impact_BO/Two_Rate/Scenario3_FINAL/200_taxa/scenario3.200.";},
)

######################################################################################################################################
######################################################################################################################################
### Some simple utility functions
######################################################################################################################################
######################################################################################################################################

GetTable <- function(prefix){
    
    files <- system(paste("ls -1 ", prefix, ".*", sep=""), intern=TRUE)
    tmp <- c()
    for(i in 1:length(files)){
        load(files[i])
        tmp <- rbind(tmp, local_results)
    }
   return(tmp)
}


GetVariance <- function(x){
    var <- mean((mean(x) - x)^2)
    return(var)
}


BadFileID <- function(x, prefix){
    files <- system(paste("ls -1 ", prefix, ".*", sep=""), intern=TRUE)
    tmp <- c()
    for(i in 1:length(x)){
        tmp <- c(tmp, files[x[i]])
    }
    return(tmp)
}


### CHECKING TO ENSURE ALL 100 REPS ARE IN THEIR RESPECTIVE FOLDERS ###
#index <- "0.05"
#tmp <- c()
#for(i in 1:100){
#    prefix <- paste("scenario4.200", index, i, "Rsave", sep=".")
#    local_results <- NULL
#    try(load(prefix))
#        if(is.null(local_results)){
#        tmp <- c(tmp, i)
#    }
#}


######################################################################################################################################
######################################################################################################################################
### Makes Table of Precision and Accuracy Scenario 1
######################################################################################################################################
######################################################################################################################################

full <- c()
for(index in c(0.01, 0.05, 0.1)){
    prefix <- paste(scen1.dir, index, sep="")
    tmp <- GetTable(prefix)
    for (i in 2:ncol(tmp)) {
        tmp[,i] <- as.numeric(tmp[,i])
    }
    full <- rbind(full, tmp)
}

results <- full
results$error.eps1 <- as.numeric(results$epsA)-as.numeric(results$true.eps1)
results$error.eps2 <- as.numeric(results$epsB)-as.numeric(results$true.eps2)
results$error.turnover1 <- as.numeric(results$turnoverA)-as.numeric(results$true.turnover1)
results$error.turnover2 <- as.numeric(results$turnoverB)-as.numeric(results$true.turnover2)
results$true.netdiv1 <- as.numeric(results$true.lambda1-results$true.mu1)
results$true.netdiv2 <- as.numeric(results$true.lambda2-results$true.mu2)
results$error.netdiv1 <- as.numeric(results$net_divA)-as.numeric(results$true.netdiv1)
results$error.netdiv2 <- as.numeric(results$net_divB)-as.numeric(results$true.netdiv2)
results$m.prop <- as.numeric(results$m.used) / as.numeric(results$m)
results$k.prop <- as.numeric(results$k.used) / as.numeric(results$k)
results$misse <- as.factor(results$misse)
###CHECK THIS###
results$psi.error <- as.numeric(results$psi) - as.numeric(results$true.psi)
################

results_summary <- results %>% group_by(true.psi, misse) %>% summarise(mean_turnover_error1=mean(error.turnover1), median_turnover_error1=median(error.turnover1), var_turnover1=GetVariance(turnoverA-mean(turnoverA)), mean_turnover_error2=mean(error.turnover2), median_turnover_error2=median(error.turnover2), var_turnover2=GetVariance(turnoverB), mean_eps_error1=mean(error.eps1), median_eps_error1=median(error.eps1), var_eps1=GetVariance(epsA), mean_eps_error2=mean(error.eps2), median_eps_error2=median(error.eps2), var_eps2=GetVariance(epsB), mean_netdiv_error1=mean(error.netdiv1), median_netdiv_error1=median(error.netdiv1), var_netdiv1=GetVariance(net_divA), mean_netdiv_error2=mean(error.netdiv2), median_netdiv_error2=median(error.netdiv2), var_netdiv2=GetVariance(net_divB), m_prop_median=median(m.prop), m_median=median(m.used), k_prop_median=median(k.prop), k_median=median(k.used), mean.psi.error=mean(psi.error), median.psi.error=median(psi.error), var.psi=GetVariance(psi))

file.name <- paste0(out2.dir, "Results_summaryScen1.200.csv")
write.csv(results_summary, file=file.name)


######################################################################################################################################
######################################################################################################################################
### Makes Boxplot Scenario 1
######################################################################################################################################
######################################################################################################################################


pdf(paste0(out.dir, "FigureS1.pdf"), width=12, height=12)
layout(matrix(c(1:9), nrow=3, byrow=TRUE))

tmp <- results[results$true.psi == 0.01 & results$misse == "extant",]
turnoverA <- as.numeric(tmp[,"turnoverA"])
turnoverB <- as.numeric(tmp[,"turnoverB"])
netdivA <- as.numeric(tmp[,"net_divA"])
netdivB <- as.numeric(tmp[,"net_divB"])
epsA <- as.numeric(tmp[,"epsA"])
epsB <- as.numeric(tmp[,"epsB"])

boxplot.organized <- data.frame(est=as.numeric(turnoverA), rate="aturnA")
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(turnoverB),rate="aturnB"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(epsA),rate="bepsA"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(epsB),rate="bepsB"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(netdivA),rate="cnet_divA"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(netdivB),rate="cnet_divB"))
boxplot.organized$misse <- rep("extant", dim(boxplot.organized)[1])

psi.set <- c(0.01)
misse.set <- c("full_k", "half_k", "zero_k", "strat_int")
for(psi.index in psi.set){
    for(misse.index in misse.set){
        tmp <- results[results$true.psi == psi.index & results$misse == misse.index,]
        turnoverA <- as.numeric(tmp[,"turnoverA"])
        turnoverB <- as.numeric(tmp[,"turnoverB"])
        netdivA <- as.numeric(tmp[,"net_divA"])
        netdivB <- as.numeric(tmp[,"net_divB"])
        epsA <- as.numeric(tmp[,"epsA"])
        epsB <- as.numeric(tmp[,"epsB"])
        if(misse.index == "strat_int"){
            misse.index = "fstrat_int"
        }
        boxplot.organized.tmp <- data.frame(est=as.numeric(turnoverA), rate="aturnA")
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(turnoverB),rate="aturnB"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(netdivA),rate="cnet_divA"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(netdivB),rate="cnet_divB"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(epsA),rate="bepsA"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(epsB),rate="bepsB"))
        boxplot.organized.tmp$misse <- rep(misse.index, dim(boxplot.organized.tmp)[1])
        boxplot.organized <- rbind(boxplot.organized, boxplot.organized.tmp)
    }
}

### ROW 1 ###
cols <- magma(6)
col.vec <- c(cols[2], cols[3], cols[4], cols[5], 0)
boxplot(boxplot.organized$est~boxplot.organized$misse*boxplot.organized$rate, at=c(1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17, 19,20,21,22,23, 25,26,27,28,29, 31,32,33,34,35), col=c(rep(col.vec,4), rep(rep(0,5),5)), border=c(rep(TRUE,20), rep(FALSE,10)), axes=FALSE, ylim=c(0,2), xlab="", ylab="", main=expression(psi==0.01), lty=1, pch=19, cex=.25, outcol=c(rep("gray", 20), rep(0,10)))
segments(0.5, results$true.turnover1[1], 5.5, results$true.turnover1[1], col="red", lwd=1.5)
segments(6.5, results$true.turnover2[1], 11.5, results$true.turnover2[1], col="red", lwd=1.5)
segments(12.5, results$true.eps1[1], 17.5, results$true.eps1[1], col="red", lwd=1.5)
segments(18.5, results$true.eps2[1], 23.5, results$true.eps2[1], col="red", lwd=1.5)
par(tck=.01)
axis(2, at = seq(0, 2, by = .5), las=1, lwd=1, cex.axis=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = c(3, 9, 15, 21, 27, 33), las=1, lwd=1, cex.axis=1, labels=c(expression(tau[A]), expression(tau[B]), expression(epsilon[A]), expression(epsilon[B]), expression(r[A]), expression(r[B])), mgp=c(.75,.5,0))
title(ylab="Scenario 1", line=2.5, cex=2)
#title(xlab="Parameter", line=2.5)
legend("topleft", c("Extant only", "Stratigraphic ranges", "All fossils used", expression(Half~of~italic(k)~fossils~used), expression(Only~italic(m)~fossils~used)), fill=col.vec, box.col=0)
mtext("a)",side=3, line=0, adj=0, cex=1.5)

par(new=TRUE)
boxplot(boxplot.organized$est~boxplot.organized$misse*boxplot.organized$rate, at=c(1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17, 19,20,21,22,23, 25,26,27,28,29, 31,32,33,34,35), col=c(rep(0,20), col.vec, col.vec), border=c(rep(FALSE,20), rep(TRUE,10)), axes=FALSE, ylim=c(-.12,.12), xlab="", ylab="", lty=1, pch=19, cex=.25, outcol=c(rep(0, 20), rep("gray",10)))
par(tck=.01)
axis(4, at = seq(-.12,0.12, by = .04), las=1, lwd=1, cex.axis=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = c(3, 9, 15, 21, 27, 33), las=1, lwd=1, cex.axis=1, labels=FALSE, mgp=c(.75,.5,0))
segments(24.5, results$true.netdiv1[1], 29.5, results$true.netdiv1[1], col="red", lwd=1.5)
segments(30.5, results$true.netdiv2[1], 35.5, results$true.netdiv2[1], col="red", lwd=1.5)
abline(v=24, col="gray")


tmp <- results[results$true.psi == 0.01 & results$misse == "extant",]
turnoverA <- as.numeric(tmp[,"turnoverA"])
turnoverB <- as.numeric(tmp[,"turnoverB"])
netdivA <- as.numeric(tmp[,"net_divA"])
netdivB <- as.numeric(tmp[,"net_divB"])
epsA <- as.numeric(tmp[,"epsA"])
epsB <- as.numeric(tmp[,"epsB"])

boxplot.organized <- data.frame(est=as.numeric(turnoverA), rate="aturnA")
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(turnoverB),rate="aturnB"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(epsA),rate="bepsA"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(epsB),rate="bepsB"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(netdivA),rate="cnet_divA"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(netdivB),rate="cnet_divB"))
boxplot.organized$misse <- rep("extant", dim(boxplot.organized)[1])

psi.set <- c(0.05)
misse.set <- c("full_k", "half_k", "zero_k", "strat_int")
for(psi.index in psi.set){
    for(misse.index in misse.set){
        tmp <- results[results$true.psi == psi.index & results$misse == misse.index,]
        turnoverA <- as.numeric(tmp[,"turnoverA"])
        turnoverB <- as.numeric(tmp[,"turnoverB"])
        netdivA <- as.numeric(tmp[,"net_divA"])
        netdivB <- as.numeric(tmp[,"net_divB"])
        epsA <- as.numeric(tmp[,"epsA"])
        epsB <- as.numeric(tmp[,"epsB"])
        if(misse.index == "strat_int"){
            misse.index = "fstrat_int"
        }
        boxplot.organized.tmp <- data.frame(est=as.numeric(turnoverA), rate="aturnA")
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(turnoverB),rate="aturnB"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(netdivA),rate="cnet_divA"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(netdivB),rate="cnet_divB"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(epsA),rate="bepsA"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(epsB),rate="bepsB"))
        boxplot.organized.tmp$misse <- rep(misse.index, dim(boxplot.organized.tmp)[1])
        boxplot.organized <- rbind(boxplot.organized, boxplot.organized.tmp)
    }
}


cols <- magma(6)
col.vec <- c(cols[2], cols[3], cols[4], cols[5], 0)
boxplot(boxplot.organized$est~boxplot.organized$misse*boxplot.organized$rate, at=c(1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17, 19,20,21,22,23, 25,26,27,28,29, 31,32,33,34,35), col=c(rep(col.vec,4), rep(rep(0,5),5)), border=c(rep(TRUE,20), rep(FALSE,10)), axes=FALSE, ylim=c(0,2), xlab="", ylab="", main=expression(psi==0.05), lty=1, pch=19, cex=.25, outcol=c(rep("gray", 20), rep(0,10)))
segments(0.5, results$true.turnover1[1], 5.5, results$true.turnover1[1], col="red", lwd=1.5)
segments(6.5, results$true.turnover2[1], 11.5, results$true.turnover2[1], col="red", lwd=1.5)
segments(12.5, results$true.eps1[1], 17.5, results$true.eps1[1], col="red", lwd=1.5)
segments(18.5, results$true.eps2[1], 23.5, results$true.eps2[1], col="red", lwd=1.5)
par(tck=.01)
axis(2, at = seq(0, 2, by = .5), las=1, lwd=1, cex.axis=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = c(3, 9, 15, 21, 27, 33), las=1, lwd=1, cex.axis=1, labels=c(expression(tau[A]), expression(tau[B]), expression(epsilon[A]), expression(epsilon[B]), expression(r[A]), expression(r[B])), mgp=c(.75,.5,0))
mtext("b)",side=3, line=0, adj=0, cex=1.5)

par(new=TRUE)
boxplot(boxplot.organized$est~boxplot.organized$misse*boxplot.organized$rate, at=c(1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17, 19,20,21,22,23, 25,26,27,28,29, 31,32,33,34,35), col=c(rep(0,20), col.vec, col.vec), border=c(rep(FALSE,20), rep(TRUE,10)), axes=FALSE, ylim=c(-.12,.12), xlab="", ylab="", lty=1, pch=19, cex=.25, outcol=c(rep(0, 20), rep("gray",10)))
par(tck=.01)
axis(4, at = seq(-.12,0.12, by = .04), las=1, lwd=1, cex.axis=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = c(3, 9, 15, 21, 27, 33), las=1, lwd=1, cex.axis=1, labels=FALSE, mgp=c(.75,.5,0))
segments(24.5, results$true.netdiv1[1], 29.5, results$true.netdiv1[1], col="red", lwd=1.5)
segments(30.5, results$true.netdiv2[1], 35.5, results$true.netdiv2[1], col="red", lwd=1.5)
abline(v=24, col="gray")


tmp <- results[results$true.psi == 0.01 & results$misse == "extant",]
turnoverA <- as.numeric(tmp[,"turnoverA"])
turnoverB <- as.numeric(tmp[,"turnoverB"])
netdivA <- as.numeric(tmp[,"net_divA"])
netdivB <- as.numeric(tmp[,"net_divB"])
epsA <- as.numeric(tmp[,"epsA"])
epsB <- as.numeric(tmp[,"epsB"])

boxplot.organized <- data.frame(est=as.numeric(turnoverA), rate="aturnA")
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(turnoverB),rate="aturnB"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(epsA),rate="bepsA"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(epsB),rate="bepsB"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(netdivA),rate="cnet_divA"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(netdivB),rate="cnet_divB"))
boxplot.organized$misse <- rep("extant", dim(boxplot.organized)[1])

psi.set <- c(0.1)
misse.set <- c("full_k", "half_k", "zero_k", "strat_int")
for(psi.index in psi.set){
    for(misse.index in misse.set){
        tmp <- results[results$true.psi == psi.index & results$misse == misse.index,]
        turnoverA <- as.numeric(tmp[,"turnoverA"])
        turnoverB <- as.numeric(tmp[,"turnoverB"])
        netdivA <- as.numeric(tmp[,"net_divA"])
        netdivB <- as.numeric(tmp[,"net_divB"])
        epsA <- as.numeric(tmp[,"epsA"])
        epsB <- as.numeric(tmp[,"epsB"])
        if(misse.index == "strat_int"){
            misse.index = "fstrat_int"
        }
        boxplot.organized.tmp <- data.frame(est=as.numeric(turnoverA), rate="aturnA")
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(turnoverB),rate="aturnB"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(netdivA),rate="cnet_divA"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(netdivB),rate="cnet_divB"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(epsA),rate="bepsA"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(epsB),rate="bepsB"))
        boxplot.organized.tmp$misse <- rep(misse.index, dim(boxplot.organized.tmp)[1])
        boxplot.organized <- rbind(boxplot.organized, boxplot.organized.tmp)
    }
}

par(mar=c(5.1, 4.1, 4.1, 3.1))
cols <- magma(6)
col.vec <- c(cols[2], cols[3], cols[4], cols[5], 0)
boxplot(boxplot.organized$est~boxplot.organized$misse*boxplot.organized$rate, at=c(1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17, 19,20,21,22,23, 25,26,27,28,29, 31,32,33,34,35), col=c(rep(col.vec,4), rep(rep(0,5),5)), border=c(rep(TRUE,20), rep(FALSE,10)), axes=FALSE, ylim=c(0,2), xlab="", ylab="", main=expression(psi==0.10), lty=1, pch=19, cex=.25, outcol=c(rep("gray", 20), rep(0,10)))
segments(0.5, results$true.turnover1[1], 5.5, results$true.turnover1[1], col="red", lwd=1.5)
segments(6.5, results$true.turnover2[1], 11.5, results$true.turnover2[1], col="red", lwd=1.5)
segments(12.5, results$true.eps1[1], 17.5, results$true.eps1[1], col="red", lwd=1.5)
segments(18.5, results$true.eps2[1], 23.5, results$true.eps2[1], col="red", lwd=1.5)
par(tck=.01)
axis(2, at = seq(0, 2, by = .5), las=1, lwd=1, cex.axis=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = c(3, 9, 15, 21, 27, 33), las=1, lwd=1, cex.axis=1, labels=c(expression(tau[A]), expression(tau[B]), expression(epsilon[A]), expression(epsilon[B]), expression(r[A]), expression(r[B])), mgp=c(.75,.5,0))
mtext("c)",side=3, line=0, adj=0, cex=1.5)

par(new=TRUE)
boxplot(boxplot.organized$est~boxplot.organized$misse*boxplot.organized$rate, at=c(1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17, 19,20,21,22,23, 25,26,27,28,29, 31,32,33,34,35), col=c(rep(0,20), col.vec, col.vec), border=c(rep(FALSE,20), rep(TRUE,10)), axes=FALSE, ylim=c(-.12,.12), xlab="", ylab="", lty=1, pch=19, cex=.25, outcol=c(rep(0, 20), rep("gray",10)))
par(tck=.01)
axis(4, at = seq(-.12,0.12, by = .04), las=1, lwd=1, cex.axis=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = c(3, 9, 15, 21, 27, 33), las=1, lwd=1, cex.axis=1, labels=FALSE, mgp=c(.75,.5,0))
segments(24.5, results$true.netdiv1[1], 29.5, results$true.netdiv1[1], col="red", lwd=1.5)
segments(30.5, results$true.netdiv2[1], 35.5, results$true.netdiv2[1], col="red", lwd=1.5)
abline(v=24, col="gray")



######################################################################################################################################
######################################################################################################################################
### Makes Table of Precision and Accuracy Scenario 2
######################################################################################################################################
######################################################################################################################################

full <- c()
for(index in c(0.01, 0.05, 0.1)){
    prefix <- paste(scen2.dir, index, sep="")
    tmp <- GetTable(prefix)
    for (i in 2:ncol(tmp)) {
        tmp[,i] <- as.numeric(tmp[,i])
    }
    full <- rbind(full, tmp)
}

results <- full
results$error.eps1 <- as.numeric(results$epsA)-as.numeric(results$true.eps1)
results$error.eps2 <- as.numeric(results$epsB)-as.numeric(results$true.eps2)
results$error.turnover1 <- as.numeric(results$turnoverA)-as.numeric(results$true.turnover1)
results$error.turnover2 <- as.numeric(results$turnoverB)-as.numeric(results$true.turnover2)
results$true.netdiv1 <- as.numeric(results$true.lambda1-results$true.mu1)
results$true.netdiv2 <- as.numeric(results$true.lambda2-results$true.mu2)
results$error.netdiv1 <- as.numeric(results$net_divA)-as.numeric(results$true.netdiv1)
results$error.netdiv2 <- as.numeric(results$net_divB)-as.numeric(results$true.netdiv2)
results$m.prop <- as.numeric(results$m.used) / as.numeric(results$m)
results$k.prop <- as.numeric(results$k.used) / as.numeric(results$k)
results$misse <- as.factor(results$misse)
###CHECK THIS###
results$psi.error <- as.numeric(results$psi) - as.numeric(results$true.psi)
################

results_summary <- results %>% group_by(true.psi, misse) %>% summarise(mean_turnover_error1=mean(error.turnover1), median_turnover_error1=median(error.turnover1), var_turnover1=GetVariance(turnoverA-mean(turnoverA)), mean_turnover_error2=mean(error.turnover2), median_turnover_error2=median(error.turnover2), var_turnover2=GetVariance(turnoverB), mean_eps_error1=mean(error.eps1), median_eps_error1=median(error.eps1), var_eps1=GetVariance(epsA), mean_eps_error2=mean(error.eps2), median_eps_error2=median(error.eps2), var_eps2=GetVariance(epsB), mean_netdiv_error1=mean(error.netdiv1), median_netdiv_error1=median(error.netdiv1), var_netdiv1=GetVariance(net_divA), mean_netdiv_error2=mean(error.netdiv2), median_netdiv_error2=median(error.netdiv2), var_netdiv2=GetVariance(net_divB), m_prop_median=median(m.prop), m_median=median(m.used), k_prop_median=median(k.prop), k_median=median(k.used), mean.psi.error=mean(psi.error), median.psi.error=median(psi.error), var.psi=GetVariance(psi))

file.name <- paste0(out2.dir, "Results_summaryScen2.200.csv")
write.csv(results_summary, file=file.name)


######################################################################################################################################
######################################################################################################################################
### Makes Boxplot Scenario 2
######################################################################################################################################
######################################################################################################################################

tmp <- results[results$true.psi == 0.01 & results$misse == "extant",]
turnoverA <- as.numeric(tmp[,"turnoverA"])
turnoverB <- as.numeric(tmp[,"turnoverB"])
netdivA <- as.numeric(tmp[,"net_divA"])
netdivB <- as.numeric(tmp[,"net_divB"])
epsA <- as.numeric(tmp[,"epsA"])
epsB <- as.numeric(tmp[,"epsB"])

boxplot.organized <- data.frame(est=as.numeric(turnoverA), rate="aturnA")
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(turnoverB),rate="aturnB"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(epsA),rate="bepsA"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(epsB),rate="bepsB"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(netdivA),rate="cnet_divA"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(netdivB),rate="cnet_divB"))
boxplot.organized$misse <- rep("extant", dim(boxplot.organized)[1])

psi.set <- c(0.01)
misse.set <- c("full_k", "half_k", "zero_k", "strat_int")
for(psi.index in psi.set){
    for(misse.index in misse.set){
        tmp <- results[results$true.psi == psi.index & results$misse == misse.index,]
        turnoverA <- as.numeric(tmp[,"turnoverA"])
        turnoverB <- as.numeric(tmp[,"turnoverB"])
        netdivA <- as.numeric(tmp[,"net_divA"])
        netdivB <- as.numeric(tmp[,"net_divB"])
        epsA <- as.numeric(tmp[,"epsA"])
        epsB <- as.numeric(tmp[,"epsB"])
        boxplot.organized.tmp <- data.frame(est=as.numeric(turnoverA), rate="aturnA")
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(turnoverB),rate="aturnB"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(netdivA),rate="cnet_divA"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(netdivB),rate="cnet_divB"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(epsA),rate="bepsA"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(epsB),rate="bepsB"))
        boxplot.organized.tmp$misse <- rep(misse.index, dim(boxplot.organized.tmp)[1])
        boxplot.organized <- rbind(boxplot.organized, boxplot.organized.tmp)
    }
}

### ROW 2 ###
cols <- magma(6)
col.vec <- c(cols[2], cols[3], cols[4], cols[5], 0)
boxplot(boxplot.organized$est~boxplot.organized$misse*boxplot.organized$rate, at=c(1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17, 19,20,21,22,23, 25,26,27,28,29, 31,32,33,34,35), col=c(rep(col.vec,4), rep(rep(0,5),5)), border=c(rep(TRUE,20), rep(FALSE,10)), axes=FALSE, ylim=c(0,2), xlab="", ylab="", main=expression(psi==0.01), lty=1, pch=19, cex=.25, outcol=c(rep("gray", 20), rep(0,10)))
segments(0.5, results$true.turnover1[1], 5.5, results$true.turnover1[1], col="red", lwd=1.5)
segments(6.5, results$true.turnover2[1], 11.5, results$true.turnover2[1], col="red", lwd=1.5)
segments(12.5, results$true.eps1[1], 17.5, results$true.eps1[1], col="red", lwd=1.5)
segments(18.5, results$true.eps2[1], 23.5, results$true.eps2[1], col="red", lwd=1.5)
par(tck=.01)
axis(2, at = seq(0, 2, by = .5), las=1, lwd=1, cex.axis=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = c(3, 9, 15, 21, 27, 33), las=1, lwd=1, cex.axis=1, labels=c(expression(tau[A]), expression(tau[B]), expression(epsilon[A]), expression(epsilon[B]), expression(r[A]), expression(r[B])), mgp=c(.75,.5,0))
title(ylab="Scenario 2", line=2.5, cex=2)
#title(xlab="Parameter", line=2.5)
mtext("d)",side=3, line=0, adj=0, cex=1.5)

par(new=TRUE)
boxplot(boxplot.organized$est~boxplot.organized$misse*boxplot.organized$rate, at=c(1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17, 19,20,21,22,23, 25,26,27,28,29, 31,32,33,34,35), col=c(rep(0,20), col.vec, col.vec), border=c(rep(FALSE,20), rep(TRUE,10)), axes=FALSE, ylim=c(-.12,.12), xlab="", ylab="", lty=1, pch=19, cex=.25, outcol=c(rep(0, 20), rep("gray",10)))
par(tck=.01)
axis(4, at = seq(-.12,0.12, by = .04), las=1, lwd=1, cex.axis=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = c(3, 9, 15, 21, 27, 33), las=1, lwd=1, cex.axis=1, labels=FALSE, mgp=c(.75,.5,0))
segments(24.5, results$true.netdiv1[1], 29.5, results$true.netdiv1[1], col="red", lwd=1.5)
segments(30.5, results$true.netdiv2[1], 35.5, results$true.netdiv2[1], col="red", lwd=1.5)
abline(v=24, col="gray")


tmp <- results[results$true.psi == 0.01 & results$misse == "extant",]
turnoverA <- as.numeric(tmp[,"turnoverA"])
turnoverB <- as.numeric(tmp[,"turnoverB"])
netdivA <- as.numeric(tmp[,"net_divA"])
netdivB <- as.numeric(tmp[,"net_divB"])
epsA <- as.numeric(tmp[,"epsA"])
epsB <- as.numeric(tmp[,"epsB"])


boxplot.organized <- data.frame(est=as.numeric(turnoverA), rate="aturnA")
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(turnoverB),rate="aturnB"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(epsA),rate="bepsA"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(epsB),rate="bepsB"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(netdivA),rate="cnet_divA"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(netdivB),rate="cnet_divB"))
boxplot.organized$misse <- rep("extant", dim(boxplot.organized)[1])

psi.set <- c(0.05)
misse.set <- c("full_k", "half_k", "zero_k", "strat_int")
for(psi.index in psi.set){
    for(misse.index in misse.set){
        tmp <- results[results$true.psi == psi.index & results$misse == misse.index,]
        turnoverA <- as.numeric(tmp[,"turnoverA"])
        turnoverB <- as.numeric(tmp[,"turnoverB"])
        netdivA <- as.numeric(tmp[,"net_divA"])
        netdivB <- as.numeric(tmp[,"net_divB"])
        epsA <- as.numeric(tmp[,"epsA"])
        epsB <- as.numeric(tmp[,"epsB"])
        if(misse.index == "strat_int"){
           misse.index = "fstrat_int"
        }
        boxplot.organized.tmp <- data.frame(est=as.numeric(turnoverA), rate="aturnA")
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(turnoverB),rate="aturnB"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(netdivA),rate="cnet_divA"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(netdivB),rate="cnet_divB"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(epsA),rate="bepsA"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(epsB),rate="bepsB"))
        boxplot.organized.tmp$misse <- rep(misse.index, dim(boxplot.organized.tmp)[1])
        boxplot.organized <- rbind(boxplot.organized, boxplot.organized.tmp)
    }
}

cols <- magma(6)
col.vec <- c(cols[2], cols[3], cols[4], cols[5], 0)
boxplot(boxplot.organized$est~boxplot.organized$misse*boxplot.organized$rate, at=c(1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17, 19,20,21,22,23, 25,26,27,28,29, 31,32,33,34,35), col=c(rep(col.vec,4), rep(rep(0,5),5)), border=c(rep(TRUE,20), rep(FALSE,10)), axes=FALSE, ylim=c(0,2), xlab="", ylab="", main=expression(psi==0.05), lty=1, pch=19, cex=.25, outcol=c(rep("gray", 20), rep(0,10)))
segments(0.5, results$true.turnover1[1], 5.5, results$true.turnover1[1], col="red", lwd=1.5)
segments(6.5, results$true.turnover2[1], 11.5, results$true.turnover2[1], col="red", lwd=1.5)
segments(12.5, results$true.eps1[1], 17.5, results$true.eps1[1], col="red", lwd=1.5)
segments(18.5, results$true.eps2[1], 23.5, results$true.eps2[1], col="red", lwd=1.5)
par(tck=.01)
axis(2, at = seq(0, 2, by = .5), las=1, lwd=1, cex.axis=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = c(3, 9, 15, 21, 27, 33), las=1, lwd=1, cex.axis=1, labels=c(expression(tau[A]), expression(tau[B]), expression(epsilon[A]), expression(epsilon[B]), expression(r[A]), expression(r[B])), mgp=c(.75,.5,0))
mtext("e)",side=3, line=0, adj=0, cex=1.5)

par(new=TRUE)
boxplot(boxplot.organized$est~boxplot.organized$misse*boxplot.organized$rate, at=c(1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17, 19,20,21,22,23, 25,26,27,28,29, 31,32,33,34,35), col=c(rep(0,20), col.vec, col.vec), border=c(rep(FALSE,20), rep(TRUE,10)), axes=FALSE, ylim=c(-.12,.12), xlab="", ylab="", lty=1, pch=19, cex=.25, outcol=c(rep(0, 20), rep("gray",10)))
par(tck=.01)
axis(4, at = seq(-.12,0.12, by = .04), las=1, lwd=1, cex.axis=1, labels=TRUE, mgp=c(.75,.5,0))#
axis(1, at = c(3, 9, 15, 21, 27, 33), las=1, lwd=1, cex.axis=1, labels=FALSE, mgp=c(.75,.5,0))
segments(24.5, results$true.netdiv1[1], 29.5, results$true.netdiv1[1], col="red", lwd=1.5)
segments(30.5, results$true.netdiv2[1], 35.5, results$true.netdiv2[1], col="red", lwd=1.5)
abline(v=24, col="gray")


tmp <- results[results$true.psi == 0.01 & results$misse == "extant",]
turnoverA <- as.numeric(tmp[,"turnoverA"])
turnoverB <- as.numeric(tmp[,"turnoverB"])
netdivA <- as.numeric(tmp[,"net_divA"])
netdivB <- as.numeric(tmp[,"net_divB"])
epsA <- as.numeric(tmp[,"epsA"])
epsB <- as.numeric(tmp[,"epsB"])


boxplot.organized <- data.frame(est=as.numeric(turnoverA), rate="aturnA")
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(turnoverB),rate="aturnB"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(epsA),rate="bepsA"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(epsB),rate="bepsB"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(netdivA),rate="cnet_divA"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(netdivB),rate="cnet_divB"))
boxplot.organized$misse <- rep("extant", dim(boxplot.organized)[1])

psi.set <- c(0.1)
misse.set <- c("full_k", "half_k", "zero_k", "strat_int")
for(psi.index in psi.set){
    for(misse.index in misse.set){
        tmp <- results[results$true.psi == psi.index & results$misse == misse.index,]
        turnoverA <- as.numeric(tmp[,"turnoverA"])
        turnoverB <- as.numeric(tmp[,"turnoverB"])
        netdivA <- as.numeric(tmp[,"net_divA"])
        netdivB <- as.numeric(tmp[,"net_divB"])
        epsA <- as.numeric(tmp[,"epsA"])
        epsB <- as.numeric(tmp[,"epsB"])
        if(misse.index == "strat_int"){
            misse.index = "fstrat_int"
        }
        boxplot.organized.tmp <- data.frame(est=as.numeric(turnoverA), rate="aturnA")
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(turnoverB),rate="aturnB"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(netdivA),rate="cnet_divA"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(netdivB),rate="cnet_divB"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(epsA),rate="bepsA"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(epsB),rate="bepsB"))
        boxplot.organized.tmp$misse <- rep(misse.index, dim(boxplot.organized.tmp)[1])
        boxplot.organized <- rbind(boxplot.organized, boxplot.organized.tmp)
    }
}

### ROW 2 ###
par(mar=c(5.1, 4.1, 4.1, 3.1))
cols <- magma(6)
col.vec <- c(cols[2], cols[3], cols[4], cols[5], 0)
boxplot(boxplot.organized$est~boxplot.organized$misse*boxplot.organized$rate, at=c(1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17, 19,20,21,22,23, 25,26,27,28,29, 31,32,33,34,35), col=c(rep(col.vec,4), rep(rep(0,5),5)), border=c(rep(TRUE,20), rep(FALSE,10)), axes=FALSE, ylim=c(0,2), xlab="", ylab="", main=expression(psi==0.10), lty=1, pch=19, cex=.25, outcol=c(rep("gray", 20), rep(0,10)))
segments(0.5, results$true.turnover1[1], 5.5, results$true.turnover1[1], col="red", lwd=1.5)
segments(6.5, results$true.turnover2[1], 11.5, results$true.turnover2[1], col="red", lwd=1.5)
segments(12.5, results$true.eps1[1], 17.5, results$true.eps1[1], col="red", lwd=1.5)
segments(18.5, results$true.eps2[1], 23.5, results$true.eps2[1], col="red", lwd=1.5)
par(tck=.01)
axis(2, at = seq(0, 2, by = .5), las=1, lwd=1, cex.axis=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = c(3, 9, 15, 21, 27, 33), las=1, lwd=1, cex.axis=1, labels=c(expression(tau[A]), expression(tau[B]), expression(epsilon[A]), expression(epsilon[B]), expression(r[A]), expression(r[B])), mgp=c(.75,.5,0))
mtext("f)",side=3, line=0, adj=0, cex=1.5)

par(new=TRUE)
boxplot(boxplot.organized$est~boxplot.organized$misse*boxplot.organized$rate, at=c(1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17, 19,20,21,22,23, 25,26,27,28,29, 31,32,33,34,35), col=c(rep(0,20), col.vec, col.vec), border=c(rep(FALSE,20), rep(TRUE,10)), axes=FALSE, ylim=c(-.12,.12), xlab="", ylab="", lty=1, pch=19, cex=.25, outcol=c(rep(0, 20), rep("gray",10)))
par(tck=.01)
axis(4, at = seq(-.12,0.12, by = .04), las=1, lwd=1, cex.axis=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = c(3, 9, 15, 21, 27, 33), las=1, lwd=1, cex.axis=1, labels=FALSE, mgp=c(.75,.5,0))
segments(24.5, results$true.netdiv1[1], 29.5, results$true.netdiv1[1], col="red", lwd=1.5)
segments(30.5, results$true.netdiv2[1], 35.5, results$true.netdiv2[1], col="red", lwd=1.5)
abline(v=24, col="gray")



######################################################################################################################################
######################################################################################################################################
### Makes Table of Precision and Accuracy Scenario 3
######################################################################################################################################
######################################################################################################################################

full <- c()
for(index in c(0.01, 0.05, 0.1)){
    prefix <- paste(scen3.dir, index, sep="")
    tmp <- GetTable(prefix)
    for (i in 2:ncol(tmp)) {
        tmp[,i] <- as.numeric(tmp[,i])
    }
    full <- rbind(full, tmp)
}

results <- full
results$error.eps1 <- as.numeric(results$epsA)-as.numeric(results$true.eps1)
results$error.eps2 <- as.numeric(results$epsB)-as.numeric(results$true.eps2)
results$error.turnover1 <- as.numeric(results$turnoverA)-as.numeric(results$true.turnover1)
results$error.turnover2 <- as.numeric(results$turnoverB)-as.numeric(results$true.turnover2)
results$true.netdiv1 <- as.numeric(results$true.lambda1-results$true.mu1)
results$true.netdiv2 <- as.numeric(results$true.lambda2-results$true.mu2)
results$error.netdiv1 <- as.numeric(results$net_divA)-as.numeric(results$true.netdiv1)
results$error.netdiv2 <- as.numeric(results$net_divB)-as.numeric(results$true.netdiv2)
results$m.prop <- as.numeric(results$m.used) / as.numeric(results$m)
results$k.prop <- as.numeric(results$k.used) / as.numeric(results$k)
results$misse <- as.factor(results$misse)
###CHECK THIS###
results$psi.error <- as.numeric(results$psi) - as.numeric(results$true.psi)
################

results_summary <- results %>% group_by(true.psi, misse) %>% summarise(mean_turnover_error1=mean(error.turnover1), median_turnover_error1=median(error.turnover1), var_turnover1=GetVariance(turnoverA-mean(turnoverA)), mean_turnover_error2=mean(error.turnover2), median_turnover_error2=median(error.turnover2), var_turnover2=GetVariance(turnoverB), mean_eps_error1=mean(error.eps1), median_eps_error1=median(error.eps1), var_eps1=GetVariance(epsA), mean_eps_error2=mean(error.eps2), median_eps_error2=median(error.eps2), var_eps2=GetVariance(epsB), mean_netdiv_error1=mean(error.netdiv1), median_netdiv_error1=median(error.netdiv1), var_netdiv1=GetVariance(net_divA), mean_netdiv_error2=mean(error.netdiv2), median_netdiv_error2=median(error.netdiv2), var_netdiv2=GetVariance(net_divB), m_prop_median=median(m.prop), m_median=median(m.used), k_prop_median=median(k.prop), k_median=median(k.used), mean.psi.error=mean(psi.error), median.psi.error=median(psi.error), var.psi=GetVariance(psi))

file.name <- paste0(out2.dir, "Results_summaryScen3.200.csv")
write.csv(results_summary, file=file.name)


######################################################################################################################################
######################################################################################################################################
### Makes Boxplot Scenario 3
######################################################################################################################################
######################################################################################################################################

tmp <- results[results$true.psi == 0.01 & results$misse == "extant",]
turnoverA <- as.numeric(tmp[,"turnoverA"])
turnoverB <- as.numeric(tmp[,"turnoverB"])
netdivA <- as.numeric(tmp[,"net_divA"])
netdivB <- as.numeric(tmp[,"net_divB"])
epsA <- as.numeric(tmp[,"epsA"])
epsB <- as.numeric(tmp[,"epsB"])

boxplot.organized <- data.frame(est=as.numeric(turnoverA), rate="aturnA")
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(turnoverB),rate="aturnB"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(epsA),rate="bepsA"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(epsB),rate="bepsB"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(netdivA),rate="cnet_divA"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(netdivB),rate="cnet_divB"))
boxplot.organized$misse <- rep("extant", dim(boxplot.organized)[1])

psi.set <- c(0.01)
misse.set <- c("full_k", "half_k", "zero_k", "strat_int")
for(psi.index in psi.set){
    for(misse.index in misse.set){
        tmp <- results[results$true.psi == psi.index & results$misse == misse.index,]
        turnoverA <- as.numeric(tmp[,"turnoverA"])
        turnoverB <- as.numeric(tmp[,"turnoverB"])
        netdivA <- as.numeric(tmp[,"net_divA"])
        netdivB <- as.numeric(tmp[,"net_divB"])
        epsA <- as.numeric(tmp[,"epsA"])
        epsB <- as.numeric(tmp[,"epsB"])
        if(misse.index == "strat_int"){
            misse.index = "fstrat_int"
        }
        boxplot.organized.tmp <- data.frame(est=as.numeric(turnoverA), rate="aturnA")
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(turnoverB),rate="aturnB"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(netdivA),rate="cnet_divA"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(netdivB),rate="cnet_divB"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(epsA),rate="bepsA"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(epsB),rate="bepsB"))
        boxplot.organized.tmp$misse <- rep(misse.index, dim(boxplot.organized.tmp)[1])
        boxplot.organized <- rbind(boxplot.organized, boxplot.organized.tmp)
    }
}

### ROW 3 ###
cols <- magma(6)
col.vec <- c(cols[2], cols[3], cols[4], cols[5], 0)
boxplot(boxplot.organized$est~boxplot.organized$misse*boxplot.organized$rate, at=c(1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17, 19,20,21,22,23, 25,26,27,28,29, 31,32,33,34,35), col=c(rep(col.vec,4), rep(rep(0,5),5)), border=c(rep(TRUE,20), rep(FALSE,10)), axes=FALSE, ylim=c(0,2), xlab="", ylab="", main=expression(psi==0.01), lty=1, pch=19, cex=.25, outcol=c(rep("gray", 20), rep(0,10)))
segments(0.5, results$true.turnover1[1], 5.5, results$true.turnover1[1], col="red", lwd=1.5)
segments(6.5, results$true.turnover2[1], 11.5, results$true.turnover2[1], col="red", lwd=1.5)
segments(12.5, results$true.eps1[1], 17.5, results$true.eps1[1], col="red", lwd=1.5)
segments(18.5, results$true.eps2[1], 23.5, results$true.eps2[1], col="red", lwd=1.5)
par(tck=.01)
axis(2, at = seq(0, 2, by = .5), las=1, lwd=1, cex.axis=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = c(3, 9, 15, 21, 27, 33), las=1, lwd=1, cex.axis=1, labels=c(expression(tau[A]), expression(tau[B]), expression(epsilon[A]), expression(epsilon[B]), expression(r[A]), expression(r[B])), mgp=c(.75,.5,0))
title(ylab="Scenario 3", line=2.5, cex=2)
title(xlab="Parameter", line=2.5)
mtext("g)",side=3, line=0, adj=0, cex=1.5)

par(new=TRUE)
boxplot(boxplot.organized$est~boxplot.organized$misse*boxplot.organized$rate, at=c(1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17, 19,20,21,22,23, 25,26,27,28,29, 31,32,33,34,35), col=c(rep(0,20), col.vec, col.vec), border=c(rep(FALSE,20), rep(TRUE,10)), axes=FALSE, ylim=c(-.12,.12), xlab="", ylab="", lty=1, pch=19, cex=.25, outcol=c(rep(0, 20), rep("gray",10)))
par(tck=.01)
axis(4, at = seq(-.12,0.12, by = .04), las=1, lwd=1, cex.axis=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = c(3, 9, 15, 21, 27, 33), las=1, lwd=1, cex.axis=1, labels=FALSE, mgp=c(.75,.5,0))
segments(24.5, results$true.netdiv1[1], 29.5, results$true.netdiv1[1], col="red", lwd=1.5)
segments(30.5, results$true.netdiv2[1], 35.5, results$true.netdiv2[1], col="red", lwd=1.5)
abline(v=24, col="gray")


tmp <- results[results$true.psi == 0.01 & results$misse == "extant",]
turnoverA <- as.numeric(tmp[,"turnoverA"])
turnoverB <- as.numeric(tmp[,"turnoverB"])
netdivA <- as.numeric(tmp[,"net_divA"])
netdivB <- as.numeric(tmp[,"net_divB"])
epsA <- as.numeric(tmp[,"epsA"])
epsB <- as.numeric(tmp[,"epsB"])

boxplot.organized <- data.frame(est=as.numeric(turnoverA), rate="aturnA")
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(turnoverB),rate="aturnB"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(epsA),rate="bepsA"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(epsB),rate="bepsB"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(netdivA),rate="cnet_divA"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(netdivB),rate="cnet_divB"))
boxplot.organized$misse <- rep("extant", dim(boxplot.organized)[1])

psi.set <- c(0.05)
misse.set <- c("full_k", "half_k", "zero_k", "strat_int")
for(psi.index in psi.set){
    for(misse.index in misse.set){
        tmp <- results[results$true.psi == psi.index & results$misse == misse.index,]
        turnoverA <- as.numeric(tmp[,"turnoverA"])
        turnoverB <- as.numeric(tmp[,"turnoverB"])
        netdivA <- as.numeric(tmp[,"net_divA"])
        netdivB <- as.numeric(tmp[,"net_divB"])
        epsA <- as.numeric(tmp[,"epsA"])
        epsB <- as.numeric(tmp[,"epsB"])
        if(misse.index == "strat_int"){
           misse.index = "fstrat_int"
        }
        boxplot.organized.tmp <- data.frame(est=as.numeric(turnoverA), rate="aturnA")
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(turnoverB),rate="aturnB"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(netdivA),rate="cnet_divA"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(netdivB),rate="cnet_divB"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(epsA),rate="bepsA"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(epsB),rate="bepsB"))
        boxplot.organized.tmp$misse <- rep(misse.index, dim(boxplot.organized.tmp)[1])
        boxplot.organized <- rbind(boxplot.organized, boxplot.organized.tmp)
    }
}


cols <- magma(6)
col.vec <- c(cols[2], cols[3], cols[4], cols[5], 0)
boxplot(boxplot.organized$est~boxplot.organized$misse*boxplot.organized$rate, at=c(1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17, 19,20,21,22,23, 25,26,27,28,29, 31,32,33,34,35), col=c(rep(col.vec,4), rep(rep(0,5),5)), border=c(rep(TRUE,20), rep(FALSE,10)), axes=FALSE, ylim=c(0,2), xlab="", ylab="", main=expression(psi==0.05), lty=1, pch=19, cex=.25, outcol=c(rep("gray", 20), rep(0,10)))
segments(0.5, results$true.turnover1[1], 5.5, results$true.turnover1[1], col="red", lwd=1.5)
segments(6.5, results$true.turnover2[1], 11.5, results$true.turnover2[1], col="red", lwd=1.5)
segments(12.5, results$true.eps1[1], 17.5, results$true.eps1[1], col="red", lwd=1.5)
segments(18.5, results$true.eps2[1], 23.5, results$true.eps2[1], col="red", lwd=1.5)
par(tck=.01)
axis(2, at = seq(0, 2, by = .5), las=1, lwd=1, cex.axis=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = c(3, 9, 15, 21, 27, 33), las=1, lwd=1, cex.axis=1, labels=c(expression(tau[A]), expression(tau[B]), expression(epsilon[A]), expression(epsilon[B]), expression(r[A]), expression(r[B])), mgp=c(.75,.5,0))
title(xlab="Parameter", line=2.5)
mtext("h)",side=3, line=0, adj=0, cex=1.5)

par(new=TRUE)
boxplot(boxplot.organized$est~boxplot.organized$misse*boxplot.organized$rate, at=c(1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17, 19,20,21,22,23, 25,26,27,28,29, 31,32,33,34,35), col=c(rep(0,20), col.vec, col.vec), border=c(rep(FALSE,20), rep(TRUE,10)), axes=FALSE, ylim=c(-.12,.12), xlab="", ylab="", lty=1, pch=19, cex=.25, outcol=c(rep(0, 20), rep("gray",10)))
par(tck=.01)
axis(4, at = seq(-.12,0.12, by = .04), las=1, lwd=1, cex.axis=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = c(3, 9, 15, 21, 27, 33), las=1, lwd=1, cex.axis=1, labels=FALSE, mgp=c(.75,.5,0))
segments(24.5, results$true.netdiv1[1], 29.5, results$true.netdiv1[1], col="red", lwd=1.5)
segments(30.5, results$true.netdiv2[1], 35.5, results$true.netdiv2[1], col="red", lwd=1.5)
abline(v=24, col="gray")


tmp <- results[results$true.psi == 0.01 & results$misse == "extant",]
turnoverA <- as.numeric(tmp[,"turnoverA"])
turnoverB <- as.numeric(tmp[,"turnoverB"])
netdivA <- as.numeric(tmp[,"net_divA"])
netdivB <- as.numeric(tmp[,"net_divB"])
epsA <- as.numeric(tmp[,"epsA"])
epsB <- as.numeric(tmp[,"epsB"])

boxplot.organized <- data.frame(est=as.numeric(turnoverA), rate="aturnA")
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(turnoverB),rate="aturnB"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(epsA),rate="bepsA"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(epsB),rate="bepsB"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(netdivA),rate="cnet_divA"))
boxplot.organized <- rbind(boxplot.organized, data.frame(est=as.numeric(netdivB),rate="cnet_divB"))
boxplot.organized$misse <- rep("extant", dim(boxplot.organized)[1])

psi.set <- c(0.1)
misse.set <- c("full_k", "half_k", "zero_k", "strat_int")
for(psi.index in psi.set){
    for(misse.index in misse.set){
        tmp <- results[results$true.psi == psi.index & results$misse == misse.index,]
        turnoverA <- as.numeric(tmp[,"turnoverA"])
        turnoverB <- as.numeric(tmp[,"turnoverB"])
        netdivA <- as.numeric(tmp[,"net_divA"])
        netdivB <- as.numeric(tmp[,"net_divB"])
        epsA <- as.numeric(tmp[,"epsA"])
        epsB <- as.numeric(tmp[,"epsB"])
        if(misse.index == "strat_int"){
            misse.index = "fstrat_int"
        }
        boxplot.organized.tmp <- data.frame(est=as.numeric(turnoverA), rate="aturnA")
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(turnoverB),rate="aturnB"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(netdivA),rate="cnet_divA"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(netdivB),rate="cnet_divB"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(epsA),rate="bepsA"))
        boxplot.organized.tmp <- rbind(boxplot.organized.tmp, data.frame(est=as.numeric(epsB),rate="bepsB"))
        boxplot.organized.tmp$misse <- rep(misse.index, dim(boxplot.organized.tmp)[1])
        boxplot.organized <- rbind(boxplot.organized, boxplot.organized.tmp)
    }
}

par(mar=c(5.1, 4.1, 4.1, 3.1))
cols <- magma(6)
col.vec <- c(cols[2], cols[3], cols[4], cols[5], 0)
boxplot(boxplot.organized$est~boxplot.organized$misse*boxplot.organized$rate, at=c(1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17, 19,20,21,22,23, 25,26,27,28,29, 31,32,33,34,35), col=c(rep(col.vec,4), rep(rep(0,5),5)), border=c(rep(TRUE,20), rep(FALSE,10)), axes=FALSE, ylim=c(0,2), xlab="", ylab="", main=expression(psi==0.10), lty=1, pch=19, cex=.25, outcol=c(rep("gray", 20), rep(0,10)))
segments(0.5, results$true.turnover1[1], 5.5, results$true.turnover1[1], col="red", lwd=1.5)
segments(6.5, results$true.turnover2[1], 11.5, results$true.turnover2[1], col="red", lwd=1.5)
segments(12.5, results$true.eps1[1], 17.5, results$true.eps1[1], col="red", lwd=1.5)
segments(18.5, results$true.eps2[1], 23.5, results$true.eps2[1], col="red", lwd=1.5)
par(tck=.01)
axis(2, at = seq(0, 2, by = .5), las=1, lwd=1, cex.axis=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = c(3, 9, 15, 21, 27, 33), las=1, lwd=1, cex.axis=1, labels=c(expression(tau[A]), expression(tau[B]), expression(epsilon[A]), expression(epsilon[B]), expression(r[A]), expression(r[B])), mgp=c(.75,.5,0))
#title(ylab="Scenario 2", line=2.5, cex=2)
title(xlab="Parameter", line=2.5)
mtext("i)",side=3, line=0, adj=0, cex=1.5)

par(new=TRUE)
boxplot(boxplot.organized$est~boxplot.organized$misse*boxplot.organized$rate, at=c(1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17, 19,20,21,22,23, 25,26,27,28,29, 31,32,33,34,35), col=c(rep(0,20), col.vec, col.vec), border=c(rep(FALSE,20), rep(TRUE,10)), axes=FALSE, ylim=c(-.12,.12), xlab="", ylab="", lty=1, pch=19, cex=.25, outcol=c(rep(0, 20), rep("gray",10)))
par(tck=.01)
axis(4, at = seq(-.12,0.12, by = .04), las=1, lwd=1, cex.axis=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = c(3, 9, 15, 21, 27, 33), las=1, lwd=1, cex.axis=1, labels=FALSE, mgp=c(.75,.5,0))
segments(24.5, results$true.netdiv1[1], 29.5, results$true.netdiv1[1], col="red", lwd=1.5)
segments(30.5, results$true.netdiv2[1], 35.5, results$true.netdiv2[1], col="red", lwd=1.5)
abline(v=24, col="gray")

dev.off()







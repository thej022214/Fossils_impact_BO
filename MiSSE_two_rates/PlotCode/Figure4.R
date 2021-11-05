

library(hisse)
library(picante)
library(phytools)
library(viridis)
library(paleotree)

userID  <-  Sys.info()['user']

#set dirs based on userID
switch(userID,
"jeremybeaulieu" = { # assume run from dir selon_work/
    out.dir <- "/Users/jeremybeaulieu/Fossils_impact_BO/Tables_Figures/";
    out2.dir <- "/Users/jeremybeaulieu/Fossils_impact_BO/MiSSE_two_rates/Bias_CV_xls/";
    scen4.dir <- "/Users/jeremybeaulieu/Fossils_impact_BO/MiSSE_two_rates/Scenario4_FINAL/";},
"jmbeauli" = { # assume run from dir selon_work/
    out.dir <- "/Users/jmbeauli/Fossils_impact_BO/Tables_Figures/";
    out2.dir <- "/Users/jmbeauli/Fossils_impact_BO/Two_Rate/Bias_CV/";
    scen4.dir <- "/Users/jmbeauli/Fossils_impact_BO/Two_Rate/Scenario1_FINAL/200_taxa/scenario1.200.";},
)



######################################################################################################################################
######################################################################################################################################
### The set of 10 random trees
######################################################################################################################################
######################################################################################################################################

set.seed(42)
max.age.true <- c()
max.age.sampled <- c()
max.age.extant <- c()
transitions.raw <- c()
ten.percent <- c()
state.freqs <- c()
ntax=200
true.psi=0.05
tree.set <- as.list(1:2)

for(index in 1:2){
    print(index)
    #Realistic with mass extinction and heritable sampling
    pars <- c(0.5, 0.2, 0.7, 0.9)
    turnover <- c(0.3, 0.6, 1.5, 1.0)
    eps <- c(0.7, 0.9, 0.95, 0.8)
    transition <- 0.1*matrix(c(NA, 0.5, 0.2, 0.6, 0.1, NA, 0.3, 0.2, 0.1, 0.4, NA, 0.8, 0.1, 0.2, 0.3, NA),nrow=4)
    mass.extinction.heights <- c(40,70)
    mass.extinction.magnitudes <- c(0.7, 0.9)
    
    sim.tab <- list()
    sim.tab$n.surviving = 0
    while(sim.tab$n.surviving==0) {
        sim.tab <- hisse::SimulateHisse(turnover.rates=turnover, eps.values=eps, transition.rates=transition, max.taxa=ntax/.75, x0=0, mass.extinction.heights=mass.extinction.heights, mass.extinction.magnitudes=mass.extinction.magnitudes,nstart=2)
    }
    phy <- hisse::SimToPhylo(sim.tab, include.extinct=TRUE)
    phy.full <- geiger::drop.extinct(phy, tol=.Machine$double.eps^.50)
    
    #    # now biased samples down to ntax extant taxa
    #sample.prob.raw <- abs(geiger::sim.char(phy.full, par=1, model="BM"))[,1,1]
    #sample.prob.raw[which(is.na(sample.prob.raw))] <- median(sample.prob.raw, na.rm=TRUE )
    #sample.prob <- sample.prob.raw/max(sample.prob.raw)
    #while(ape::Ntip(geiger::drop.extinct(phy, tol=.Machine$double.eps^.50))>ntax) {
    #    to.kill <- sample(names(sample.prob), size=1, prob=sample.prob)
    #    phy <- ape::drop.tip(phy, to.kill)
    #    phy$tip.state <- phy$tip.state[-which(names(phy$tip.state)==to.kill)]
    #    sample.prob <- sample.prob[-which(names(sample.prob)==to.kill)]
    #}
    if(Ntip(phy)>ntax){
        f <- hisse:::GetFossils(phy, psi=true.psi)
        if(sum(f$has_sampled_descendant) > 0){
            try(test.tol <- geiger::drop.extinct(phy, tol=.Machine$double.eps^.50))
            if(!is.null(test.tol)){
                pp <- hisse:::ProcessSimSample(phy, f)
                sf <- 0.75
                all.dead <- FALSE
                max.age.true <- c(max.age.true, max(dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))[-c(1:Ntip(phy))]))
                max.age.sampled <- c(max.age.sampled, max(dateNodes(pp$phy, rootAge=max(node.depth.edgelength(pp$phy)))[-c(1:Ntip(pp$phy))]))
                max.age.extant <- c(max.age.extant, max(dateNodes(test.tol, rootAge=max(node.depth.edgelength(test.tol)))[-c(1:Ntip(test.tol))]))
                transitions.raw <- rbind(transitions.raw, sim.tab$transition.counts[!is.na(sim.tab$transition.counts)])
                ten.percent <- c(ten.percent, sum(sim.tab$transition.counts[!is.na(sim.tab$transition.counts)]) / length(phy$edge.length))
                true.tip.states <- phy$tip.state[names(phy$tip.state) %in% test.tol$tip.label]
                state.freqs <- rbind(state.freqs, table(true.tip.states))
                #quartz()
                #ltt.plot(phy)
                tree.set[[index]] <- phy
            }
        }
    }
}


######################################################################################################################################
######################################################################################################################################
### Tree and Fossil Portion of the Figure
######################################################################################################################################
######################################################################################################################################

### Sim tree figure:
#phy <- tree.set[[2]]
#f <- hisse:::GetFossils(phy, psi=0.05)
#phy.full <- geiger::drop.extinct(phy, tol=.Machine$double.eps^.50)
#sample.prob.raw <- abs(geiger::sim.char(phy.full, par=1, model="BM"))[,1,1]
#sample.prob.raw[which(is.na(sample.prob.raw))] <- median(sample.prob.raw, na.rm=TRUE )
#sample.prob <- sample.prob.raw/max(sample.prob.raw)
#while(ape::Ntip(geiger::drop.extinct(phy, tol=.Machine$double.eps^.50))>ntax) {
#    to.kill <- sample(names(sample.prob), size=1, prob=sample.prob)
#    phy <- ape::drop.tip(phy, to.kill)
#    phy$tip.state <- phy$tip.state[-which(names(phy$tip.state)==to.kill)]
#    sample.prob <- sample.prob[-which(names(sample.prob)==to.kill)]
#}

#phy.extant.only <- geiger::drop.extinct(phy, tol=.Machine$double.eps^.50)
#phy <- tree.set[[2]]
#trait <- cbind(phy$tip.label, 0)
#trait[which(phy$tip.label %in% phy.extant.only$tip.label),2] <- 1
#trait.final <- data.frame(taxon=trait[,1], trait <- as.numeric(trait[,2]))

#rownames(trait.final) <- trait.final[,1]
#phy0 <- picante::match.phylo.data(phy, trait.final)$phy
#trait0 <- picante::match.phylo.data(phy, trait.final)$data
#mode <- trait0[,2]
#names(mode) <- trait0[,1]


### Split the table into m and k for the points ###
#extinct.samples <- f[which(f$fossiltype_long=="extinct_terminal" | f$fossiltype_long=="extinct_internal"),]
#k.samples.tmp <- extinct.samples[which(extinct.samples$has_sampled_descendant == TRUE),]
#extinct.samples <- extinct.samples[which(extinct.samples$has_sampled_descendant == FALSE),]

#k.samples <- f[which(f$fossiltype_long == "surviving_terminal" | f$fossiltype_long == "surviving_internal"),]
#k.samples <- rbind(k.samples, k.samples.tmp)
###################################################

#par(mfcol=c(2,3),mar=c(4,4.5,0.5,0.5), oma=c(1.5,2,1,1))

#colors_states <- c("transparent", viridis(4)[1])
# Plot states at tips
#plot(ladderize(phy, right=FALSE), show.tip.label=F, edge.width=0.3, adj=1, cex=0.08)
#abline(v=40, col="red", lty=3, lwd=1)
#abline(v=70, col="red", lty=3, lwd=1)
#par(fg="transparent")
#AddFossilPoints(f, ladderize(phy), pch=19, cex=0.8, col=viridis(4)[4])
#hisse:::AddFossilPoints(ladderize(phy, right=FALSE), extinct.samples, pch=19, cex=0.8, col="#0D79F2")
#hisse:::AddFossilPoints(ladderize(phy, right=FALSE), k.samples, pch=19, cex=0.8, col="#F25E0D")
#tiplabels(pie=to.matrix(mode, sort(unique(mode))), piecol=colors_states, cex=0.25, lwd=0.3, frame = "n")
#par(fg="black")
#tiplabels(pch=1, bg=tip.cols, adj=1, cex=0.1, width = 0.1)
#axisPhylo(cex.axis=.75, lwd=.75)

######################################################################################################################################
######################################################################################################################################
### LTT Portion of the Figure -- PANEL A
######################################################################################################################################
######################################################################################################################################

pdf(paste0(out.dir, "Figure4.pdf"), width=12, height=4)
layout(matrix(c(1:3), nrow=1, byrow=TRUE))

phy <- tree.set[[2]]
ss <-ape::ltt.plot.coords(phy)
plot(ss, axes=FALSE, ylab="", xlab="", ylim=c(0,300), xlim=c(-120,0), col="white")
lines(ss)
par(tck=.01)
axis(2, at = seq(0, 300, by = 50), las=1, lwd=1, cex.axis=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(-120, 0, by = 20), las=1, lwd=1, cex.axis=1, labels=c(120,100,80,60,40,20,0), mgp=c(.75,.5,0))
title(ylab="N", line=2.5, cex=2)
title(xlab="Time", line=2.5)
mtext("a)",side=3, line=0, adj=0, cex=1.5)


######################################################################################################################################
######################################################################################################################################
### The boxplots -- PANELS B AND C
######################################################################################################################################
######################################################################################################################################

turnover <- c(0.3, 0.6, 1.5, 1.0)
eps <- c(0.7, 0.9, 0.95, 0.8)
transition <- 0.1*matrix(c(NA, 0.5, 0.2, 0.6, 0.1, NA, 0.3, 0.2, 0.1, 0.4, NA, 0.8, 0.1, 0.2, 0.3, NA),nrow=4)
mass.extinction.heights <- c(40,70)
mass.extinction.magnitudes <- c(0.7, 0.9)
cols <- magma(6)
col.vec <- c(cols[2], cols[3], cols[4], cols[5], 0)

GetActualDifference <- function(x) {
    x <- strsplit(x, " ")[[1]]
    return(as.numeric(x[1])-as.numeric(x[2]))
}


LoopActualDifference <- function(v) {
    v_result <- v
    for (i in seq_along(v)) {
        v_result[i] <- GetActualDifference(v[i])
    }
    return(v_result)
}


MatrixBreakdown <- function(compare_values){
    compare_values_mat <- matrix(compare_values, nrow=4, byrow=TRUE)
    compare.list <- as.list(1:2)
    compare.list[[1]] <- diag(compare_values_mat)
    mirror.comparisons <- c()
    for(i in 2:dim(compare_values_mat)[1]){
        for(j in 1:(i-1)){
            if(i != j){
                mirror.comparisons  <- rbind(mirror.comparisons , c(compare_values_mat[i,j], compare_values_mat[j,i]))
            }
        }
    }
    for(mirror.index in 1:dim(mirror.comparisons)[1]){
        compare.list[[mirror.index+1]] <- mirror.comparisons[mirror.index,]
    }
    return(compare.list)
}



GetSError <- function(data){
    interval <- seq(0,1, by=0.01)
    res <- numeric(length(interval))
    for(index in 1:length(interval)){
        tmp <- quantile(data, probs=interval[index])
        if(tmp < 0.5){
            res[index] <- 1
        }
    }
    s.val <- max(which(res==1))
    if(is.finite(s.val)){
        return(interval[s.val])
    }else{
        return(1)
    }
}


######################################################################################################################################
######################################################################################################################################
### Extant -- ROW 1
######################################################################################################################################
######################################################################################################################################

setwd(scen4.dir)
files.fullset <- list.files(pattern="Rsave")
true.psi <- "1"
keep <- c()
for(file.index in 1:length(files.fullset)){
    tmp <- strsplit(files.fullset[file.index], "[.]")
    if(tmp[[1]][4] == true.psi){
        keep <- c(keep, file.index)
    }
}
files <- files.fullset[keep]


all_results_full <- data.frame()
all_results <- data.frame()
all_freqs <- data.frame()
for(i in seq_along(files)) {
    try(rm(misse.extant.tip))
    try(rm(true.tip.rates))
    load(files[i])
    if(!is.null(misse.extant.tip)) {
        true.tip.rates.extant <- true.tip.rates[names(true.tip.rates) %in% misse.extant.tip$taxon]
        misse.extant.tip$turnover.true <- turnover[1+unname(true.tip.rates[misse.extant.tip$taxon])]
        misse.extant.tip$eps.true <- eps[1+unname(true.tip.rates[misse.extant.tip$taxon])]
        misse.extant.tip$state.true <- true.tip.rates[misse.extant.tip$taxon]
        misse.extant.tip$file.index <- i
        misse.extant.tip$file.name <- files[i]
        all_results <- rbind(all_results, misse.extant.tip)
        all_freqs <- rbind(all_freqs, table(true.tip.rates.extant))
    }
}
all_results$true_minus_est_eps <- all_results$eps.true-all_results$extinct.frac
all_results$true_minus_est_turnover <- all_results$turnover.true-all_results$turnover
all_results$misse <- "extant"
all_results_full <- all_results

f1 <- subset(all_results, file.index==1)
compare_values <- apply(expand.grid(unique(f1$turnover.true), unique(f1$turnover.true)), 1, paste, collapse=" ")
compare_values.list <- MatrixBreakdown(compare_values)

pairwise_summaries <- data.frame(matrix(nrow=0, ncol=8))
colnames(pairwise_summaries) <- c("file.index", "turnover.pairing", "median.turnover.actual", "fraction.turnover.difference.actual.positive", "number.comparisons", "median.turnover.difference", "ci.25.turnover.difference", "ci.75.turnover.difference")
for (index in sequence(max(all_results$file.index))) {
    print(paste(index, " of ", max(all_results$file.index)))
    f1 <- subset(all_results, file.index==index)
    distmat.turnover <- outer (f1$turnover, f1$turnover, `-`)
    distmat.turnover.compare_true <- outer (f1$turnover.true, f1$turnover.true, `paste`)
    delta.turnover.actual <- list()
    for (difference.index in sequence(length(compare_values.list))) {
        if(difference.index == 1){
            actual_turn <- c()
            for(compare.index in sequence(length(compare_values.list[[difference.index]]))){
                actual_turn <- c(actual_turn, as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])]))
            }
        }
        else{
            compare.index <- 2
            actual_turn <- as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])])
        }
        delta.turnover.actual[[difference.index]] <- actual_turn
        pairwise_summaries[nrow(pairwise_summaries)+1,] <- c(index, compare_values.list[[difference.index]][2], median(f1$turnover), sum(actual_turn>0)/length(actual_turn), length(actual_turn), median(actual_turn), quantile(actual_turn, 0.25), quantile(actual_turn, 0.75))
    }
}

pairwise_summaries$misse <- "extant"
pairwise_summaries$true_difference_turnover <- LoopActualDifference(pairwise_summaries$turnover.pairing)
pairwise_summaries$median.turnover.difference <- as.numeric(pairwise_summaries$median.turnover.difference)
pairwise_summaries$fraction.turnover.difference.actual.positive <- as.numeric(pairwise_summaries$fraction.turnover.difference.actual.positive)

#ss <- boxplot(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, plot=FALSE)
#print(ss$stat)
extant_01 <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="GetSError"), c("Diff", "<50%"))
extant_01_med <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="median"), c("Diff", "Median"))


pairwise_summaries_full <- pairwise_summaries


######################################################################################################################################
######################################################################################################################################
### Full K sampling
######################################################################################################################################
######################################################################################################################################

all_results <- data.frame()
all_freqs <- data.frame()
for(i in seq_along(files)) {
    try(rm(misse.mk.full.tip))
    try(rm(true.tip.rates))
    load(files[i])
    if(!is.null(misse.mk.full.tip)) {
        true.tip.rates.extant <- true.tip.rates[names(true.tip.rates) %in% misse.mk.full.tip$taxon]
        misse.mk.full.tip$turnover.true <- turnover[1+unname(true.tip.rates[misse.mk.full.tip$taxon])]
        misse.mk.full.tip$eps.true <- eps[1+unname(true.tip.rates[misse.mk.full.tip$taxon])]
        misse.mk.full.tip$state.true <- true.tip.rates[misse.mk.full.tip$taxon]
        misse.mk.full.tip$file.index <- i
        misse.mk.full.tip$file.name <- files[i]
        all_results <- rbind(all_results, misse.mk.full.tip[1:200,])
        all_freqs <- rbind(all_freqs, table(true.tip.rates.extant))
    }
}
all_results$true_minus_est_eps <- all_results$eps.true-all_results$extinct.frac
all_results$true_minus_est_turnover <- all_results$turnover.true-all_results$turnover
all_results$misse <- "full_k"
all_results_full <- rbind(all_results_full, all_results)

pairwise_summaries <- data.frame(matrix(nrow=0, ncol=8))
colnames(pairwise_summaries) <- c("file.index", "turnover.pairing", "median.turnover.actual", "fraction.turnover.difference.actual.positive", "number.comparisons", "median.turnover.difference", "ci.25.turnover.difference", "ci.75.turnover.difference")
for (index in sequence(max(all_results$file.index))) {
    print(paste(index, " of ", max(all_results$file.index)))
    f1 <- subset(all_results, file.index==index)
    distmat.turnover <- outer (f1$turnover, f1$turnover, `-`)
    distmat.turnover.compare_true <- outer (f1$turnover.true, f1$turnover.true, `paste`)
    delta.turnover.actual <- list()
    for (difference.index in sequence(length(compare_values.list))) {
        if(difference.index == 1){
            actual_turn <- c()
            for(compare.index in sequence(length(compare_values.list[[difference.index]]))){
                actual_turn <- c(actual_turn, as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])]))
            }
            compare_values_it <- compare_values.list[[difference.index]][1]
        }
        else{
            compare.index = 2
            compare_values_it <- compare_values.list[[difference.index]][1]
            actual_turn <- as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])])
        }
        delta.turnover.actual[[difference.index]] <- actual_turn
        pairwise_summaries[nrow(pairwise_summaries)+1,] <- c(index, compare_values.list[[difference.index]][2], median(f1$turnover), sum(actual_turn>0)/length(actual_turn), length(actual_turn), median(actual_turn), quantile(actual_turn, 0.25),quantile(actual_turn, 0.75))
    }
}

pairwise_summaries$misse <- "full_k"
pairwise_summaries$true_difference_turnover <- LoopActualDifference(pairwise_summaries$turnover.pairing)
pairwise_summaries$median.turnover.difference <- as.numeric(pairwise_summaries$median.turnover.difference)
pairwise_summaries$fraction.turnover.difference.actual.positive <- as.numeric(pairwise_summaries$fraction.turnover.difference.actual.positive)

#ss <- boxplot(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, plot=FALSE)
#print(ss$stat)

full_k_0.01 <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="GetSError"), c("Diff", "<50%"))

full_k_0.01_med <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="median"), c("Diff", "Median"))

pairwise_summaries_full <- rbind(pairwise_summaries_full, pairwise_summaries)


######################################################################################################################################
######################################################################################################################################
### Half K sampling
######################################################################################################################################
######################################################################################################################################

all_results <- data.frame()
all_freqs <- data.frame()
for(i in seq_along(files)) {
    try(rm(misse.mk.missinghalf.tip))
    try(rm(true.tip.rates))
    load(files[i])
    if(!is.null(misse.mk.missinghalf.tip)) {
        true.tip.rates.extant <- true.tip.rates[names(true.tip.rates) %in% misse.mk.missinghalf.tip$taxon]
        misse.mk.missinghalf.tip$turnover.true <- turnover[1+unname(true.tip.rates[misse.mk.missinghalf.tip$taxon])]
        misse.mk.missinghalf.tip$eps.true <- eps[1+unname(true.tip.rates[misse.mk.missinghalf.tip$taxon])]
        misse.mk.missinghalf.tip$state.true <- true.tip.rates[misse.mk.missinghalf.tip$taxon]
        misse.mk.missinghalf.tip$file.index <- i
        misse.mk.missinghalf.tip$file.name <- files[i]
        all_results <- rbind(all_results, misse.mk.missinghalf.tip[1:200,])
        all_freqs <- rbind(all_freqs, table(true.tip.rates.extant))
    }
}
all_results$true_minus_est_eps <- all_results$eps.true-all_results$extinct.frac
all_results$true_minus_est_turnover <- all_results$turnover.true-all_results$turnover
all_results$misse <- "half_k"
all_results_full <- rbind(all_results_full, all_results)

pairwise_summaries <- data.frame(matrix(nrow=0, ncol=8))
colnames(pairwise_summaries) <- c("file.index", "turnover.pairing", "median.turnover.actual", "fraction.turnover.difference.actual.positive", "number.comparisons", "median.turnover.difference", "ci.25.turnover.difference", "ci.75.turnover.difference")
for (index in sequence(max(all_results$file.index))) {
    print(paste(index, " of ", max(all_results$file.index)))
    f1 <- subset(all_results, file.index==index)
    distmat.turnover <- outer (f1$turnover, f1$turnover, `-`)
    distmat.turnover.compare_true <- outer (f1$turnover.true, f1$turnover.true, `paste`)
    delta.turnover.actual <- list()
    for (difference.index in sequence(length(compare_values.list))) {
        if(difference.index == 1){
            actual_turn <- c()
            for(compare.index in sequence(length(compare_values.list[[difference.index]]))){
                actual_turn <- c(actual_turn, as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])]))
            }
            compare_values_it <- compare_values.list[[difference.index]][1]
        }
        else{
            compare.index = 2
            compare_values_it <- compare_values.list[[difference.index]][1]
            actual_turn <- as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])])
        }
        delta.turnover.actual[[difference.index]] <- actual_turn
        pairwise_summaries[nrow(pairwise_summaries)+1,] <- c(index, compare_values.list[[difference.index]][2], median(f1$turnover), sum(actual_turn>0)/length(actual_turn), length(actual_turn), median(actual_turn), quantile(actual_turn, 0.25),quantile(actual_turn, 0.75))
    }
}

pairwise_summaries$misse <- "half_k"
pairwise_summaries$true_difference_turnover <- LoopActualDifference(pairwise_summaries$turnover.pairing)
pairwise_summaries$median.turnover.difference <- as.numeric(pairwise_summaries$median.turnover.difference)
pairwise_summaries$fraction.turnover.difference.actual.positive <- as.numeric(pairwise_summaries$fraction.turnover.difference.actual.positive)

#ss <- boxplot(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, plot=FALSE)
#print(ss$stat)
half_k_0.01 <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="GetSError"), c("Diff", "<50%"))
half_k_0.01_med <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="median"), c("Diff", "Median"))

pairwise_summaries_full <- rbind(pairwise_summaries_full, pairwise_summaries)

######################################################################################################################################
######################################################################################################################################
### No K sampling
######################################################################################################################################
######################################################################################################################################

all_results <- data.frame()
all_freqs <- data.frame()
for(i in seq_along(files)) {
    try(rm(misse.zero_k.tip))
    try(rm(true.tip.rates))
    load(files[i])
    if(!is.null(misse.zero_k.tip)) {
        true.tip.rates.extant <- true.tip.rates[names(true.tip.rates) %in% misse.zero_k.tip$taxon]
        misse.zero_k.tip$turnover.true <- turnover[1+unname(true.tip.rates[misse.zero_k.tip$taxon])]
        misse.zero_k.tip$eps.true <- eps[1+unname(true.tip.rates[misse.zero_k.tip$taxon])]
        misse.zero_k.tip$state.true <- true.tip.rates[misse.zero_k.tip$taxon]
        misse.zero_k.tip$file.index <- i
        misse.zero_k.tip$file.name <- files[i]
        all_results <- rbind(all_results, misse.zero_k.tip[1:200,])
        all_freqs <- rbind(all_freqs, table(true.tip.rates.extant))
    }
}
all_results$true_minus_est_eps <- all_results$eps.true-all_results$extinct.frac
all_results$true_minus_est_turnover <- all_results$turnover.true-all_results$turnover
all_results$misse <- "zero_k"
all_results_full <- rbind(all_results_full, all_results)


pairwise_summaries <- data.frame(matrix(nrow=0, ncol=8))
colnames(pairwise_summaries) <- c("file.index", "turnover.pairing", "median.turnover.actual", "fraction.turnover.difference.actual.positive", "number.comparisons", "median.turnover.difference", "ci.25.turnover.difference", "ci.75.turnover.difference")
for (index in sequence(max(all_results$file.index))) {
    print(paste(index, " of ", max(all_results$file.index)))
    f1 <- subset(all_results, file.index==index)
    distmat.turnover <- outer (f1$turnover, f1$turnover, `-`)
    distmat.turnover.compare_true <- outer (f1$turnover.true, f1$turnover.true, `paste`)
    delta.turnover.actual <- list()
    for (difference.index in sequence(length(compare_values.list))) {
        if(difference.index == 1){
            actual_turn <- c()
            for(compare.index in sequence(length(compare_values.list[[difference.index]]))){
                actual_turn <- c(actual_turn, as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])]))
            }
            compare_values_it <- compare_values.list[[difference.index]][1]
        }
        else{
            compare.index = 2
            compare_values_it <- compare_values.list[[difference.index]][1]
            actual_turn <- as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])])
        }
        delta.turnover.actual[[difference.index]] <- actual_turn
        pairwise_summaries[nrow(pairwise_summaries)+1,] <- c(index, compare_values.list[[difference.index]][2], median(f1$turnover), sum(actual_turn>0)/length(actual_turn), length(actual_turn), median(actual_turn), quantile(actual_turn, 0.25),quantile(actual_turn, 0.75))
    }
}

pairwise_summaries$misse <- "zero_k"
pairwise_summaries$true_difference_turnover <- LoopActualDifference(pairwise_summaries$turnover.pairing)
pairwise_summaries$median.turnover.difference <- as.numeric(pairwise_summaries$median.turnover.difference)
pairwise_summaries$fraction.turnover.difference.actual.positive <- as.numeric(pairwise_summaries$fraction.turnover.difference.actual.positive)

#ss <- boxplot(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, plot=FALSE)
#print(ss$stat)
zero_k_0.01 <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="GetSError"), c("Diff", "<50%"))
zero_k_0.01_med <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="median"), c("Diff", "Median"))


pairwise_summaries_full <- rbind(pairwise_summaries_full, pairwise_summaries)


######################################################################################################################################
######################################################################################################################################
### Strat Intervals -- ROW 1
######################################################################################################################################
######################################################################################################################################

all_results <- data.frame()
all_freqs <- data.frame()
for(i in seq_along(files)) {
    try(rm(misse.mk.strat.tip))
    try(rm(true.tip.rates))
    load(files[i])
    if(!is.null(misse.mk.strat.tip)) {
        true.tip.rates.extant <- true.tip.rates[names(true.tip.rates) %in% misse.mk.strat.tip$taxon]
        misse.mk.strat.tip$turnover.true <- turnover[1+unname(true.tip.rates[misse.mk.strat.tip$taxon])]
        misse.mk.strat.tip$eps.true <- eps[1+unname(true.tip.rates[misse.mk.strat.tip$taxon])]
        misse.mk.strat.tip$state.true <- true.tip.rates[misse.mk.strat.tip$taxon]
        misse.mk.strat.tip$file.index <- i
        misse.mk.strat.tip$file.name <- files[i]
        all_results <- rbind(all_results, misse.mk.strat.tip[1:200,])
        all_freqs <- rbind(all_freqs, table(true.tip.rates.extant))
    }
}
all_results$true_minus_est_eps <- all_results$eps.true-all_results$extinct.frac
all_results$true_minus_est_turnover <- all_results$turnover.true-all_results$turnover
all_results$misse <- "fstrat_int"
all_results_full <- rbind(all_results_full, all_results)


pairwise_summaries <- data.frame(matrix(nrow=0, ncol=8))
colnames(pairwise_summaries) <- c("file.index", "turnover.pairing", "median.turnover.actual", "fraction.turnover.difference.actual.positive", "number.comparisons", "median.turnover.difference", "ci.25.turnover.difference", "ci.75.turnover.difference")
for (index in sequence(max(all_results$file.index))) {
    print(paste(index, " of ", max(all_results$file.index)))
    f1 <- subset(all_results, file.index==index)
    distmat.turnover <- outer (f1$turnover, f1$turnover, `-`)
    distmat.turnover.compare_true <- outer (f1$turnover.true, f1$turnover.true, `paste`)
    delta.turnover.actual <- list()
    for (difference.index in sequence(length(compare_values.list))) {
        if(difference.index == 1){
            actual_turn <- c()
            for(compare.index in sequence(length(compare_values.list[[difference.index]]))){
                actual_turn <- c(actual_turn, as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])]))
            }
            compare_values_it <- compare_values.list[[difference.index]][1]
        }
        else{
            compare.index = 2
            compare_values_it <- compare_values.list[[difference.index]][1]
            actual_turn <- as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])])
        }
        delta.turnover.actual[[difference.index]] <- actual_turn
        pairwise_summaries[nrow(pairwise_summaries)+1,] <- c(index, compare_values.list[[difference.index]][2], median(f1$turnover), sum(actual_turn>0)/length(actual_turn), length(actual_turn), median(actual_turn), quantile(actual_turn, 0.25),quantile(actual_turn, 0.75))
    }
}

pairwise_summaries$misse <- "fstrat_int"
pairwise_summaries$true_difference_turnover <- LoopActualDifference(pairwise_summaries$turnover.pairing)
pairwise_summaries$median.turnover.difference <- as.numeric(pairwise_summaries$median.turnover.difference)
pairwise_summaries$fraction.turnover.difference.actual.positive <- as.numeric(pairwise_summaries$fraction.turnover.difference.actual.positive)

#ss <- boxplot(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, plot=FALSE)
#print(ss$stat)

strat_int_0.1 <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="GetSError"), c("Diff", "<50%"))
strat_int_0.1_med <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="median"), c("Diff", "Median"))

pairwise_summaries_full <- rbind(pairwise_summaries_full, pairwise_summaries)


boxplot(all_results_full$turnover~all_results_full$misse*all_results_full$turnover.true, col=col.vec, axes=FALSE, xlab="", ylim=c(0,2.5), ylab="", main=expression(psi==0.01), at=c(1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17, 19,20,21,22,23), outline=FALSE)
segments(0.5, 0.3, 5.5, 0.3, col="red", lwd=1.5)
segments(6.5, 0.6, 11.5, 0.6, col="red", lwd=1.5)
segments(12.5, 1.0, 17.5, 1.0, col="red", lwd=1.5)
segments(18.5, 1.5, 23.5, 1.5, col="red", lwd=1.5)
par(tck=.01)
axis(2, at = seq(0, 2.5, by = .5), las=1, lwd=1, cex.axis=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = c(3, 9, 15, 21), las=1, lwd=1, cex.axis=1, labels=c(0.3, 0.6, 1.0, 1.5), mgp=c(.75,.5,0))
title(ylab=expression(Estimated~hat(tau)), line=2.5, cex=2)
title(xlab=expression(True~tau), line=2.5)
legend("topleft", c("Extant only", "Stratigraphic ranges", "All fossils used", expression(Half~of~italic(k)~fossils~used), expression(Only~italic(m)~fossils~used)), fill=col.vec, box.col=0)
mtext("b)",side=3, line=0, adj=0, cex=1.5)




#boxplot(pairwise_summaries_full$median.turnover.difference~pairwise_summaries_full$misse*pairwise_summaries_full$true_difference_turnover, col=viridis(4), axes=FALSE, ylim=c(-1,4), xlab="", ylab="", main="", outline=FALSE)
#par(tck=.01)
#axis(2, at = seq(-1, 4, by = 1), las=1, lwd=1, cex.axis=1, labels=TRUE, mgp=c(.75,.5,0))
#axis(1, at = c(2.5, 6.5, 10.5, 14.5, 18.5, 22.5, 26.5), las=1, lwd=1, cex.axis=1, labels=c(-0.5,0.3,0.4,0.5,0.7,0.9,1.2), mgp=c(.75,.5,0))
#title(ylab=expression(Estimated~hat(tau[i])-hat(tau[j])), line=2.5, cex=2)
#title(xlab=expression(True~tau[i]-tau[j]), line=2.5)
#mtext("D",side=3, line=0, adj=0, cex=1.5)


######################################################################################################################################
######################################################################################################################################
### Extant -- ROW 3
######################################################################################################################################
######################################################################################################################################

files.fullset <- list.files(pattern="Rsave")
true.psi <- "01"
keep <- c()
for(file.index in 1:length(files.fullset)){
    tmp <- strsplit(files.fullset[file.index], "[.]")
    if(tmp[[1]][4] == true.psi){
        keep <- c(keep, file.index)
    }
}
files <- files.fullset[keep]

all_results_full <- data.frame()
all_results <- data.frame()
all_freqs <- data.frame()
for(i in seq_along(files)) {
    try(rm(misse.extant.tip))
    try(rm(true.tip.rates))
    load(files[i])
    if(!is.null(misse.extant.tip)) {
        true.tip.rates.extant <- true.tip.rates[names(true.tip.rates) %in% misse.extant.tip$taxon]
        misse.extant.tip$turnover.true <- turnover[1+unname(true.tip.rates[misse.extant.tip$taxon])]
        misse.extant.tip$eps.true <- eps[1+unname(true.tip.rates[misse.extant.tip$taxon])]
        misse.extant.tip$state.true <- true.tip.rates[misse.extant.tip$taxon]
        misse.extant.tip$file.index <- i
        misse.extant.tip$file.name <- files[i]
        all_results <- rbind(all_results, misse.extant.tip)
        all_freqs <- rbind(all_freqs, table(true.tip.rates.extant))
    }
}
all_results$true_minus_est_eps <- all_results$eps.true-all_results$extinct.frac
all_results$true_minus_est_turnover <- all_results$turnover.true-all_results$turnover
all_results$misse <- "extant"
all_results_full <- all_results

f1 <- subset(all_results, file.index==1)
compare_values <- apply(expand.grid(unique(f1$turnover.true), unique(f1$turnover.true)), 1, paste, collapse=" ")
compare_values.list <- MatrixBreakdown(compare_values)

pairwise_summaries <- data.frame(matrix(nrow=0, ncol=8))
colnames(pairwise_summaries) <- c("file.index", "turnover.pairing", "median.turnover.actual", "fraction.turnover.difference.actual.positive", "number.comparisons", "median.turnover.difference", "ci.25.turnover.difference", "ci.75.turnover.difference")
for (index in sequence(max(all_results$file.index))) {
    print(paste(index, " of ", max(all_results$file.index)))
    f1 <- subset(all_results, file.index==index)
    distmat.turnover <- outer (f1$turnover, f1$turnover, `-`)
    distmat.turnover.compare_true <- outer (f1$turnover.true, f1$turnover.true, `paste`)
    delta.turnover.actual <- list()
    for (difference.index in sequence(length(compare_values.list))) {
        if(difference.index == 1){
            actual_turn <- c()
            for(compare.index in sequence(length(compare_values.list[[difference.index]]))){
                actual_turn <- c(actual_turn, as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])]))
            }
        }
        else{
            compare.index <- 2
            actual_turn <- as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])])
        }
        delta.turnover.actual[[difference.index]] <- actual_turn
        pairwise_summaries[nrow(pairwise_summaries)+1,] <- c(index, compare_values.list[[difference.index]][2], median(f1$turnover), sum(actual_turn>0)/length(actual_turn), length(actual_turn), median(actual_turn), quantile(actual_turn, 0.25), quantile(actual_turn, 0.75))
    }
}

pairwise_summaries$misse <- "extant"
pairwise_summaries$true_difference_turnover <- LoopActualDifference(pairwise_summaries$turnover.pairing)
pairwise_summaries$median.turnover.difference <- as.numeric(pairwise_summaries$median.turnover.difference)
pairwise_summaries$fraction.turnover.difference.actual.positive <- as.numeric(pairwise_summaries$fraction.turnover.difference.actual.positive)

#ss <- boxplot(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, plot=FALSE)
#print("extant_0.01")
#print(ss$stat)

pairwise_summaries_full <- pairwise_summaries


######################################################################################################################################
######################################################################################################################################
### Full K sampling -- ROW 3
######################################################################################################################################
######################################################################################################################################
files.fullset <- list.files(pattern="Rsave")
true.psi <- "1"
keep <- c()
for(file.index in 1:length(files.fullset)){
    tmp <- strsplit(files.fullset[file.index], "[.]")
    if(tmp[[1]][4] == true.psi){
        keep <- c(keep, file.index)
    }
}
files <- files.fullset[keep]

all_results <- data.frame()
all_freqs <- data.frame()
for(i in seq_along(files)) {
    try(rm(misse.mk.full.tip))
    try(rm(true.tip.rates))
    load(files[i])
    if(!is.null(misse.mk.full.tip)) {
        true.tip.rates.extant <- true.tip.rates[names(true.tip.rates) %in% misse.mk.full.tip$taxon]
        misse.mk.full.tip$turnover.true <- turnover[1+unname(true.tip.rates[misse.mk.full.tip$taxon])]
        misse.mk.full.tip$eps.true <- eps[1+unname(true.tip.rates[misse.mk.full.tip$taxon])]
        misse.mk.full.tip$state.true <- true.tip.rates[misse.mk.full.tip$taxon]
        misse.mk.full.tip$file.index <- i
        misse.mk.full.tip$file.name <- files[i]
        all_results <- rbind(all_results, misse.mk.full.tip[1:200,])
        all_freqs <- rbind(all_freqs, table(true.tip.rates.extant))
    }
}
all_results$true_minus_est_eps <- all_results$eps.true-all_results$extinct.frac
all_results$true_minus_est_turnover <- all_results$turnover.true-all_results$turnover
all_results$misse <- "full_k"
all_results_full <- rbind(all_results_full, all_results)

pairwise_summaries <- data.frame(matrix(nrow=0, ncol=8))
colnames(pairwise_summaries) <- c("file.index", "turnover.pairing", "median.turnover.actual", "fraction.turnover.difference.actual.positive", "number.comparisons", "median.turnover.difference", "ci.25.turnover.difference", "ci.75.turnover.difference")
for (index in sequence(max(all_results$file.index))) {
    print(paste(index, " of ", max(all_results$file.index)))
    f1 <- subset(all_results, file.index==index)
    distmat.turnover <- outer (f1$turnover, f1$turnover, `-`)
    distmat.turnover.compare_true <- outer (f1$turnover.true, f1$turnover.true, `paste`)
    delta.turnover.actual <- list()
    for (difference.index in sequence(length(compare_values.list))) {
        if(difference.index == 1){
            actual_turn <- c()
            for(compare.index in sequence(length(compare_values.list[[difference.index]]))){
                actual_turn <- c(actual_turn, as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])]))
            }
            compare_values_it <- compare_values.list[[difference.index]][1]
        }
        else{
            compare.index = 2
            compare_values_it <- compare_values.list[[difference.index]][1]
            actual_turn <- as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])])
        }
        delta.turnover.actual[[difference.index]] <- actual_turn
        pairwise_summaries[nrow(pairwise_summaries)+1,] <- c(index, compare_values.list[[difference.index]][2], median(f1$turnover), sum(actual_turn>0)/length(actual_turn), length(actual_turn), median(actual_turn), quantile(actual_turn, 0.25),quantile(actual_turn, 0.75))
    }
}

pairwise_summaries$misse <- "full_k"
pairwise_summaries$true_difference_turnover <- LoopActualDifference(pairwise_summaries$turnover.pairing)
pairwise_summaries$median.turnover.difference <- as.numeric(pairwise_summaries$median.turnover.difference)
pairwise_summaries$fraction.turnover.difference.actual.positive <- as.numeric(pairwise_summaries$fraction.turnover.difference.actual.positive)

#ss <- boxplot(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, plot=FALSE)
#print(ss$stat)

full_k_0.1 <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="GetSError"), c("Diff", "<50%"))
full_k_0.1_med <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="median"), c("Diff", "Median"))


pairwise_summaries_full <- rbind(pairwise_summaries_full, pairwise_summaries)


######################################################################################################################################
######################################################################################################################################
### Half K sampling -- ROW 3
######################################################################################################################################
######################################################################################################################################

all_results <- data.frame()
all_freqs <- data.frame()
for(i in seq_along(files)) {
    try(rm(misse.mk.missinghalf.tip))
    try(rm(true.tip.rates))
    load(files[i])
    if(!is.null(misse.mk.missinghalf.tip)) {
        true.tip.rates.extant <- true.tip.rates[names(true.tip.rates) %in% misse.mk.missinghalf.tip$taxon]
        misse.mk.missinghalf.tip$turnover.true <- turnover[1+unname(true.tip.rates[misse.mk.missinghalf.tip$taxon])]
        misse.mk.missinghalf.tip$eps.true <- eps[1+unname(true.tip.rates[misse.mk.missinghalf.tip$taxon])]
        misse.mk.missinghalf.tip$state.true <- true.tip.rates[misse.mk.missinghalf.tip$taxon]
        misse.mk.missinghalf.tip$file.index <- i
        misse.mk.missinghalf.tip$file.name <- files[i]
        all_results <- rbind(all_results, misse.mk.missinghalf.tip[1:200,])
        all_freqs <- rbind(all_freqs, table(true.tip.rates.extant))
    }
}
all_results$true_minus_est_eps <- all_results$eps.true-all_results$extinct.frac
all_results$true_minus_est_turnover <- all_results$turnover.true-all_results$turnover
all_results$misse <- "half_k"
all_results_full <- rbind(all_results_full, all_results)

pairwise_summaries <- data.frame(matrix(nrow=0, ncol=8))
colnames(pairwise_summaries) <- c("file.index", "turnover.pairing", "median.turnover.actual", "fraction.turnover.difference.actual.positive", "number.comparisons", "median.turnover.difference", "ci.25.turnover.difference", "ci.75.turnover.difference")
for (index in sequence(max(all_results$file.index))) {
    print(paste(index, " of ", max(all_results$file.index)))
    f1 <- subset(all_results, file.index==index)
    distmat.turnover <- outer (f1$turnover, f1$turnover, `-`)
    distmat.turnover.compare_true <- outer (f1$turnover.true, f1$turnover.true, `paste`)
    delta.turnover.actual <- list()
    for (difference.index in sequence(length(compare_values.list))) {
        if(difference.index == 1){
            actual_turn <- c()
            for(compare.index in sequence(length(compare_values.list[[difference.index]]))){
                actual_turn <- c(actual_turn, as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])]))
            }
            compare_values_it <- compare_values.list[[difference.index]][1]
        }
        else{
            compare.index = 2
            compare_values_it <- compare_values.list[[difference.index]][1]
            actual_turn <- as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])])
        }
        delta.turnover.actual[[difference.index]] <- actual_turn
        pairwise_summaries[nrow(pairwise_summaries)+1,] <- c(index, compare_values.list[[difference.index]][2], median(f1$turnover), sum(actual_turn>0)/length(actual_turn), length(actual_turn), median(actual_turn), quantile(actual_turn, 0.25),quantile(actual_turn, 0.75))
    }
}

pairwise_summaries$misse <- "half_k"
pairwise_summaries$true_difference_turnover <- LoopActualDifference(pairwise_summaries$turnover.pairing)
pairwise_summaries$median.turnover.difference <- as.numeric(pairwise_summaries$median.turnover.difference)
pairwise_summaries$fraction.turnover.difference.actual.positive <- as.numeric(pairwise_summaries$fraction.turnover.difference.actual.positive)

#ss <- boxplot(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, plot=FALSE)

#print(ss$stat)
half_k_0.1 <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="GetSError"), c("Diff", "<50%"))
half_k_0.1_med <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="median"), c("Diff", "Median"))

pairwise_summaries_full <- rbind(pairwise_summaries_full, pairwise_summaries)

######################################################################################################################################
######################################################################################################################################
### No K sampling -- ROW 3
######################################################################################################################################
######################################################################################################################################

all_results <- data.frame()
all_freqs <- data.frame()
for(i in seq_along(files)) {
    try(rm(misse.zero_k.tip))
    try(rm(true.tip.rates))
    load(files[i])
    if(!is.null(misse.zero_k.tip)) {
        true.tip.rates.extant <- true.tip.rates[names(true.tip.rates) %in% misse.zero_k.tip$taxon]
        misse.zero_k.tip$turnover.true <- turnover[1+unname(true.tip.rates[misse.zero_k.tip$taxon])]
        misse.zero_k.tip$eps.true <- eps[1+unname(true.tip.rates[misse.zero_k.tip$taxon])]
        misse.zero_k.tip$state.true <- true.tip.rates[misse.zero_k.tip$taxon]
        misse.zero_k.tip$file.index <- i
        misse.zero_k.tip$file.name <- files[i]
        all_results <- rbind(all_results, misse.zero_k.tip[1:200,])
        all_freqs <- rbind(all_freqs, table(true.tip.rates.extant))
    }
}
all_results$true_minus_est_eps <- all_results$eps.true-all_results$extinct.frac
all_results$true_minus_est_turnover <- all_results$turnover.true-all_results$turnover
all_results$misse <- "zero_k"
all_results_full <- rbind(all_results_full, all_results)


pairwise_summaries <- data.frame(matrix(nrow=0, ncol=8))
colnames(pairwise_summaries) <- c("file.index", "turnover.pairing", "median.turnover.actual", "fraction.turnover.difference.actual.positive", "number.comparisons", "median.turnover.difference", "ci.25.turnover.difference", "ci.75.turnover.difference")
for (index in sequence(max(all_results$file.index))) {
    print(paste(index, " of ", max(all_results$file.index)))
    f1 <- subset(all_results, file.index==index)
    distmat.turnover <- outer (f1$turnover, f1$turnover, `-`)
    distmat.turnover.compare_true <- outer (f1$turnover.true, f1$turnover.true, `paste`)
    delta.turnover.actual <- list()
    for (difference.index in sequence(length(compare_values.list))) {
        if(difference.index == 1){
            actual_turn <- c()
            for(compare.index in sequence(length(compare_values.list[[difference.index]]))){
                actual_turn <- c(actual_turn, as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])]))
            }
            compare_values_it <- compare_values.list[[difference.index]][1]
        }
        else{
            compare.index = 2
            compare_values_it <- compare_values.list[[difference.index]][1]
            actual_turn <- as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])])
        }
        delta.turnover.actual[[difference.index]] <- actual_turn
        pairwise_summaries[nrow(pairwise_summaries)+1,] <- c(index, compare_values.list[[difference.index]][2], median(f1$turnover), sum(actual_turn>0)/length(actual_turn), length(actual_turn), median(actual_turn), quantile(actual_turn, 0.25),quantile(actual_turn, 0.75))
    }
}

pairwise_summaries$misse <- "zero_k"
pairwise_summaries$true_difference_turnover <- LoopActualDifference(pairwise_summaries$turnover.pairing)
pairwise_summaries$median.turnover.difference <- as.numeric(pairwise_summaries$median.turnover.difference)
pairwise_summaries$fraction.turnover.difference.actual.positive <- as.numeric(pairwise_summaries$fraction.turnover.difference.actual.positive)

#ss <- boxplot(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, plot=FALSE)
#print(ss$stat)

zero_k_0.1 <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="GetSError"), c("Diff", "<50%"))
zero_k_0.1_med <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="median"), c("Diff", "Median"))

pairwise_summaries_full <- rbind(pairwise_summaries_full, pairwise_summaries)


######################################################################################################################################
######################################################################################################################################
### Strat Intervals -- ROW 3
######################################################################################################################################
######################################################################################################################################

all_results <- data.frame()
all_freqs <- data.frame()
for(i in seq_along(files)) {
    try(rm(misse.mk.strat.tip))
    try(rm(true.tip.rates))
    load(files[i])
    if(!is.null(misse.mk.strat.tip)) {
        true.tip.rates.extant <- true.tip.rates[names(true.tip.rates) %in% misse.mk.strat.tip$taxon]
        misse.mk.strat.tip$turnover.true <- turnover[1+unname(true.tip.rates[misse.mk.strat.tip$taxon])]
        misse.mk.strat.tip$eps.true <- eps[1+unname(true.tip.rates[misse.mk.strat.tip$taxon])]
        misse.mk.strat.tip$state.true <- true.tip.rates[misse.mk.strat.tip$taxon]
        misse.mk.strat.tip$file.index <- i
        misse.mk.strat.tip$file.name <- files[i]
        all_results <- rbind(all_results, misse.mk.strat.tip[1:200,])
        all_freqs <- rbind(all_freqs, table(true.tip.rates.extant))
    }
}
all_results$true_minus_est_eps <- all_results$eps.true-all_results$extinct.frac
all_results$true_minus_est_turnover <- all_results$turnover.true-all_results$turnover
all_results$misse <- "fstrat_int"
all_results_full <- rbind(all_results_full, all_results)


pairwise_summaries <- data.frame(matrix(nrow=0, ncol=8))
colnames(pairwise_summaries) <- c("file.index", "turnover.pairing", "median.turnover.actual", "fraction.turnover.difference.actual.positive", "number.comparisons", "median.turnover.difference", "ci.25.turnover.difference", "ci.75.turnover.difference")
for (index in sequence(max(all_results$file.index))) {
    print(paste(index, " of ", max(all_results$file.index)))
    f1 <- subset(all_results, file.index==index)
    distmat.turnover <- outer (f1$turnover, f1$turnover, `-`)
    distmat.turnover.compare_true <- outer (f1$turnover.true, f1$turnover.true, `paste`)
    delta.turnover.actual <- list()
    for (difference.index in sequence(length(compare_values.list))) {
        if(difference.index == 1){
            actual_turn <- c()
            for(compare.index in sequence(length(compare_values.list[[difference.index]]))){
                actual_turn <- c(actual_turn, as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])]))
            }
            compare_values_it <- compare_values.list[[difference.index]][1]
        }
        else{
            compare.index = 2
            compare_values_it <- compare_values.list[[difference.index]][1]
            actual_turn <- as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])])
        }
        delta.turnover.actual[[difference.index]] <- actual_turn
        pairwise_summaries[nrow(pairwise_summaries)+1,] <- c(index, compare_values.list[[difference.index]][2], median(f1$turnover), sum(actual_turn>0)/length(actual_turn), length(actual_turn), median(actual_turn), quantile(actual_turn, 0.25),quantile(actual_turn, 0.75))
    }
}

pairwise_summaries$misse <- "fstrat_int"
pairwise_summaries$true_difference_turnover <- LoopActualDifference(pairwise_summaries$turnover.pairing)
pairwise_summaries$median.turnover.difference <- as.numeric(pairwise_summaries$median.turnover.difference)
pairwise_summaries$fraction.turnover.difference.actual.positive <- as.numeric(pairwise_summaries$fraction.turnover.difference.actual.positive)

#ss <- boxplot(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, plot=FALSE)
#print(ss$stat)

strat_int_0.1 <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="GetSError"), c("Diff", "<50%"))
strat_int_0.1_med <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="median"), c("Diff", "Median"))

pairwise_summaries_full <- rbind(pairwise_summaries_full, pairwise_summaries)

boxplot(all_results_full$turnover~all_results_full$misse*all_results_full$turnover.true, col=col.vec, axes=FALSE, xlab="", ylim=c(0,2.5), ylab="", main=expression(psi==0.10), at=c(1,2,3,4,5, 7,8,9,10,11, 13,14,15,16,17, 19,20,21,22,23), outline=FALSE)
segments(0.5, 0.3, 5.5, 0.3, col="red", lwd=1.5)
segments(6.5, 0.6, 11.5, 0.6, col="red", lwd=1.5)
segments(12.5, 1.0, 17.5, 1.0, col="red", lwd=1.5)
segments(18.5, 1.5, 23.5, 1.5, col="red", lwd=1.5)
par(tck=.01)
axis(2, at = seq(0, 2.5, by = .5), las=1, lwd=1, cex.axis=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = c(3, 9, 15, 21), las=1, lwd=1, cex.axis=1, labels=c(0.3, 0.6, 1.0, 1.5), mgp=c(.75,.5,0))
title(ylab=expression(Estimated~hat(tau)), line=2.5, cex=2)
title(xlab=expression(True~tau), line=2.5)
mtext("c)",side=3, line=0, adj=0, cex=1.5)

dev.off()






######################################################################################################################################
######################################################################################################################################
### Full K sampling -- Psi = 0.05. Not shown in plot
######################################################################################################################################
######################################################################################################################################
files.fullset <- list.files(pattern="Rsave")
true.psi <- "05"
keep <- c()
for(file.index in 1:length(files.fullset)){
    tmp <- strsplit(files.fullset[file.index], "[.]")
    if(tmp[[1]][4] == true.psi){
        keep <- c(keep, file.index)
    }
}
files <- files.fullset[keep]

all_results <- data.frame()
all_freqs <- data.frame()
for(i in seq_along(files)) {
    try(rm(misse.mk.full.tip))
    try(rm(true.tip.rates))
    load(files[i])
    if(!is.null(misse.mk.full.tip)) {
        true.tip.rates.extant <- true.tip.rates[names(true.tip.rates) %in% misse.mk.full.tip$taxon]
        misse.mk.full.tip$turnover.true <- turnover[1+unname(true.tip.rates[misse.mk.full.tip$taxon])]
        misse.mk.full.tip$eps.true <- eps[1+unname(true.tip.rates[misse.mk.full.tip$taxon])]
        misse.mk.full.tip$state.true <- true.tip.rates[misse.mk.full.tip$taxon]
        misse.mk.full.tip$file.index <- i
        misse.mk.full.tip$file.name <- files[i]
        all_results <- rbind(all_results, misse.mk.full.tip[1:200,])
        all_freqs <- rbind(all_freqs, table(true.tip.rates.extant))
    }
}
all_results$true_minus_est_eps <- all_results$eps.true-all_results$extinct.frac
all_results$true_minus_est_turnover <- all_results$turnover.true-all_results$turnover
all_results$misse <- "full_k"
all_results_full <- rbind(all_results_full, all_results)

pairwise_summaries <- data.frame(matrix(nrow=0, ncol=8))
colnames(pairwise_summaries) <- c("file.index", "turnover.pairing", "median.turnover.actual", "fraction.turnover.difference.actual.positive", "number.comparisons", "median.turnover.difference", "ci.25.turnover.difference", "ci.75.turnover.difference")
for (index in sequence(max(all_results$file.index))) {
    print(paste(index, " of ", max(all_results$file.index)))
    f1 <- subset(all_results, file.index==index)
    distmat.turnover <- outer (f1$turnover, f1$turnover, `-`)
    distmat.turnover.compare_true <- outer (f1$turnover.true, f1$turnover.true, `paste`)
    delta.turnover.actual <- list()
    for (difference.index in sequence(length(compare_values.list))) {
        if(difference.index == 1){
            actual_turn <- c()
            for(compare.index in sequence(length(compare_values.list[[difference.index]]))){
                actual_turn <- c(actual_turn, as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])]))
            }
            compare_values_it <- compare_values.list[[difference.index]][1]
        }
        else{
            compare.index = 2
            compare_values_it <- compare_values.list[[difference.index]][1]
            actual_turn <- as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])])
        }
        delta.turnover.actual[[difference.index]] <- actual_turn
        pairwise_summaries[nrow(pairwise_summaries)+1,] <- c(index, compare_values.list[[difference.index]][2], median(f1$turnover), sum(actual_turn>0)/length(actual_turn), length(actual_turn), median(actual_turn), quantile(actual_turn, 0.25),quantile(actual_turn, 0.75))
    }
}

pairwise_summaries$misse <- "full_k"
pairwise_summaries$true_difference_turnover <- LoopActualDifference(pairwise_summaries$turnover.pairing)
pairwise_summaries$median.turnover.difference <- as.numeric(pairwise_summaries$median.turnover.difference)
pairwise_summaries$fraction.turnover.difference.actual.positive <- as.numeric(pairwise_summaries$fraction.turnover.difference.actual.positive)

#ss <- boxplot(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, plot=FALSE)
#print(ss$stat)

full_k_0.05 <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="GetSError"), c("Diff", "<50%"))
full_k_0.05_med <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="median"), c("Diff", "Median"))


pairwise_summaries_full <- rbind(pairwise_summaries_full, pairwise_summaries)


######################################################################################################################################
######################################################################################################################################
### Half K sampling -- ROW 3
######################################################################################################################################
######################################################################################################################################

all_results <- data.frame()
all_freqs <- data.frame()
for(i in seq_along(files)) {
    try(rm(misse.mk.missinghalf.tip))
    try(rm(true.tip.rates))
    load(files[i])
    if(!is.null(misse.mk.missinghalf.tip)) {
        true.tip.rates.extant <- true.tip.rates[names(true.tip.rates) %in% misse.mk.missinghalf.tip$taxon]
        misse.mk.missinghalf.tip$turnover.true <- turnover[1+unname(true.tip.rates[misse.mk.missinghalf.tip$taxon])]
        misse.mk.missinghalf.tip$eps.true <- eps[1+unname(true.tip.rates[misse.mk.missinghalf.tip$taxon])]
        misse.mk.missinghalf.tip$state.true <- true.tip.rates[misse.mk.missinghalf.tip$taxon]
        misse.mk.missinghalf.tip$file.index <- i
        misse.mk.missinghalf.tip$file.name <- files[i]
        all_results <- rbind(all_results, misse.mk.missinghalf.tip[1:200,])
        all_freqs <- rbind(all_freqs, table(true.tip.rates.extant))
    }
}
all_results$true_minus_est_eps <- all_results$eps.true-all_results$extinct.frac
all_results$true_minus_est_turnover <- all_results$turnover.true-all_results$turnover
all_results$misse <- "half_k"
all_results_full <- rbind(all_results_full, all_results)

pairwise_summaries <- data.frame(matrix(nrow=0, ncol=8))
colnames(pairwise_summaries) <- c("file.index", "turnover.pairing", "median.turnover.actual", "fraction.turnover.difference.actual.positive", "number.comparisons", "median.turnover.difference", "ci.25.turnover.difference", "ci.75.turnover.difference")
for (index in sequence(max(all_results$file.index))) {
    print(paste(index, " of ", max(all_results$file.index)))
    f1 <- subset(all_results, file.index==index)
    distmat.turnover <- outer (f1$turnover, f1$turnover, `-`)
    distmat.turnover.compare_true <- outer (f1$turnover.true, f1$turnover.true, `paste`)
    delta.turnover.actual <- list()
    for (difference.index in sequence(length(compare_values.list))) {
        if(difference.index == 1){
            actual_turn <- c()
            for(compare.index in sequence(length(compare_values.list[[difference.index]]))){
                actual_turn <- c(actual_turn, as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])]))
            }
            compare_values_it <- compare_values.list[[difference.index]][1]
        }
        else{
            compare.index = 2
            compare_values_it <- compare_values.list[[difference.index]][1]
            actual_turn <- as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])])
        }
        delta.turnover.actual[[difference.index]] <- actual_turn
        pairwise_summaries[nrow(pairwise_summaries)+1,] <- c(index, compare_values.list[[difference.index]][2], median(f1$turnover), sum(actual_turn>0)/length(actual_turn), length(actual_turn), median(actual_turn), quantile(actual_turn, 0.25),quantile(actual_turn, 0.75))
    }
}

pairwise_summaries$misse <- "half_k"
pairwise_summaries$true_difference_turnover <- LoopActualDifference(pairwise_summaries$turnover.pairing)
pairwise_summaries$median.turnover.difference <- as.numeric(pairwise_summaries$median.turnover.difference)
pairwise_summaries$fraction.turnover.difference.actual.positive <- as.numeric(pairwise_summaries$fraction.turnover.difference.actual.positive)

#ss <- boxplot(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, plot=FALSE)

#print(ss$stat)
half_k_0.05 <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="GetSError"), c("Diff", "<50%"))
half_k_0.05_med <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="median"), c("Diff", "Median"))

pairwise_summaries_full <- rbind(pairwise_summaries_full, pairwise_summaries)

######################################################################################################################################
######################################################################################################################################
### No K sampling -- ROW 3
######################################################################################################################################
######################################################################################################################################

all_results <- data.frame()
all_freqs <- data.frame()
for(i in seq_along(files)) {
    try(rm(misse.zero_k.tip))
    try(rm(true.tip.rates))
    load(files[i])
    if(!is.null(misse.zero_k.tip)) {
        true.tip.rates.extant <- true.tip.rates[names(true.tip.rates) %in% misse.zero_k.tip$taxon]
        misse.zero_k.tip$turnover.true <- turnover[1+unname(true.tip.rates[misse.zero_k.tip$taxon])]
        misse.zero_k.tip$eps.true <- eps[1+unname(true.tip.rates[misse.zero_k.tip$taxon])]
        misse.zero_k.tip$state.true <- true.tip.rates[misse.zero_k.tip$taxon]
        misse.zero_k.tip$file.index <- i
        misse.zero_k.tip$file.name <- files[i]
        all_results <- rbind(all_results, misse.zero_k.tip[1:200,])
        all_freqs <- rbind(all_freqs, table(true.tip.rates.extant))
    }
}
all_results$true_minus_est_eps <- all_results$eps.true-all_results$extinct.frac
all_results$true_minus_est_turnover <- all_results$turnover.true-all_results$turnover
all_results$misse <- "zero_k"
all_results_full <- rbind(all_results_full, all_results)


pairwise_summaries <- data.frame(matrix(nrow=0, ncol=8))
colnames(pairwise_summaries) <- c("file.index", "turnover.pairing", "median.turnover.actual", "fraction.turnover.difference.actual.positive", "number.comparisons", "median.turnover.difference", "ci.25.turnover.difference", "ci.75.turnover.difference")
for (index in sequence(max(all_results$file.index))) {
    print(paste(index, " of ", max(all_results$file.index)))
    f1 <- subset(all_results, file.index==index)
    distmat.turnover <- outer (f1$turnover, f1$turnover, `-`)
    distmat.turnover.compare_true <- outer (f1$turnover.true, f1$turnover.true, `paste`)
    delta.turnover.actual <- list()
    for (difference.index in sequence(length(compare_values.list))) {
        if(difference.index == 1){
            actual_turn <- c()
            for(compare.index in sequence(length(compare_values.list[[difference.index]]))){
                actual_turn <- c(actual_turn, as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])]))
            }
            compare_values_it <- compare_values.list[[difference.index]][1]
        }
        else{
            compare.index = 2
            compare_values_it <- compare_values.list[[difference.index]][1]
            actual_turn <- as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])])
        }
        delta.turnover.actual[[difference.index]] <- actual_turn
        pairwise_summaries[nrow(pairwise_summaries)+1,] <- c(index, compare_values.list[[difference.index]][2], median(f1$turnover), sum(actual_turn>0)/length(actual_turn), length(actual_turn), median(actual_turn), quantile(actual_turn, 0.25),quantile(actual_turn, 0.75))
    }
}

pairwise_summaries$misse <- "zero_k"
pairwise_summaries$true_difference_turnover <- LoopActualDifference(pairwise_summaries$turnover.pairing)
pairwise_summaries$median.turnover.difference <- as.numeric(pairwise_summaries$median.turnover.difference)
pairwise_summaries$fraction.turnover.difference.actual.positive <- as.numeric(pairwise_summaries$fraction.turnover.difference.actual.positive)

#ss <- boxplot(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, plot=FALSE)
#print(ss$stat)

zero_k_0.05 <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="GetSError"), c("Diff", "<50%"))
zero_k_0.05_med <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="median"), c("Diff", "Median"))

pairwise_summaries_full <- rbind(pairwise_summaries_full, pairwise_summaries)


######################################################################################################################################
######################################################################################################################################
### Strat Intervals -- ROW 3
######################################################################################################################################
######################################################################################################################################

all_results <- data.frame()
all_freqs <- data.frame()
for(i in seq_along(files)) {
    try(rm(misse.mk.strat.tip))
    try(rm(true.tip.rates))
    load(files[i])
    if(!is.null(misse.mk.strat.tip)) {
        true.tip.rates.extant <- true.tip.rates[names(true.tip.rates) %in% misse.mk.strat.tip$taxon]
        misse.mk.strat.tip$turnover.true <- turnover[1+unname(true.tip.rates[misse.mk.strat.tip$taxon])]
        misse.mk.strat.tip$eps.true <- eps[1+unname(true.tip.rates[misse.mk.strat.tip$taxon])]
        misse.mk.strat.tip$state.true <- true.tip.rates[misse.mk.strat.tip$taxon]
        misse.mk.strat.tip$file.index <- i
        misse.mk.strat.tip$file.name <- files[i]
        all_results <- rbind(all_results, misse.mk.strat.tip[1:200,])
        all_freqs <- rbind(all_freqs, table(true.tip.rates.extant))
    }
}
all_results$true_minus_est_eps <- all_results$eps.true-all_results$extinct.frac
all_results$true_minus_est_turnover <- all_results$turnover.true-all_results$turnover
all_results$misse <- "fstrat_int"
all_results_full <- rbind(all_results_full, all_results)


pairwise_summaries <- data.frame(matrix(nrow=0, ncol=8))
colnames(pairwise_summaries) <- c("file.index", "turnover.pairing", "median.turnover.actual", "fraction.turnover.difference.actual.positive", "number.comparisons", "median.turnover.difference", "ci.25.turnover.difference", "ci.75.turnover.difference")
for (index in sequence(max(all_results$file.index))) {
    print(paste(index, " of ", max(all_results$file.index)))
    f1 <- subset(all_results, file.index==index)
    distmat.turnover <- outer (f1$turnover, f1$turnover, `-`)
    distmat.turnover.compare_true <- outer (f1$turnover.true, f1$turnover.true, `paste`)
    delta.turnover.actual <- list()
    for (difference.index in sequence(length(compare_values.list))) {
        if(difference.index == 1){
            actual_turn <- c()
            for(compare.index in sequence(length(compare_values.list[[difference.index]]))){
                actual_turn <- c(actual_turn, as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])]))
            }
            compare_values_it <- compare_values.list[[difference.index]][1]
        }
        else{
            compare.index = 2
            compare_values_it <- compare_values.list[[difference.index]][1]
            actual_turn <- as.numeric(distmat.turnover[which(distmat.turnover.compare_true==compare_values.list[[difference.index]][compare.index])])
        }
        delta.turnover.actual[[difference.index]] <- actual_turn
        pairwise_summaries[nrow(pairwise_summaries)+1,] <- c(index, compare_values.list[[difference.index]][2], median(f1$turnover), sum(actual_turn>0)/length(actual_turn), length(actual_turn), median(actual_turn), quantile(actual_turn, 0.25),quantile(actual_turn, 0.75))
    }
}

pairwise_summaries$misse <- "fstrat_int"
pairwise_summaries$true_difference_turnover <- LoopActualDifference(pairwise_summaries$turnover.pairing)
pairwise_summaries$median.turnover.difference <- as.numeric(pairwise_summaries$median.turnover.difference)
pairwise_summaries$fraction.turnover.difference.actual.positive <- as.numeric(pairwise_summaries$fraction.turnover.difference.actual.positive)

#ss <- boxplot(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, plot=FALSE)
#print(ss$stat)

strat_int_0.05 <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="GetSError"), c("Diff", "<50%"))
strat_int_0.05_med <- setNames(aggregate(pairwise_summaries$fraction.turnover.difference.actual.positive~pairwise_summaries$true_difference_turnover, FUN="median"), c("Diff", "Median"))

pairwise_summaries_full <- rbind(pairwise_summaries_full, pairwise_summaries)




library(tidyverse)
library(Metrics)
library(MuMIn)
library(geiger)
library(phangorn)
library(TreeSim)
library(beepr)
library(ggplot2)
library("scales")
library(bde)
library(ggrepel)


######################################################################################################################################
######################################################################################################################################
### FUNCTIONS FOR SIMULATING TREES TO EXAMINE K VS M SAMPLING
######################################################################################################################################
######################################################################################################################################
GetSummary <- function(ntip, true.mu, true.lambda, true.turnover, true.eps, true.psi) {
    phy <- TreeSim::sim.bd.taxa(n = ntip, numbsim = 1, lambda = true.lambda, mu = true.mu)[[1]]
    fossils <- hisse:::GetFossils(phy, psi=true.psi)
    fossil.set <- hisse::ProcessSimSample(phy=phy, f=fossils)
    local_results <- data.frame(mu=true.mu, lambda=true.lambda, eps=true.eps, turnover=true.turnover, netdiv=true.lambda-true.mu, psi=true.psi, ntip=ntip, m.type=Ntip(fossil.set$phy)-ntip, k.type=dim(fossil.set$k.samples)[1], k.ratio=dim(fossil.set$k.samples)[1]/((Ntip(fossil.set$phy)-ntip)+dim(fossil.set$k.samples)[1]), stringsAsFactors = FALSE)
    return(local_results)
}


results <- data.frame()
true.turnover.vector <- seq(.1, 1.0, 0.1)
true.eps.vector <- seq(0.05, 0.95, 0.1)
true.psi.vector <- c(0.01, 0.05, 0.10, 0.20)
grid <- expand.grid(true.turnover.vector, true.eps.vector, true.psi.vector)
ntip.vector=c(200)
nrep <- 20
for(grid_index in 1:dim(grid)[1]){
    print(grid_index)
    true.turnover <- grid[grid_index,1]
    true.eps <- grid[grid_index,2]
    true.psi <- grid[grid_index,3]
    true.lambda <- true.turnover/ (1 + true.eps)
    true.mu <- (true.turnover * true.eps) / (1 + true.eps)
    ntip <- ntip.vector[1]
    for (rep_index in sequence(nrep)) {
        local_results <- data.frame()
        good <- FALSE
        while(good == FALSE){
            print(rep_index)
            try(local_results <- GetSummary(ntip=ntip, true.mu=true.mu, true.lambda=true.lambda, true.turnover=true.turnover, true.eps=true.eps, true.psi=true.psi))
            if(nrow(local_results)>0) {
                results <- plyr::rbind.fill(results, local_results)
                good <- TRUE
            }
        }
    }
}
write.csv(results, file="m_vs_kFINAL.csv")
beep(8)



######################################################################################################################################
######################################################################################################################################
### Generate Figure 3
######################################################################################################################################
######################################################################################################################################

userID  <-  Sys.info()['user']

#set dirs based on userID
switch(userID,
"jeremybeaulieu" = { # assume run from dir selon_work/
    contour.dir <- "/Users/jeremybeaulieu/Fossils_impact_BO/Single_Rate/Contour/Surfaces/";
    quick.sim.dir <- "/Users/jeremybeaulieu/Fossils_impact_BO/Single_Rate/Contour/QuickStratSim/";
    out.dir <- "/Users/jeremybeaulieu/Fossils_impact_BO/Tables_Figures/";},
)


data <- read.delim("m_vs_kFINAL.txt")
data$k_over_k_plus_m <- data$k.type/(data$m.type+data$k.type)
data$m_over_k <- data$m.type/data$k.type
data$k_over_m <- data$k.type/data$m.type
data$k_over_m[!is.finite(data$k_over_m)] <- 1000*max(data$k_over_m[is.finite(data$k_over_m)], na.rm=TRUE)
data$k_over_m_breaks <- cut(data$k_over_m,breaks = c(-Inf,1/10,  0.8, 1.2,  10, Inf),right = FALSE)
eps_values <- sort(unique(data$eps))
turnover_values <- sort(unique(data$turnover))
psi_values <- sort(unique(data$psi))

data_hist <- hist(data$k_over_k_plus_m, freq=TRUE, breaks=99)
data_hist$percentage <- 100*data_hist$counts/sum(data_hist$counts)
data_hist_df <- data.frame(mids=data_hist$mids, percentage=data_hist$percentage)

empirical_km <- rbind(
data.frame(paper="  Empirical: Zhang et al. 2016 Median", k=0, m=45),
data.frame(paper="  Empirical: Gavryushkina et al. 2017", k=18, m=17),
data.frame(paper="  Empirical: Slater et al. 2017", k=1, m=70),
data.frame(paper="  Empirical: Pyron 2017", k=0, m=62)
)
empirical_km$k_over_k_plus_m <- empirical_km$k/(empirical_km$k+empirical_km$m)
empirical_km$k_over_k_plus_m[1] <- 0.01 #since the median was slightly above 0
empirical_km$height <- NA
for (i in sequence(nrow(empirical_km))) {
    empirical_km$height[i] <- data_hist$percentage[which.min(abs(data_hist$mid - empirical_km$k_over_k_plus_m[i]))]
}


ggplot(data, aes(x=k_over_k_plus_m)) + geom_histogram(bins=100) + scale_y_continuous(labels = scales::percent)

g <- ggplot(data_hist_df, aes(x=mids, y=percentage)) + geom_bar(stat="identity") + xlab("Fraction of fossils that are k-type") + ylab("Percentage of simulated datasets")
g <- g + geom_point(data=empirical_km, aes(col="red", x=k_over_k_plus_m, y=height), show.legend = FALSE) +  geom_text_repel(data=empirical_km, aes(label=paper, x=k_over_k_plus_m, y=height, angle=90), min.segment.length = 0, vjust=1)
g <- g + geom_segment(aes(x=0.5, xend=0.75, y=0, yend=0, colour='darkblue'),show.legend = FALSE) + theme_light()
print(g)
ggsave(paste0(out.dir, "Figure3.pdf"), width=6, height=6)


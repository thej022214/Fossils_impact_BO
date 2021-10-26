
library(hisse)
library(paleotree)
library(interp)
library(lhs)
library(TreeSim)
library(data.table)
######################################################################################################################################
######################################################################################################################################
### FUNCTIONS FOR GETTING CONTOUR COMPONENTS
######################################################################################################################################
######################################################################################################################################

contourSearchPoints <- function(variables, lower, upper, nreps){
    #Creates a latin square design of the parameter space:
    X <- randomLHS(nreps, variables)
    param.points <- matrix(0, nrow=nreps, ncol=variables)
    param.points[,1] <- qunif(X[,1], lower[1], upper[1])
    param.points[,2] <- qunif(X[,2], lower[2], upper[2])
    return(param.points)
}

#c1 and c2 are taken from Stadler 2010, pg. 398
c1 <- function(lambda, mu, psi){
    return(abs(sqrt(((lambda - mu - psi)^2)+(4*lambda*psi))))
}

c2 <- function(lambda, mu, psi, rho){
    return(-(lambda-mu-(2*lambda*rho)-psi)/c1(lambda, mu, psi))
}

#This is Master equation from Stalder 2010, pg. 398.
p_0 <- function(x, lambda, mu, psi, rho){
    c1_val <- c1(lambda=lambda, mu=mu, psi=psi)
    c2_val <- c2(lambda=lambda, mu=mu, psi=psi, rho=rho)
    return(log((lambda + mu + psi + c1_val * ((exp(-c1_val*x) * (1-c2_val)-(1+c2_val))/(exp(-c1_val*x) * (1-c2_val)+(1+c2_val))))/(2*lambda)))
}


#This is Master equation from Stadler 2010, pg. 398.
p_one <- function(x,lambda, mu, psi, rho){
    c1_val <- c1(lambda=lambda, mu=mu, psi=psi)
    c2_val <- c2(lambda=lambda, mu=mu, psi=psi, rho=rho)
    return(log((4*rho)/(2*(1-c2_val^2) + exp(-c1_val*x)*(1-c2_val)^2 + exp(c1_val*x) * (1+c2_val)^2)))
}


#Test 1 -- Constant case, no phi
GetLambdaMuNoPsi <- function(x, phy, convert=TRUE){
    #x <- exp(x)
    if(convert == TRUE){
        lambda <- x[1] / (1 + x[2])
        mu <- (x[2] * x[1]) / (1 + x[2])
    }else{
        lambda <- x[1]
        mu <- x[2]
    }
    gen <- hisse:::FindGenerations(phy)
    dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy, f=1, hidden.states=3)
    nb.tip <- Ntip(phy)
    nb.node <- phy$Nnode
    model.vec <- c(lambda+mu, mu/lambda, rep(0,51))
    cache = hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=1, fixed.eps=NULL, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-500), ode.eps=0)#
    cache$psi = 0
    loglik <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, fossil.taxa=NULL)
    
    return(-loglik)
}


#Test 2 -- Constant case, estimate phi
GetPsiRestFix <- function(x, turn, ef, rho, n, m, k, x_times, y_times, type="stad", sim.obj=NULL, k.samples, convert=TRUE){
    psi <- x
    if(convert == TRUE){
        lambda <- turn / (1 + ef)
        mu <- (ef * turn) / (1 + ef)
    }else{
        lambda <- turn
        mu <- ef
    }
    if(type == "stad"){
        loglik <- NA
        try(loglik <- (((n+m-2) * log(lambda)) + ((k+m) * log(psi))) - log(1-exp(p_0(max(x_times),lambda,mu,psi=0,rho)))*2 + p_one(max(x_times), lambda, mu, psi, rho) + sum(p_one(x_times, lambda,mu,psi,rho)) + (sum(p_0(y_times,lambda,mu,psi,rho)) - sum(p_one(y_times,lambda,mu,psi,rho))), silent=TRUE)
        if(is.na(loglik[1])){
            return(10000000)
        }else{
            return(-loglik)
        }
    }
    if(type == "hisse"){
        if(!is.null(k.samples)){
            phy <- hisse:::AddKNodes(sim.obj$phy, sim.obj$k.samples)
            nb.tip <- Ntip(phy)
            nb.node <- phy$Nnode
            dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy, f=1, hidden.states=1)
            model.vec <- c(lambda+mu, mu/lambda, rep(0,51))
            cache = hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=1, fixed.eps=NULL, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-500), ode.eps=0)#
            cache$psi = psi
            fossil.taxa <- which(dat.tab$branch.type == 1)
            gen <- hisse:::FindGenerations(phy)
            k.samples <- hisse:::GetKSampleMRCA(phy, k.samples)
            
            loglik <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, fossil.taxa=fossil.taxa, node=k.samples$node, fix.type=k.samples$type)
            return(-loglik)
        }else{
            phy <- sim.obj$phy
            nb.tip <- Ntip(phy)
            nb.node <- phy$Nnode
            dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy, f=1, hidden.states=1)
            model.vec <- c(lambda+mu, mu/lambda, rep(0,51))
            cache <- hisse:::ParametersToPassMiSSE(model.vec=model.vec, hidden.states=1, fixed.eps=NULL, nb.tip=nb.tip, nb.node=nb.node, bad.likelihood=exp(-500), ode.eps=0)#
            cache$psi <- psi
            fossil.taxa <- which(dat.tab$branch.type == 1)
            gen <- hisse:::FindGenerations(phy)

            loglik <- hisse:::DownPassMisse(dat.tab=dat.tab, cache=cache, gen=gen, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, fossil.taxa=fossil.taxa, node=NULL, fix.type=NULL)
            return(-loglik)
        }
    }
}

#Function that calls the test1 function, GetLambdaMuNoPsi()
contourSearchMiSSEfos1 <- function(sim.obj, param.points, psi, rho, nreps, type, convert) {
    
    phy <- sim.obj$phy
    dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy, f=rho, hidden.states=1)
    gen <- hisse:::FindGenerations(phy)
    
    res <- matrix(,nreps,3)
    for(nrep.index in 1:nreps){
        print(nrep.index)
        fixed.pars <- c(param.points[nrep.index,1], param.points[nrep.index,2])
        out <- GetLambdaMuNoPsi(x=fixed.pars, phy=phy, convert=convert)
        res[nrep.index,] <- c(-out, fixed.pars[1], fixed.pars[2])
    }
    return(res)
}

#Function that calls the test2 function, GetPsiRestFix()
contourSearchMiSSEfos2 <- function(sim.obj, param.points, rho, nreps, type, drop.k, convert) {
    
    phy <- sim.obj$phy
    sample.points <- sim.obj$k.samples
    
    dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy, f=rho, hidden.states=1)
    gen <- hisse:::FindGenerations(phy)
    fossil.taxa <- which(dat.tab$branch.type == 1)
    #This is for starting values:
    fossil.ages <- dat.tab$TipwardAge[which(dat.tab$branch.type == 1)]
    split.times <- dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))[-c(1:Ntip(phy))]
    m <- length(fossil.taxa)
    if(drop.k == TRUE){
        k <- 0
        k.samples <- NULL
    }else{
        k <- length(sample.points$timefrompresent)
        k.samples <- sample.points
    }
    res <- matrix(,nreps,3)
    for(nrep.index in 1:nreps){
        print(nrep.index)
        fixed.pars <- c(param.points[nrep.index,1], param.points[nrep.index,2])
        out <- optimize(GetPsiRestFix, interval=c(0,0.2), turn=fixed.pars[1], ef=fixed.pars[2], rho=rho, n=Ntip(phy)-m, m=m, k=k, x_times=split.times, y_times=fossil.ages, type=type, sim.obj=sim.obj, k.samples=k.samples, convert=convert)
        res[nrep.index,] <- c(-out$objective, fixed.pars[1], fixed.pars[2])
    }
    return(res)
}

#Function getting contour data
GetContour <- function(sim.obj, focal.params=c("turnover", "extinct. frac."), psi=NULL, focal.params.lower=c(0,0), focal.params.upper=c(5,5), nreps=1000, rho=1, type=c("stad", "hisse"), test.no=c("test1", "test2"), drop.k=FALSE, convert=TRUE){
    
    param.points <- contourSearchPoints(variables=2, lower=focal.params.lower, upper=focal.params.upper, nreps=nreps)
    if(test.no == "test1"){
        surface.data <- contourSearchMiSSEfos1(sim.obj=sim.obj, param.points=param.points, rho=rho, nreps=nreps, type=type, convert=convert)
    }
    if(test.no == "test2"){
        surface.data <- contourSearchMiSSEfos2(sim.obj=sim.obj, param.points=param.points, rho=rho, nreps=nreps, type=type, drop.k=drop.k, convert=convert)
    }

    obj <- list(surface.data=surface.data, focal.params=focal.params, focal.params.lower=focal.params.lower, focal.params.upper=focal.params.upper)
    
    return(obj)
}

#Function that plots contour surface
plot.contour <- function(x, mle.point=NULL, levels=c(0:20*0.1), xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL, col=grey.colors(21, start=0, end=1), ...){
    
    if(is.null(xlab)){
        xlab = x$focal.params[1]
    }
    if(is.null(ylab)){
        ylab = x$focal.params[2]
    }
    if(is.null(xlim)){
        xlim = c(x$focal.params.lower[1], x$focal.params.upper[1], 1)
    }
    if(is.null(ylim)){
        ylim = c(x$focal.params.lower[2], x$focal.params.upper[2], 1)
    }
    
    mydata <- data.frame(x=matrix(x$surface.data[,2],ncol=1),y=matrix(x$surface.data[,3],ncol=1),z=matrix(x$surface.data[,1],ncol=1))
    mydata$z <- (-1)*mydata$z
    mydata$z <- mydata$z-min(mydata$z)
    
    interp.res <- interp(x=mydata$x, y=mydata$y, z=mydata$z, xo=seq(min(mydata$x), max(mydata$x),length = 400), yo=seq(min(mydata$y), max(mydata$y),length = 400), duplicate=FALSE)
    plot(NA,xlab="", ylab="", frame=FALSE, axes=FALSE, xaxs="i", yaxs="i", ylim=ylim[1:2], xlim=xlim[1:2], ...)
    .filled.contour(interp.res$x, interp.res$y, interp.res$z, levels=levels, col=col)
    par(tck=.01)
    axis(2, at = seq(ylim[1], ylim[2], by = ylim[3]), las=1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
    axis(1, at = seq(xlim[1], xlim[2], by = xlim[3]), las=1, lwd=1, labels=TRUE, mgp=c(.75,.5,0))
    if(!is.null(mle.point)){
        points(x=mle.point[1], y=mle.point[2], pch=19, col="blue")
    }
    title(xlab=xlab, line=2.5)
    title(ylab=ylab, line=2)
}


######################################################################################################################################
######################################################################################################################################
### Generate Figure 1 Data
######################################################################################################################################
######################################################################################################################################

### Getting the tree for constant BD case###
#set.seed(42)
#done=FALSE
#while(done == FALSE){
#    try(sim.tab <- hisse::SimulateHisse(turnover=rep(.4+.3,2), eps=rep(.3/.4,2), max.taxa=200, x0=0, transition.rates=matrix(c(NA, 0.005, 0.005, NA), nrow=2), nstart=2))
#    if(sim.tab$n.surviving > 0){
#        try(phy <- SimToPhylo(sim.tab, include.extinct=TRUE))
#        done = TRUE
#    }
#}
#f <- hisse:::GetFossils(phy, psi=0.005)
#pp <- hisse:::ProcessSimSample(phy, f)
#MiSSE(pp$phy, turnover=1, eps=1, includes.fossils=TRUE, k.sample=pp$k.samples, sann=FALSE)

#phy.extant <- geiger::drop.extinct(phy, tol=.Machine$double.eps^.50)
#MiSSE(phy.extant, turnover=1, eps=1, sann=FALSE)

#save(phy, f, pp, file="simTreeContour005.Rsave")
########################

setDTthreads(threads=1)

load("simTreeContour005.Rsave")

pp.red <- pp
pp.red$phy <- geiger::drop.extinct(pp$phy, tol=.Machine$double.eps^.50)
cont.info1 <- GetContour(pp.red, focal.params.lower=c(0,0), focal.params.upper=c(2, 1), nreps=5000, type="hisse", test.no="test1")
cont.info1$surface.data <- cont.info1$surface.data[is.finite(cont.info1$surface.data[,1]),]
cont.info1$focal.params <- c(expression(tau), expression(epsilon))
save(cont.info1, file="test1_5000points005.Rsave")

cont.info2 <- GetContour(pp, focal.params.lower=c(0,0), focal.params.upper=c(2, 1), nreps=5000, type="hisse", test.no="test2")
cont.info2$surface.data <- cont.info2$surface.data[is.finite(cont.info2$surface.data[,1]),]
cont.info2$focal.params <- c(expression(tau), expression(epsilon))
save(cont.info2, file="test2_5000points005.Rsave")

cont.info3 <- GetContour(pp, focal.params.lower=c(0,0), focal.params.upper=c(2, 1), nreps=5000, type="hisse", test.no="test2", drop.k=TRUE)
cont.info3$surface.data <- cont.info3$surface.data[is.finite(cont.info3$surface.data[,1]),]
cont.info3$focal.params <- c(expression(tau), expression(epsilon))
save(cont.info3, file="test3_5000points005.Rsave")


#f <- hisse:::GetFossils(phy, psi=0.05)
#pp <- hisse:::ProcessSimSample(phy, f)
#save(phy, f, pp, file="simTreeContour050.Rsave")
#load("simTreeContour050.Rsave")

#cont.info4 <- GetContour(pp, focal.params.lower=c(0,0), focal.params.upper=c(2, 1), nreps=5000, type="hisse", test.no="test2")
#cont.info4$surface.data <- cont.info4$surface.data[is.finite(cont.info4$surface.data[,1]),]
#cont.info4$focal.params <- c(expression(tau), expression(epsilon))
#save(cont.info4, file="test2_5000points050.Rsave")

#cont.info5 <- GetContour(pp, focal.params.lower=c(0,0), focal.params.upper=c(2, 1), nreps=5000, type="hisse", test.no="test2", drop.k=TRUE)
#cont.info5$surface.data <- cont.info5$surface.data[is.finite(cont.info5$surface.data[,1]),]
#cont.info5$focal.params <- c(expression(tau), expression(epsilon))
#save(cont.info5, file="test3_5000points050.Rsave")


########################

# Rate shift -- turnover in regime A and B.
#               net div  in regime A and B.

########################


######################################################################################################################################
######################################################################################################################################
### Generate Figure 1 Plot
######################################################################################################################################
######################################################################################################################################

#pdf("Figure1.pdf", width=12, height=8)
#par(mfcol=c(2,3),mar=c(4,4.5,0.5,0.5), oma=c(1.5,2,1,1))

#load("test1_5000points005.Rsave")
#plot.contour(cont.info1, mle.point=c(0.4+0.3, 0.3/0.4), levels=c(0:200*0.01), ylim=c(0,1,.1), xlim=c(0,1.6,.2), col=grey.colors(21, start=0, end=1))
#mtext("A",side=3, line=0, adj=0, cex=1)
#text(.1, .9, expression(extant~only), pos=4)

#plot(0,type='n',axes=FALSE,ann=FALSE)

#load("test2_5000points005.Rsave")
#plot.contour(cont.info2, mle.point=c(0.4+0.3, 0.3/0.4), levels=c(0:200*0.01), ylim=c(0,1,.1), xlim=c(0,1.6,.2), col=grey.colors(21, start=0, end=1))
#mtext("B",side=3, line=0, adj=0, cex=1)
#text(.1, .9, expression(psi~estimated~m+k), pos=4)

#load("test2_5000points040.Rsave")
#plot.contour(cont.info4, mle.point=c(0.4+0.3, 0.3/0.4), levels=c(0:200*0.01), ylim=c(0,1,.1), xlim=c(0,1.6,.2), col=grey.colors(21, start=0, end=1))
#mtext("D",side=3, line=0, adj=0, cex=1)
#text(.1, .9, expression(psi~estimated~m+k), pos=4)

#load("test3_5000points005.Rsave")
#plot.contour(cont.info3, mle.point=c(0.4+0.3, 0.3/0.4), levels=c(0:200*0.01), ylim=c(0,1,.1), xlim=c(0,1.6,.2), col=grey.colors(21, start=0, end=1))
#mtext("C",side=3, line=0, adj=0, cex=1)
#text(.1, .9, expression(psi~estimated~m~only), pos=4)

#load("test3_5000points040.Rsave")
#plot.contour(cont.info5, mle.point=c(0.4+0.3, 0.3/0.4), levels=c(0:200*0.01), ylim=c(0,1,.1), xlim=c(0,1.6,.2), col=grey.colors(21, start=0, end=1))
#mtext("E",side=3, line=0, adj=0, cex=1)
#text(.1, .9, expression(psi~estimated~m~only), pos=4)

#dev.off()




#pdf("Figure2.pdf",  height=6, width = 7)
#par(mfcol=c(2,2),mar=c(4,4,0.5,0.5), oma=c(1.5,2,1,1))
#plot(results$m, results$k, ylim=c(0,200), ylab="", xlab="", pch=19, xlim=c(0,200), axes=FALSE, cex=0.25)
#title(ylab = expression(italic("k")~fossils), line=2.5)
#title(xlab = expression(italic("m")~fossils), line=2)
#par(tck=.01)
#axis(2, at = seq(0, 200, by = 50), las =1, lwd=1, labels=TRUE, mgp=c(1,.5,0))
#axis(1, at = seq(0, 200, by = 50), las =1, lwd=1, labels=TRUE, mgp=c(1,.5,0))
#abline(0, 1, lty=2, col="blue")
#dev.off()






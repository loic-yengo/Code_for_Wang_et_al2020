library(ggplot2)
library(reshape)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")

load("Figure_3.Rdata")

d$ExpLoa    <- 100* d$pred_loa_II/d$loa
d$ExpLoa_se  <- 1.96*100*sqrt( (d$RA_se^2) * (d$pred_loa_II^2)/(d$loa^4) )
d1 <- d[d$Sig_pval == 1,]

nt   <- 8
Cols <- c("coral1","dodgerblue","maroon","khaki3","violet","goldenrod","seagreen","lightblue")
f <- function(pop,legend=FALSE){
  u  <- d[which(d$pop==pop),]
  l  <- which(u[,"sigLoss"]==1) 
  ym <- 100#max(100,max(u[,"RA"] + u[,"RA_se"]))
  if(legend){
    #ym <- 150
  }
  Ax <- barplot(c(t(cbind(0,u[l,c("pred_RA_II","RA")]))),col=rep(Cols[l],each=3),axes=FALSE,
                ylab="Relative accuracy (%)",
                density = rep(c(0,100,40),nt),border=0,main=pop,ylim=c(0,ym))
  kx <- seq(3,length(Ax),by=3)
  ax <- sapply(kx,function(k) mean(Ax[(k-1):k]))
  axis(1,at=ax,labels=u[l,"trait"],las=2)
  axis(2,at=seq(0,100,by=20))
  for(k in 1:length(Ax)){
    arrows(Ax[kx[k]],u[l,"RA"][k]-u[l,"RA_se"][k],
           Ax[kx[k]],u[l,"RA"][k]+u[l,"RA_se"][k],
           col=Cols[k],code=3,len=0.05,angle=90)
  }
  abline(h=c(50,100),col="grey",lty=2)
  if(legend){
    legend(0,90,legend=c("Predicted (LD+MAF)","Observed"),horiz = TRUE,fill="black",density = c(100,40),
           box.lty=0,x.intersp = 0.5)
  }
}

fx <- function(pop){
  u    <- d[which(d$pop==pop),]
  l    <- which(u[,"sigLoss"]==1)  
  ymax <- +160 #max(u[,"ExpLoa"] + u[,"ExpLoa_se"],na.rm=T)
  ymin <- -10  #min(u[,"ExpLoa"] - u[,"ExpLoa_se"],na.rm=T)
  par(mar=c(5,6,3,3))
  Ax <- barplot(u[l,"ExpLoa"],col=Cols[l],axes=FALSE,
                ylab="Proportion of LOA explained\nby LD and MAF (%)",
                border=0,main="",ylim=c(ymin,ymax))
  axis(1,at=Ax,labels=u[l,"trait"],las=2)
  axis(2,at=seq(0,100,by=50))
  for(k in 1:length(Ax)){
    arrows(Ax[k],u[l,"ExpLoa"][k]-u[l,"ExpLoa_se"][k],
           Ax[k],u[l,"ExpLoa"][k]+u[l,"ExpLoa_se"][k],
           col=Cols[l][k],code=3,len=0.05,angle=90)
  }
  abline(h=c(50,100),col="grey",lty=2)
}


png("Figure_3.png",width=3000,height=2000,res=300)
op <- par(mfrow=c(2,3))
f("SAS");mtext("A",at=0,side=3,line=+1,font=2,cex=1.25)
f("EAS");mtext("B",at=0,side=3,line=+1,font=2,cex=1.25)
f("AFR",legend=TRUE);mtext("C",at=0,side=3,line=+1,font=2,cex=1.25)
fx("SAS");mtext("D",at=0,side=3,line=-3,font=2,cex=1.25)
fx("EAS");mtext("E",at=0,side=3,line=-3,font=2,cex=1.25)
fx("AFR");mtext("F",at=0,side=3,line=-3,font=2,cex=1.25)
par(op)
dev.off()
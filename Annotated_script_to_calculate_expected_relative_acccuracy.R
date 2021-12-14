## On my computer
setwd("~/Desktop/Papers/EA42021/Transferrability/")
e    <- read.table("UKB_GWAS/EUR.ldcor.pairs.gz",h=T,stringsAsFactors = F) # output from ldcorpair with an EUR sample
a    <- read.table("UKB_GWAS/AFR.ldcor.pairs.gz",h=T,stringsAsFactors = F) # output from ldcorpair with an AFR sample

e$ID <- paste0(e$SNP1,":",e$SNP2) # create an single ID for pairs of SNPs
a$ID <- paste0(a$SNP1,":",a$SNP2) # create an single ID for pairs of SNPs

e    <- e[which(e$R_LD_SQ>0.45),] # filter pairs with an LD squared correlation > 0.45
l    <- intersect(e$ID,a$ID)
a    <- a[which(a$ID%in%l),]
e    <- e[which(e$ID%in%l),]
l    <- which(a$A1_1==e$A1_1 & a$A2_1==e$A2_1 & a$A1_2==e$A1_2 & a$A2_2==e$A2_2 & a$SNP1!=a$SNP2 & e$SNP1!=e$SNP2) # QC pairs on consistency of alleles
a    <- a[l,]
e    <- e[l,]

# reading sumstatistics
# format: SNP | ALLELE | BETA
gwas <- read.table("UKB_GWAS/507_GWS.rsid.sumstat",h=T,stringsAsFactors = F,sep="\t") 
rownames(gwas) <- gwas$SNP
bg             <- gwas$BETA
names(bg)      <- gwas$SNP


# R1 statistic (contribution of LD differences) is a ratio = num / den
# First calculate all elements required to calculate num and den
cxy   <- aggregate(a$R_LD*e$R_LD,by=list(a$SNP1),FUN=mean)
vx    <- aggregate(e$R_LD_SQ_CORR,by=list(e$SNP1),FUN=mean)
vy    <- aggregate(a$R_LD_SQ_CORR,by=list(a$SNP1),FUN=mean)
hx    <- aggregate(e[,"FREQ_A1_1"]*(1-e[,"FREQ_A1_1"]),by=list(e$SNP1),FUN=mean)
hy    <- aggregate(a[,"FREQ_A1_1"]*(1-a[,"FREQ_A1_1"]),by=list(a$SNP1),FUN=mean)

## Calculate num and den
num   <- sum(cxy[,2] * sqrt(hy[,2] / hx[,2]))
den   <- sum(vx[,2])

R1    <- (num/den)^2
R2    <- mean(bg[hx[,1]]*bg[hx[,1]]*hx[,2],na.rm = TRUE) / mean(bg[hy[,1]]*bg[hy[,1]]*hy[,2],na.rm=TRUE)

RAafr <- R1*R2
RAafr # Expected relative accuracy due to MAF and LD differences between ancestries

## This function "f" is to defined to calculate leave-one-chromosome-out jackkniffe standard error
f <- function(k){
  A     <- a[-which(a$CHR==k),]
  E     <- e[-which(e$CHR==k),]
  cxy   <- aggregate(A$R_LD*E$R_LD,by=list(A$SNP1),FUN=mean)
  vx    <- aggregate(E$R_LD_SQ_CORR,by=list(E$SNP1),FUN=mean)
  vy    <- aggregate(A$R_LD_SQ_CORR,by=list(A$SNP1),FUN=mean)
  hx    <- aggregate(E[,"FREQ_A1_1"]*(1-E[,"FREQ_A1_1"]),by=list(E$SNP1),FUN=mean)
  hy    <- aggregate(A[,"FREQ_A1_1"]*(1-A[,"FREQ_A1_1"]),by=list(A$SNP1),FUN=mean)
  num   <- sum(cxy[,2] * sqrt(hy[,2] / hx[,2]))
  den   <- sum(vx[,2])
  R1    <- (num/den)^2
  R2    <- sum(bg[hx[,1]]*bg[hx[,1]]*hx[,2]) / sum(bg[hy[,1]]*bg[hy[,1]]*hy[,2])
  RAafr <- R1*R2
  RAafr
}

RAset <- sapply(1:22,f) # calculating the RAafr after excluding each of the 22 autosomal chromosomes
SE    <- sqrt( sum( (1-1/22) *( (RAset-RAafr)^2 ) ) ) # Jackkniffe standard error
CI    <- RAafr + 1.96*SE*c(-1,1) # Approximate 95% confidence interval


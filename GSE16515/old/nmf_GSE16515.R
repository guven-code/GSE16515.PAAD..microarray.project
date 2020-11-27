library(NMF)
library(tictoc)

# DEdata=read.csv("GSE16515_DE.csv")
# 
# sampleNames<-colnames(DEdata[,-1]) #this gives just numbers of the expressed genes
# sampleNames[1:3] #checking patients 1 to 3
# 
# probeIDs = DEdata[,1] # okay finally this gives the geneNames
# head(probeIDs)
# 
# #BRCA=BRCA[,-1] #this drops the first column; geneNames
# #newBRCA=BRCA[1:100,1:100]
# 
# colnames(DEdata)=probeIDs
# 
# DEdata.mat<-data.matrix(DEdata,rownames.force = NA)
# 
# DEdata.mat=DEdata.mat[,-1]
# rownames(DEdata.mat)=probeIDs


#saveRDS(DEdata.mat,"DEdataNMF.rds")
DEdataNMF<-readRDS("DEdataNMF.rds")


ptm <- proc.time()
rankSurvey.DEdata=nmfEstimateRank(DEdataNMF,r=10:50,nrun=30,.options='t',seed=12345)
saveRDS(rankSurvey.DEdata,"rankSurvey10_50nrun30.Rds")
proc.time() - ptm

rankSurvey<-readRDS("rankSurvey2_45nrun30.Rds")
randomRankSurvey<-readRDS("randomRankSurvey2_45nrun30.Rds")

# It requested for one r=3 for instance 8.5 min to run 
# all the factorizations sequentially, 4.8 min when using 
# multi-core computation alone (see the hardware specification above), 
# and 2.5 min when distributed over 4 quadri-core nodes on a HPC cluster . 
# Besides the memory used by the R session itself, a single NMF 
# run of Brunet's algorithm (with r = 3) required on average 25Mb.
# The panCancer data set and the fitted rank-3 factorization used 850Kb and 450Kb 
# respectively. 

#1.4 GHz Intel Core i5
#4 GB 1600 MHz DDR3


L=length(rankSurvey$measures$rank)

x=rankSurvey.DEdata$measures$rank
y=rankSurvey.DEdata$measures$rss

z=rankSurvey.DEdata$measures$silhouette.basis
d=rankSurvey.DEdata$measures$dispersion

pdf("nmf.plots/2plots.pdf")
plot(x,d,type="l",col="red")
par(new=TRUE)
plot(x,z,type="l",col="blue")
dev.off()
#plot(rankSurvey,randomRankSurvey)
pdf("nmf.plots/4plots.pdf")
par(mfrow=c(2,2)) 
plot(rankSurvey,"cophenetic")
plot(rankSurvey,"silhouette")
plot(rankSurvey,"rss")
dev.off()

plot(x,z,type="l")#,ylim=c(3,10))
lo <- loess(z~x) #Local Polynomial Regression Fitting
xl <- seq(min(x),max(x), (max(x) - min(x))/1000)
out = predict(lo,xl) #linear regression 
lines(xl, out, col='red', lwd=2)

infl <- c(FALSE, diff(diff(out)>0)!=0)

xl[infl ]

round(xl[infl ])

points(xl[infl ], out[infl ], col="blue")






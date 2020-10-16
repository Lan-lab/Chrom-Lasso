library(MASS)
x<-as.matrix(read.table("randomSamples.combined"))
m<-p.adjust(x[,4], method ="fdr", n = length(x[,4]))
n<-p.adjust(x[,4], method ="BY", n = length(x[,4]))
q<-cbind(x,m,n)
write.matrix(q,file="randomSamples.combined.fdr",sep="\t")
y<-x[x[,3]>=0,]
dim(y)
s<-p.adjust(y[,4], method ="fdr", n = length(y[,4]))
t<-p.adjust(y[,4], method ="BY", n = length(y[,4]))
o<-cbind(y,s,t)
write.matrix(o,file="randomSamples.combined.posFdr",sep="\t")

#> x<-as.matrix(read.table("randomSamples_chr22"))
#> hist(x[,4])
#> z<-p.adjust(x[,4], method ="fdr", n = length(x[,4]))
#> hist(z)
#> length(z[z<0.01])
#[1] 28926
#> length(z[z<0.05])
#[1] 43908
#> length(x[,4])
#[1] 180737
#> savehistory("fdrFromRandomSamples.r")
#> z<-p.adjust(x[,4], method ="BY", n = length(x[,4]))
#> hist(z)
#> length(z[z<0.05])
#[1] 24434
#> length(z[z<0.01])
#[1] 19862
#

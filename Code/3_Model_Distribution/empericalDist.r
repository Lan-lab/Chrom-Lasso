library(MASS)
pdf("Mouse.empericalDist.bin100.pdf")

s<-read.table("empericalDist.bin100")
allFrag<-s[1:10000,2]
x<-read.table("bin100.withNonloop")
#a<-log(x[1:1000,2])

b<-log(x[1:1000,1]*100+50)
print(sum(as.numeric(x[,2])))
c<-x[1:1000,2]
excludeFRpair<-x[1:10000,2]
FRpair<-(allFrag-excludeFRpair)/allFrag
plot(FRpair)
print(mean(FRpair[15:10000]))
c<-c/3*4
#c<-s[1:1000,2]
#d<-log(x[1:1000,1]*100+50)

y<-read.table("csDistanceDistr")
meanNum<-mean(y[,2])
c=c/y[,2]

glmfit<-glm(c~b, family="poisson")
print (glmfit)
print(coef(glmfit)[1])
print(coef(glmfit)[2])

mu<- (coef(glmfit)[1])
beta<- (coef(glmfit)[2])

#lmfit<-lm(a~b)
#print (lmfit)
#mu<- (coef(lmfit)[1])
#beta<- (coef(lmfit)[2])
head(c)

plot(b,log(c))
lines(c(0,-mu/beta),c(mu,0))
plot(b,c,ylim=c(0,10))
plot(b,log(c),xlim=c(4,12))
#lines(smooth.spline(b[11:1000],log(c[11:1000])))
#smspl<-(smooth.spline(b[11:1000],log(c[11:1000])))
#print(smspl$fit)
#summary(smspl)
lmobj<-lm(log(c) ~ poly(b, 7, raw=TRUE))
#lines(seq(4,12,0.01),predict.lm(lmobj,seq(4,12,0.01)))
lines(b, fitted(lmobj), col='red', type='b')
coefNames<-matrix(0.0, length(coef(lmobj)),1)
print(coefNames)
coefNames[1,1]<-0
for(i in 2:length(coef(lmobj)))
{
    coefNames[i,1]<-substr(names(coef(lmobj))[i],nchar(names(coef(lmobj))[i]),nchar(names(coef(lmobj))[i]))
}
print(coefNames)
coefMatrix<-cbind(coefNames,coef(lmobj))
print(coefMatrix)
write.matrix(coefMatrix, file="PolyCoef", sep=" ")
print(summary(lmobj))
print(lmobj)
#co<-coef(lmobj)
#co[1]+co[2]*7+co[3]*7^2+co[4]*7^3+co[5]*7^4+co[6]*7^5+co[7]*7^6+co[8]*7^7 
dev.off()

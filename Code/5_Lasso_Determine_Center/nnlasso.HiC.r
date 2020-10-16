library(lars)
library(glmnet)
library(glmpath)
library(covTest)
library(nnlasso)

args <- commandArgs(trailingOnly = TRUE)
regionFreq<-readLines(args[1])
distMat<-readLines(args[2])
totalRegions_freq<-length(regionFreq)/2
totalLines_dist<-length(distMat)
j=1 #distMatCounter

#pdf("test.pdf")
#for(i in 1:100)
for(i in 1:totalRegions_freq)
{
    regionHeader<-as.numeric((strsplit(regionFreq[2*i-1],"\t"))[[1]])
    freq<-as.numeric((strsplit(regionFreq[2*i],"\t"))[[1]])
    print(regionHeader[2])
    freqSize<-length(freq)

    distMatHeader<-as.numeric((strsplit(distMat[j],"\t"))[[1]])
    j=j+1

    dist<-matrix(0.0,0,freqSize)
    if(distMatHeader[3]==regionHeader[2])
    {
	dataFlag<-TRUE
	while(dataFlag)
	{
	    distData<-as.numeric((strsplit(distMat[j],"\t"))[[1]])
	    if(length(distData)==3 || j>totalLines_dist)
	    {
		dataFlag<-FALSE
	    }
	    else
	    {
		dist<-rbind(dist,distData)
		#print(dim(dist))
		j=j+1
	    }
	}
    }
    else
    {
	print("ERROR! Region number not equal to test number.")
    }
    dist<-t(dist)
    #print(dim(dist))
    #print(length(freq))
    #print(freq)
    #print(dist)
    nnlassoobj<-nnlasso(dist, freq, family="normal",intercept=FALSE)
    #print(coef(nnlassoobj))
    #cvobj<-cv.nnlasso(dist,freq,family="normal", k=5,nlambda=50,tau=1,plot=FALSE, errorbars=FALSE)
    cvobj<-tryCatch({
	cvobj<-cv.nnlasso(dist,freq,family="normal", k=5,nlambda=50,tau=1,plot=FALSE, errorbars=FALSE)
	optCoef<-predict(nnlassoobj,mode="lambda",at=cvobj$lambda)
	cat(optCoef[2:length(optCoef)])
	cat("\n")
    }, error=function(err){
	cvobj<-cv.nnlasso(dist,freq,family="normal", k=5,nlambda=50,tau=0.5,plot=FALSE, errorbars=FALSE)
	optCoef<-predict(nnlassoobj,mode="lambda",at=cvobj$lambda)
	cat(optCoef[2:length(optCoef)])
	cat("\n")
	return(cvobj)
    }, finally={})
    #print(cvobj$lambda)
    #optCoef<-predict(nnlassoobj,mode="lambda",at=cvobj$lambda)
    #print(optCoef[2:length(optCoef)],max.levels=(length(optCoef)-1),digits=3)
    #print(optCoef[2:length(optCoef)],digits=3,width=(length(optCoef)-1))
    #cat(optCoef[2:length(optCoef)])
    #cat("\n")
    #cross validation, get nnlassoobj, then predict
    #cv.nnlasso(x,y,family="normal",k=5)
}

#add loop structure to catch error recursively, until no error occur: xxx.recur.r

#dev.off()

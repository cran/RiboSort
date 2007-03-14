`speciesMDS` <-
function(x,dissimilarity="euclidean",type="non-metric",labels=TRUE)
{
distance=dist(x, method=dissimilarity)
if (min(distance)<=0)
{
distance=distance+abs(min(distance))+1
}

if(type=="classical")
{
mds=cmdscale(distance,k=2)
x1=mds[,1]
y1=mds[,2]
}
if(type=="sammon")
{
library(MASS)
mds=sammon(distance,k=2)
x1=mds$points[,1]
y1=mds$points[,2]
}
if(type=="non-metric")
{
library(MASS)
mds=isoMDS(distance,k=2)
x1=mds$points[,1]
y1=mds$points[,2]
}
p=plot(x1,y1, xlim=c(extendrange(x1,f=0.05)[1], extendrange(x1,f=0.1)[2]), xlab="Axis 1", ylab="Axis 2", type="p", pch=16,main=c(dissimilarity,type))
length=extendrange(x1,f=0.1)[2]-extendrange(x1,f=0.05)[1]
if (labels==TRUE)
{
text(x1+(length/20), y1, labels=rownames(x))
}
return(mds)
}


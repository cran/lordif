plotDIF <-
function(obj,labels=c("Reference","Focal"),cexp=0.8) {
sumpp<-function(pp) {
ws<-rowSums(pp*(col(pp)-1))
return(ws)
}
ndif<-sum(obj$flag)
maxcat<-ncol(obj$ipar.sparse)
if(ndif<1) stop("no DIF item present")
windows(record=T)
theta<-seq(obj$options$minTheta,obj$options$maxTheta,obj$options$inc)
difitems<-(1:obj$ni)[obj$flag]
itemnames<-row.names(obj$ipar.sparse)
gpar<-array(NA,c(ndif,maxcat,obj$ng))
cpar<-as.matrix(obj$ipar.sparse[1:(obj$ni-ndif),])
pp<-array(NA,c(length(theta),ndif,maxcat,obj$ng))
gtheta<-split(obj$calib.sparse$theta,obj$group)
gdensity<-matrix(0,length(theta),obj$ng)
for (i in 1:obj$ng) {
gdensity[,i]<-density(unlist(gtheta[names(table(obj$group))[i]]),n=length(theta),from=obj$options$minTheta,to=obj$options$maxTheta,bw=.25)$y
}
plot(theta,gdensity[,1],type="l",xlab="theta",ylab="Density",ylim=c(0,max(gdensity)),lty=1,col=1,main="Trait Distributions")
for (g in 2:obj$ng) {
lines(theta,gdensity[,g],lty=g,col=g)
}
legend("topright",labels,lty=1:obj$ng,col=1:obj$ng,cex=0.7,bg="white")
par(mfrow=c(2,2))
for (i in 1:length(difitems)) {
ncat<-obj$ncat[difitems[i]]
plot(theta,seq(0,ncat-1,along.with=theta),type="n",xlab="theta",ylab="Item Score",main=paste("Item True Score Functions - Item ",difitems[i],sep=""))
for (g in 1:obj$ng) {
gpar[i,,g]<-unlist(obj$ipar.sparse[which(itemnames==paste("I",difitems[i],".",g,sep="")),])
pp[,i,1:ncat,g]<-probgrm(theta,gpar[i,1,g],gpar[i,2:ncat,g])
lines(theta,sumpp(pp[,i,1:ncat,g]),lty=g,col=g)
}
legend("bottomright",labels,lty=1:obj$ng,col=1:obj$ng,cex=0.7,bg="white")
chi12<-paste(obj$stats[difitems[i],"df12"],")=",obj$stats[difitems[i],"chi12"],sep="")
pseudo12<-obj$stats[difitems[i],paste("pseudo12.",obj$options$pseudo.R2,sep="")]
beta12<-round(obj$stats[difitems[i],"beta12"],4)
chi13<-paste(obj$stats[difitems[i],"df13"],")=",obj$stats[difitems[i],"chi13"],sep="")
pseudo13<-obj$stats[difitems[i],paste("pseudo13.",obj$options$pseudo.R2,sep="")]
chi23<-paste(obj$stats[difitems[i],"df23"],")=",obj$stats[difitems[i],"chi23"],sep="")
pseudo23<-obj$stats[difitems[i],paste("pseudo23.",obj$options$pseudo.R2,sep="")]
text(min(theta),ncat-1,substitute(paste("Pr(",chi[12]^2,",",chi12,",",R[12]^2,"=",pseudo12,",",Delta,"(",beta[1],")=",beta12,sep="")),adj=c(0,1),cex=cexp)
text(min(theta),(ncat-1)*.9,substitute(paste("Pr(",chi[13]^2,",",chi13,",",R[13]^2,"=",pseudo13,sep="")),adj=c(0,1),cex=cexp)
text(min(theta),(ncat-1)*.8,substitute(paste("Pr(",chi[23]^2,",",chi23,",",R[23]^2,"=",pseudo23,sep="")),adj=c(0,1),cex=cexp)
plot(theta,seq(0,ncat-1,along.with=theta),type="n",xlab="theta",ylab="Item Score",main="Differences in Item True Score Functions")
for (g in 2:obj$ng) {
lines(theta,abs(sumpp(pp[,i,1:ncat,1])-sumpp(pp[,i,1:ncat,g])),lty=g,col=g)
}
plot(theta,seq(0,1,along.with=theta),type="n",xlab="theta",ylab="Probability",main="Item Response Functions")
for (g in 1:obj$ng) {
for (k in 1:ncat) {
lines(theta,pp[,i,k,g],lty=g,cex=0.1,col=g)
}
}
for (g in 1:obj$ng) {
text(obj$options$minTheta,.8-(g-1)*par()$cxy[2],paste(round(gpar[i,,g][!is.na(gpar[i,,g])],2),collapse=", "),col=g,adj=c(0,0),cex=cexp)
for (k in 2:ncat) {
if (!is.na(gpar[i,k,g])) text(gpar[i,k,g],0,"|",col=g)
}
}
plot(theta,seq(0,ncat-1,along.with=theta),type="n",xlab="theta",ylab="Size",main="Impact (Weighted by Density)")
for (g in 2:obj$ng) {
lines(theta,gdensity[,g]*abs(sumpp(pp[,i,1:ncat,1])-sumpp(pp[,i,1:ncat,g])),lty=g,col=g)
}
}
par(mfrow=c(1,2))
plot(theta,seq(0,sum(!is.na(obj$ipar))-obj$ni,along=theta),xlab="theta",ylab="TCC",type="n",main="All Items")
for (g in 1:obj$ng) {
apar<-rbind(cpar,gpar[,,g])
lines(theta,tcc(apar[,1],apar[,-1],theta),lty=g,col=g)
}
legend("bottomright",labels,lty=1:obj$ng,col=1:obj$ng,cex=cexp,bg="white")
plot(theta,seq(0,sum(!is.na(gpar[,,1]))-ndif,along=theta),xlab="theta",ylab="TCC",type="n",main="DIF Items")
for (g in 1:obj$ng) {
lines(theta,tcc(gpar[,1,g],matrix(gpar[,-1,g],nrow=ndif,byrow=T),theta),lty=g,col=g)
}
legend("bottomright",labels,lty=1:obj$ng,col=1:obj$ng,cex=cexp,bg="white")
layout(matrix(c(1,2),ncol=2),widths=c(1,2))
boxplot(obj$calib$theta-obj$calib.sparse$theta,col = "light grey")
difference<-obj$calib$theta-obj$calib.sparse$theta
plot(obj$calib$theta,difference,type="n",xlab="initial theta",ylab="initial - purified")
abline(h=0)
abline(h=mean(obj$calib$theta-obj$calib.sparse$theta),lty=2)
for (i in 1:obj$ng) {
points(obj$calib$theta[obj$group==as.numeric(names(table(obj$group))[i])],difference[obj$group==as.numeric(names(table(obj$group))[i])],col=i,pch=i)
}
legend("bottomright",labels,pch=1:obj$ng,col=1:obj$ng,cex=cexp,bg="white")
}


plotMC <-
function(obj,mfrow=c(3,1),cexp=1.0) {
nr<-obj$nr
Item<-1:dim(obj$cutoff)[1]
windows(record=T)
par(mfrow=mfrow)
par(mar=c(2,5,1,2)+0.1)
max.chi<-max(pretty(c(obj$cutoff$chi12,obj$cutoff$chi13,obj$cutoff$chi23)))
plot(Item,obj$cutoff$chi12,ylab=substitute(paste("Pr(",chi[12]^2,")",sep="")),ylim=c(0,max.chi),type="b")
abline(h=obj$alpha,col="red")
plot(Item,obj$cutoff$chi13,ylab=substitute(paste("Pr(",chi[13]^2,")",sep="")),ylim=c(0,max.chi),type="b")
abline(h=obj$alpha,col="red")
plot(Item,obj$cutoff$chi23,ylab=substitute(paste("Pr(",chi[23]^2,")",sep="")),ylim=c(0,max.chi),type="b")
abline(h=obj$alpha,col="red")
max.R2<-max(pretty(c(obj$cutoff$pseudo13.McFadden,obj$cutoff$pseudo13.Nagelkerke,obj$cutoff$pseudo13.CoxSnell)))
plot(Item,obj$cutoff$pseudo12.McFadden,ylab=substitute(paste(R[2]^2," - ",R[1]^2,sep="")),type="b",col="black",lty=1,ylim=c(0,max.R2))
points(Item,obj$cutoff$pseudo12.Nagelkerke,type="b",col="blue",lty=2,pch=2)
points(Item,obj$cutoff$pseudo12.CoxSnell,type="b",col="red",lty=3,pch=3)
legend("topleft",c("McFadden","Nagelkerke","Cox & Snell"),lty=1:3,col=c("black","blue","red"),pch=1:3,bg="white",cex=cexp)
plot(Item,obj$cutoff$pseudo13.McFadden,ylab=substitute(paste(R[3]^2," - ",R[1]^2,sep="")),type="b",col="black",lty=1,ylim=c(0,max.R2))
points(Item,obj$cutoff$pseudo13.Nagelkerke,type="b",col="blue",lty=2,pch=2)
points(Item,obj$cutoff$pseudo13.CoxSnell,type="b",col="red",lty=3,pch=3)
plot(Item,obj$cutoff$pseudo23.McFadden,ylab=substitute(paste(R[3]^2," - ",R[2]^2,sep="")),type="b",col="black",lty=1,ylim=c(0,max.R2))
points(Item,obj$cutoff$pseudo23.Nagelkerke,type="b",col="blue",lty=2,pch=2)
points(Item,obj$cutoff$pseudo23.CoxSnell,type="b",col="red",lty=3,pch=3)
max.beta<-max(pretty(obj$cutoff$beta12))
plot(Item,obj$cutoff$beta12,ylab=substitute(Delta(beta[1])),type="b",ylim=c(0,max.beta))
}


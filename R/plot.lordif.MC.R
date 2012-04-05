plot.lordif.MC <-
function(x,mfrow=c(3,1),...) {
if (class(x)!="lordif.MC") stop(paste(deparse(substitute(x)),"must be of class lordif.MC"))
nr<-x$nr
Item<-1:dim(x$cutoff)[1]
if (Sys.info()["sysname"]=="Windows") dev.new(record=T) else par(ask=T)
par(mfrow=mfrow)
par(mar=c(2,5,1,2)+0.1)
max.chi<-max(pretty(c(x$cutoff$chi12,x$cutoff$chi13,x$cutoff$chi23)))
plot(Item,x$cutoff$chi12,ylab=substitute(paste("Pr(",chi[12]^2,")",sep="")),ylim=c(0,max.chi),type="b",...)
abline(h=x$alpha,col="red")
plot(Item,x$cutoff$chi13,ylab=substitute(paste("Pr(",chi[13]^2,")",sep="")),ylim=c(0,max.chi),type="b",...)
abline(h=x$alpha,col="red")
plot(Item,x$cutoff$chi23,ylab=substitute(paste("Pr(",chi[23]^2,")",sep="")),ylim=c(0,max.chi),type="b",...)
abline(h=x$alpha,col="red")
max.R2<-max(pretty(c(x$cutoff$pseudo13.McFadden,x$cutoff$pseudo13.Nagelkerke,x$cutoff$pseudo13.CoxSnell)))
plot(Item,x$cutoff$pseudo12.McFadden,ylab=substitute(paste(R[2]^2," - ",R[1]^2,sep="")),type="b",col="black",lty=1,ylim=c(0,max.R2))
points(Item,x$cutoff$pseudo12.Nagelkerke,type="b",col="blue",lty=2,pch=2)
points(Item,x$cutoff$pseudo12.CoxSnell,type="b",col="red",lty=3,pch=3)
legend("topleft",c("McFadden","Nagelkerke","Cox & Snell"),lty=1:3,col=c("black","blue","red"),pch=1:3,bg="white")
plot(Item,x$cutoff$pseudo13.McFadden,ylab=substitute(paste(R[3]^2," - ",R[1]^2,sep="")),type="b",col="black",lty=1,ylim=c(0,max.R2))
points(Item,x$cutoff$pseudo13.Nagelkerke,type="b",col="blue",lty=2,pch=2)
points(Item,x$cutoff$pseudo13.CoxSnell,type="b",col="red",lty=3,pch=3)
plot(Item,x$cutoff$pseudo23.McFadden,ylab=substitute(paste(R[3]^2," - ",R[2]^2,sep="")),type="b",col="black",lty=1,ylim=c(0,max.R2))
points(Item,x$cutoff$pseudo23.Nagelkerke,type="b",col="blue",lty=2,pch=2)
points(Item,x$cutoff$pseudo23.CoxSnell,type="b",col="red",lty=3,pch=3)
max.beta<-max(pretty(x$cutoff$beta12))
plot(Item,x$cutoff$beta12,ylab=substitute(Delta(beta[1])),type="b",ylim=c(0,max.beta),...)
}

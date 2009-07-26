`extract` <-
function(ipar) {
ncat<-unlist(lapply(ipar$coefficients,length))
maxCat<-max(ncat)
ni<-length(ipar$coefficients)
out.ipar<-data.frame(matrix(NA,ni,maxCat))
names(out.ipar)<-c("a",paste("cb",1:(maxCat-1),sep=""))
for (i in 1:ni) {
out.ipar[i,1]<-ipar$coefficients[[i]][ncat[i]]
for (j in 1:(ncat[i]-1)) {
out.ipar[i,j+1]<-ipar$coefficients[[i]][j]/out.ipar[i,1]
}
}
return(out.ipar)
}


calcprob <-
function(ipar,theta){
ni<-nrow(ipar) 
maxCat<-ncol(ipar) 
NCAT<-apply(ipar,1,function (x) sum(!is.na(x))) 
DISC<-ipar[,"a"] 
CB<-ipar[,paste("cb",1:(maxCat-1),sep=""),drop=F] 
pp<-array(NA,c(length(theta),ni,maxCat)) 
for (i in 1:ni) {
pp[,i,1:(NCAT[i])]<-probgrm(theta,DISC[i],CB[i,])
}
return(pp)
}


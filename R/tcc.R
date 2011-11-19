tcc <-
function(a,cb,theta) {
ni<-length(a)
T<-numeric(length(theta))
for (i in 1:ni) {
T<-T+probgrm(theta,a[i],cb[i,])%*%seq(0,sum(!is.na(cb[i,])));
}
return(T)
}


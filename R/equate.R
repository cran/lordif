equate <-
function(ipar.to,ipar.from,theta) {
SL<-function(AK) {
a.to<-ipar.to$a
a.from<-ipar.from$a
cb.to<-ipar.to[-1]
cb.from<-ipar.from[-1]
a.from<-a.from/AK[1]
cb.from<-cb.from*AK[1]+AK[2]
return(sum(tcc(a.to,cb.to,theta)-tcc(a.from,cb.from,theta))^2)
}
return (nlminb(c(1.0,0.0),SL,lower=c(0.5,-2.0),upper=c(2.0,2.0))$par)
}

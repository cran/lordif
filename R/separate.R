`separate` <-
function(resp,flag,gr) {
resp.nodif<-resp[,!flag] 
resp.dif<-resp[,flag] 
nobs<-length(gr) 
gr.freq<-table(gr) 
ng<-length(gr.freq) 
gr.label<-names(gr.freq) 
ndif<-ncol(resp.dif) 
dif.items<-names(resp.dif) 
sparse.resp<-data.frame(matrix(NA,nobs,ncol(resp.dif)*ng)) 
colnames(sparse.resp)<-paste(rep(dif.items,rep(ng,ndif)),".",rep(1:ng,ndif),sep="") 
for (i in 1:ndif) {
for (j in 1:ng) {
sparse.resp[gr==gr.label[j],ng*i-ng+j]<-resp.dif[gr==gr.label[j],i] 
}
}
out<-data.frame(resp.nodif,sparse.resp)
return(out)
}


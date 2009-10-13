lordif <-
function(resp.data,group,selection=NULL,criterion="Chisqr",pseudo.R2="McFadden",alpha=0.01,beta.change=0.1,R2.change=0.02,
              maxIter=10,minCell=5,minTheta=-4.0,maxTheta=4.0,inc=0.1) {
if (!(toupper(criterion) %in% c("CHISQR","R2","BETA"))){warning("criterion must be one of the following: \"Chisqr\", \"R2\", or \"Beta\", will be reset to default");criterion<-"Chisqr"}
if (!(toupper(pseudo.R2) %in% c("MCFADDEN","NAGELKERKE","COXSNELL"))) {warning("pseudo.R2 must be one of the following \"McFadden\", \"Nagelkerke\", or \"CoxSnell\", will be reset to default");pseudo.R2<-"McFadden"}
if (alpha<=0 & alpha>1) {warning ("alpha must be > 0 & < 1, will be reset to default");alpha<-.01}
if (beta.change<=0 & beta.change>=1) {warning ("beta.change must be > 0 & < 1, will be reset to default");beta.change<-.1};
if (R2.change<=0 & R2.change>=1) {warning("R2.change must be > 0 & < 1, will be reset to default");R2.change<-.035}
if (maxIter<1) {warning("maxIter must be >= 1, will be reset to default");maxIter<5}
if (minCell<1) {warning("minCell must be >= 1, will be reset to 5");minCell<-5}
if (minTheta>=maxTheta) {warning("minTheta must be < maxTheta, will be reset to default");minTheta<--4;maxTheta<-4}
if (inc<=0) {warning("inc must be > 0, will be reset to default");inc<-.1}
resp.data<-resp.data[!is.na(group),]
group<-group[!is.na(group)]
nobs<-nrow(resp.data)
if (length(selection)==0) selection<-1:ncol(resp.data) 
tni<-ncol(resp.data)
ni<-length(selection)
if (nobs != length(group)) stop ("nobs in response matrix and group vector are not congruent") 
if (ni<5) stop("number of items must be at least 5") 
if (prod(is.element(selection,1:tni))!=1) stop("selection is not in the total set") 
compare<-function(x,table) {
for (i in 1:nrow(table)) {
if (all(table[i,]==x)) return(TRUE)
}
return(FALSE)
}
resp.recoded<-data.frame(matrix(NA,nobs,ni)) 
names(resp.recoded)<-paste("I",selection,sep="") 
for (i in 1:ni) {
resp.recoded[,i]<-collapse(resp.data[,selection[i]],group,minCell) 
}
ncat<-as.numeric(apply(resp.recoded,2,max,na.rm=T))
ng<-length(table(group))
calib<-grm(resp.recoded,constrained=FALSE,IRT.param=TRUE,control=list(iter.qN=50))
ipar<-extract(calib)
theta.grid<-seq(minTheta,maxTheta,inc)
theta<-calctheta(ipar,resp.recoded,theta.grid)
meanraw<-apply(resp.recoded,1,mean,na.rm=T)
outraw<-rundif(selection,resp.recoded,meanraw,group,criterion,alpha,beta.change,pseudo.R2,R2.change)
out<-rundif(selection,resp.recoded,theta$theta,group,criterion,alpha,beta.change,pseudo.R2,R2.change)
item.labels<-names(resp.recoded)
row.names(ipar)<-item.labels 
flag.matrix<-rbind(logical(ni),!logical(ni))
flags<-out$flag 
iter<-1
ndif<-sum(flags)
cat(paste("Iteration ",iter," : ",ndif," items flagged for DIF (",paste(selection[flags],sep="",collapse=","),")\n",sep=""))
if (ndif==ni) {
warning("all items got flagged for DIF - stopping\n")
}
sparse.matrix<-NULL;calib.sparse<-NULL;ipar.sparse<-NULL;theta.sparse<-NULL
if (!compare(flags,flag.matrix)) {
repeat {
iter<-iter+1
flag.matrix<-rbind(flag.matrix,flags)
sparse.matrix<-separate(resp.recoded,flags,group)
calib.sparse<-grm(sparse.matrix,constrained=FALSE,IRT.param=TRUE,control=list(iter.qN=50)) 
ipar.sparse<-extract(calib.sparse) 
eqconst<-equate(ipar[!flags,],ipar.sparse[1:sum(!flags),],theta.grid)
ipar.sparse[,1]<-ipar.sparse[,1]/eqconst[1]
ipar.sparse[,2:ncol(ipar.sparse)]<-ipar.sparse[,2:ncol(ipar.sparse)]*eqconst[1]+eqconst[2]
theta.sparse<-calctheta(ipar.sparse,sparse.matrix,theta.grid)
pre.flags<-flags
out<-rundif(selection,resp.recoded,theta.sparse$theta,group,criterion,alpha,beta.change,pseudo.R2,R2.change)
flags<-out$flag
ndif<-sum(flags)
cat(paste("Iteration ",iter," : ",ndif," items flagged for DIF (",paste(selection[flags],sep="",collapse=","),")\n",sep=""))
if (ndif==ni) {
warning("all items got flagged for DIF - stopping\n")
break
}
if (compare(flags,flag.matrix) | iter == maxIter) {
if (!all(pre.flags==flags)) {
sparse.matrix<-separate(resp.recoded,flags,group)
calib.sparse<-grm(sparse.matrix,constrained=FALSE,IRT.param=TRUE,control=list(iter.qN=50)) 
ipar.sparse<-extract(calib.sparse) 
eqconst<-equate(ipar[!flags,],ipar.sparse[1:sum(!flags),],theta.grid)
ipar.sparse[,1]<-ipar.sparse[,1]/eqconst[1]
ipar.sparse[,2:ncol(ipar.sparse)]<-ipar.sparse[,2:ncol(ipar.sparse)]*eqconst[1]+eqconst[2]
theta.sparse<-calctheta(ipar.sparse,sparse.matrix,theta.grid)
}
break
}
}
if (!compare(flags,flag.matrix) & (iter == maxIter)) {
sparse.matrix<-separate(resp.recoded,flags,group)
calib.sparse<-grm(sparse.matrix,constrained=FALSE,IRT.param=TRUE,control=list(iter.qN=50)) 
ipar.sparse<-extract(calib.sparse)
eqconst<-equate(ipar[!flags,],ipar.sparse[1:sum(!flags),],theta.grid)
ipar.sparse[,1]<-ipar.sparse[,1]/eqconst[1]
ipar.sparse[,2:ncol(ipar.sparse)]<-ipar.sparse[,2:ncol(ipar.sparse)]*eqconst[1]+eqconst[2]
theta.sparse<-calctheta(ipar.sparse,sparse.matrix,theta.grid) 
}
row.names(ipar.sparse)<-names(sparse.matrix)
}
return(list(options=list(criterion=criterion,pseudo.R2=pseudo.R2,alpha=alpha,beta.change=beta.change,R2.change=R2.change,
              maxIter=maxIter,minCell=minCell,minTheta=minTheta,maxTheta=maxTheta,inc=inc),selection=selection,stats=out$stats,
              flag=out$flag,recoded=resp.recoded,group=group,ng=ng,ni=ni,ncat=ncat,calib=theta,calib.sparse=theta.sparse,
              iteration=iter,ipar=ipar,ipar.sparse=ipar.sparse,stats.raw=outraw$stats,meanraw=meanraw,flag.raw=outraw$flag))
}


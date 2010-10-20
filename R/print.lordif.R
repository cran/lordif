print.lordif <-
function(x, ...) {
cat("Call:\n")
print(x$call)
cat("\n")
cat(paste("  Number of DIF groups:",x$ng,"\n\n"))
cat(paste("  Number of items flagged for DIF:",sum(x$flag),"of",x$ni,"\n\n"))
cat(paste("  Items flagged:",paste(which(x$flag),collapse=", "),"\n\n"))
cat(paste("  Number of iterations for purification:",x$iteration,"of",x$options$maxIter,"\n\n"))
cat(paste("  Detection criterion:",x$options$criterion,"\n\n"))
threshold<-switch(toupper(x$options$criterion),
"CHISQR"=paste("alpha =",x$options$alpha),
"R2"=paste("R-square change >=",x$options$R2.change),
"BETA"=paste("Beta change >=",x$options$beta.change))
cat(paste("  Threshold:",threshold,"\n\n"))
if (sum(x$flag)>0) {
out<-switch(toupper(x$options$criterion),
"CHISQR"=c("chi12","chi13","chi23"),
            "R2"=switch(toupper(x$options$pseudo.R2),
"MCFADDEN"=c("pseudo12.McFadden","pseudo13.McFadden","pseudo23.McFadden"),
"NAGELKERKE"=c("pseudo12.Nagelkerke","pseudo13.Nagelkerke","pseudo23.Nagelkerke"),
"COXSNELL"=c("pseudo12.CoxSnell","pseudo13.CoxSnell","pseudo23.CoxSnell")),
"BETA"="beta12")

print(x$stats[c("item","ncat",out)])
}
invisible(x)
}


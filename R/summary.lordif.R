summary.lordif <-
function(object, ...) {
cat("Call:\n")
print(object$call)
cat("\n")
print(object$options[c("criterion","alpha","pseudo.R2","R2.change","beta.change","maxIter","minCell")])
print(object[c("stats","flag","flag.raw")])
invisible(object)
}

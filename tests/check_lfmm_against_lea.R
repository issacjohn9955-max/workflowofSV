# Verify the batched second-stage regression and GC against LEA on the demo.
arg <- commandArgs(trailingOnly=TRUE)
stopifnot(length(arg)==1L)
root <- arg[1]
if (nzchar(Sys.getenv("SV_R_LIBRARY"))) .libPaths(c(Sys.getenv("SV_R_LIBRARY"),.libPaths()))
p <- file.path(root,"prepared")
g <- as.matrix(read.delim(file.path(p,"lfmm_dosage.tsv.gz"),row.names=1,check.names=FALSE))
d <- read.delim(file.path(p,"individual_design.tsv"),row.names=1)
truth <- jsonlite::fromJSON(file.path(p,"settings.json"))
env <- truth$environments
for (k in c(3,5,10)) {
    x <- scale(as.matrix(d[,c(env,paste0("PC",1:k),"lon","lat")]))
    fit <- LEA::lfmm2(input=g,env=x,K=truth$lfmm_k,lambda=1e-5,effect.sizes=FALSE)
    native <- LEA::lfmm2.test(fit,input=g,env=x,linear=TRUE,genomic.control=TRUE)
    tbl <- read.delim(file.path(root,paste0("lfmm_pc",k),"associations.tsv.gz"))
    stopifnot(identical(tbl$locus,rep(colnames(g),length(env))))
    a <- as.vector(t(native$pvalues[seq_along(env),,drop=FALSE]))
    b <- as.vector(t(native$zscores[seq_along(env),,drop=FALSE]))
    stopifnot(max(abs(a-tbl$p_gc))<1e-9, max(abs(b-tbl$statistic))<1e-8)
}
cat("PASS: LFMM2 statistics and genomic-control P agree with LEA for PC3, PC5 and PC10.\n")

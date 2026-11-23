# Each process fits one method and one PC count on frozen matrices.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==4L, args[1] %in% c("lfmm","rda"))
method <- args[1]
pc <- as.integer(args[2])
prepared <- args[3]
out <- args[4]
stopifnot(pc %in% c(3L,5L,10L))
lib <- Sys.getenv("SV_R_LIBRARY")
if (nzchar(lib)) .libPaths(c(lib,.libPaths()))
cfg <- jsonlite::fromJSON(file.path(prepared,"settings.json"))
write_tsv <- function(data, path) {
    p <- file.path(out,path)
    connection <- if (endsWith(path,".gz")) gzfile(p,"wt") else file(p,"wt")
    on.exit(close(connection))
    write.table(data,connection,sep="\t",quote=FALSE,row.names=FALSE,na="NA")
}
read_frozen <- function(prefix, matrix_name) {
    classes <- if (cfg$rda_blocked) c(id="character",block="character") else c(id="character")
    meta <- read.delim(file.path(prepared,paste0(prefix,"_design.tsv")),row.names=1,check.names=FALSE,
                      colClasses=classes)
    y <- as.matrix(read.delim(file.path(prepared,matrix_name),row.names=1,check.names=FALSE,
                             colClasses=c(id="character")))
    stopifnot(identical(rownames(y),rownames(meta)),!anyDuplicated(colnames(y)),all(is.finite(y)))
    list(y=y,meta=meta)
}
env_names <- cfg$environments
pc_names <- paste0("PC",seq_len(pc))
geo <- if (cfg$adjust_geography) c("lon","lat") else character(0)
cov_names <- c(pc_names,geo)
design <- function(meta) {
    x <- scale(as.matrix(meta[,c(env_names,cov_names),drop=FALSE]))
    if (!all(is.finite(x)) || qr(cbind(1,x))$rank!=1+ncol(x))
        stop("Design contains a constant or redundant environment/covariate.")
    x
}

if (method=="lfmm") {
    d <- read_frozen("individual","lfmm_dosage.tsv.gz")
    x <- design(d$meta)
    g <- d$y
    fit <- LEA::lfmm2(input=g,env=x,K=cfg$lfmm_k,lambda=1e-5,effect.sizes=FALSE)
    model <- cbind(Intercept=1,x,fit@U)
    qr_model <- qr(model)
    df <- nrow(model)-ncol(model)
    if (qr_model$rank != ncol(model) || df <= 0) stop("LFMM2 test design is rank deficient or saturated.")
    stderr_unit <- numeric(ncol(model))
    stderr_unit[qr_model$pivot] <- sqrt(diag(chol2inv(qr.R(qr_model))))
    m <- ncol(g)
    indices <- seq_along(env_names)+1L
    coefficients <- statistics <- matrix(NA_real_,length(env_names),m)
    for (first in seq.int(1L,m,by=2048L)) {
        j <- first:min(first+2047L,m)
        y <- g[,j,drop=FALSE]
        beta <- qr.coef(qr_model,y)
        rss <- colSums(qr.resid(qr_model,y)^2)
        tss <- colSums(sweep(y,2,colMeans(y))^2)
        ok <- is.finite(rss)&rss>1e-12*tss
        if (any(ok)) {
            coefficients[,j[ok]] <- beta[indices,ok,drop=FALSE]
            errors <- outer(stderr_unit[indices],sqrt(rss[ok]/df))
            statistics[,j[ok]] <- beta[indices,ok,drop=FALSE]/errors
        }
    }
    results <- diagnostics <- vector("list",length(env_names))
    for (i in seq_along(env_names)) {
        stat <- statistics[i,]
        good <- is.finite(stat)
        inflation <- if (any(good)) median(stat[good]^2)/qchisq(.5,1) else NA_real_
        p_gc <- rep(NA_real_,m)
        reasons <- ifelse(good,"tested","degenerate_residual")
        if (is.finite(inflation) && inflation>0) {
            p_gc[good] <- pchisq(stat[good]^2/inflation,1,lower.tail=FALSE)
        } else {
            reasons[good] <- "invalid_genomic_control"
        }
        q <- p.adjust(p_gc,"BH")
        results[[i]] <- data.frame(locus=colnames(g),environment=env_names[i],pc_count=pc,
            latent_k=cfg$lfmm_k,beta=coefficients[i,],statistic=stat,p_ols=2*pt(-abs(stat),df),
            p_gc=p_gc,q=q,status=reasons)
        diagnostics[[i]] <- data.frame(environment=env_names[i],pc_count=pc,latent_k=cfg$lfmm_k,
            n=nrow(g),residual_df=df,loci=m,tested=sum(is.finite(q)),gc_lambda=inflation,
            significant=sum(q<cfg$alpha,na.rm=TRUE))
    }
    results <- do.call(rbind,results)
    write_tsv(results,"associations.tsv.gz")
    write_tsv(results[which(results$q<cfg$alpha),],"discoveries.tsv")
    write_tsv(do.call(rbind,diagnostics),"diagnostics.tsv")
} else {
    d <- read_frozen("group","rda_frequency.tsv.gz")
    x <- as.data.frame(design(d$meta))
    Y <- d$y
    form <- as.formula(paste("Y ~",paste(env_names,collapse=" + "),
                             "+ Condition(",paste(cov_names,collapse=" + "),")"))
    fit <- vegan::rda(form,data=x,scale=TRUE)
    df <- nrow(x)-1-ncol(x)
    stopifnot(df>0,fit$CCA$rank>0L)
    set.seed(cfg$seed)
    restriction <- if (cfg$rda_blocked) factor(d$meta$block) else NULL
    control <- permute::how(nperm=cfg$permutations,blocks=restriction)
    permutation_matrix <- permute::shuffleSet(nrow(Y),nset=cfg$permutations,control=control)
    if (nrow(permutation_matrix)<19L) stop("Too few valid group permutations. Review the block design.")
    if (cfg$rda_blocked) {
        for (i in seq_len(nrow(permutation_matrix)))
            stopifnot(all(restriction[permutation_matrix[i,]]==restriction))
    }
    global <- anova(fit,model="reduced",permutations=permutation_matrix,parallel=1)
    marginal <- anova(fit,by="margin",model="reduced",permutations=permutation_matrix,parallel=1)
    by_axis <- anova(fit,by="axis",model="reduced",permutations=permutation_matrix,parallel=1)
    axes <- names(fit$CCA$eig)
    p_axis <- by_axis[match(axes,rownames(by_axis)),"Pr(>F)"]
    q_axis <- p.adjust(p_axis,"BH")
    p_global <- global[1,"Pr(>F)"]
    eligible <- is.finite(q_axis)&q_axis<cfg$alpha&p_global<cfg$alpha
    r2 <- vegan::RsquareAdj(fit)
    write_tsv(data.frame(pc_count=pc,groups=nrow(Y),loci=ncol(Y),residual_df=df,
        permutations=nrow(permutation_matrix),blocked=cfg$rda_blocked,global_p=p_global,
        pseudo_f=global[1,"F"],constrained_fraction=fit$CCA$tot.chi/fit$tot.chi,
        adjusted_r2=r2$adj.r.squared),"global.tsv")
    write_tsv(data.frame(axis=axes,inertia=as.numeric(fit$CCA$eig),p=p_axis,q=q_axis,eligible=eligible),
              "axes.tsv")
    write_tsv(data.frame(environment=env_names,p=marginal[match(env_names,rownames(marginal)),"Pr(>F)"]),
              "marginal_terms.tsv")
    choices <- seq_along(axes)
    scores <- function(which) as.matrix(vegan::scores(fit,display=which,choices=choices,scaling=2))
    loadings <- scores("species")
    loci <- rownames(loadings)
    z <- scale(loadings)
    peak <- if (any(eligible)) apply(abs(z[,eligible,drop=FALSE]),1,max) else rep(NA_real_,length(loci))
    result <- data.frame(locus=loci,max_abs_loading_z=peak,candidate=is.finite(peak)&peak>cfg$loading_z)
    write_tsv(cbind(result,as.data.frame(loadings)),"locus_loadings.tsv.gz")
    write_tsv(result[result$candidate,],"exploratory_candidates.tsv")
    si <- scores("sites")
    bp <- scores("bp")
    write_tsv(data.frame(group=rownames(si),si),"site_scores.tsv")
    write_tsv(data.frame(term=rownames(bp),bp),"environment_arrows.tsv")
    write_tsv(data.frame(permutation_matrix),"permutations.tsv.gz")
}
writeLines(capture.output(sessionInfo()),file.path(out,"session_info.txt"))

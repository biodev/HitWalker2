library(Matrix)

#Rscript --vanilla network/perform_rwr.R /var/www/hitwalker2_inst/static/network/data/9606.protein.links.v9.1.mm 0.4 9606.ENSP00000373700 9606.ENSP00000381066 9606.ENSP00000335153

cargs <- commandArgs(TRUE)

mm.file.base <- cargs[1]
graph.thresh <- as.numeric(cargs[2])
hit.prots <- cargs[3:length(cargs)]

init.mat <- readMM(paste0(mm.file.base, ".mtx"))

which.zero <- init.mat@x < graph.thresh

cur.mat <- new("dsTMatrix",  i = init.mat@i[which.zero==F], j = init.mat@j[which.zero==F], x = init.mat@x[which.zero==F], Dim=init.mat@Dim, uplo=init.mat@uplo)

new.graph.names <- read.delim(paste0(mm.file.base, ".names"), header=FALSE, stringsAsFactors=FALSE)

rownames(cur.mat) <- new.graph.names[,1]
colnames(cur.mat) <- new.graph.names[,1]

d.sub.mat <- Diagonal(x=1/sqrt(rowSums(cur.mat)))
    
graph.sp.mat <- d.sub.mat %*% cur.mat %*% d.sub.mat

seed.prots <- rep(1, length(hit.prots))
names(seed.prots) <- hit.prots

residue <- 1
iter <- 1

#c=.3, threshold=1e-10, maxit=100
c.val <- .3
threshold <- 1e-10
maxit <- 100
verbose <- T

#probability of restart...
prox.vector <- rep(0, nrow(graph.sp.mat))
names(prox.vector) <- rownames(graph.sp.mat)

seed.genes.in.graph <- intersect(names(seed.prots), names(prox.vector))

if (length(seed.genes.in.graph) == 0)
{
    default.fail("ERROR: None of the seeds were found in the supplied graph")
}

if (verbose == TRUE) cat (paste("Using", length(seed.genes.in.graph), "seed genes\n"))

#make sure the supplied values sum to 1

if(sum(seed.prots[seed.genes.in.graph]) != 1)
{
    if (verbose == TRUE) cat("Sum of selected seed.prots not equal to one, renormalizing...\n")
    seed.prots[seed.genes.in.graph] <- seed.prots[seed.genes.in.graph]/sum(seed.prots[seed.genes.in.graph])
}

prox.vector[seed.genes.in.graph] <- seed.prots[seed.genes.in.graph] #if the weights are being uniformly applied... 1/length(seed.genes.in.graph)

restart.vector <- prox.vector

while (residue > threshold && iter < maxit)
{
    if (verbose == TRUE) cat (paste("iter:", iter, "residue:", residue, "\n"))
    old.prox.vector <-  prox.vector
    prox.vector <- as.numeric((1-c.val)*graph.sp.mat %*% prox.vector + (c.val*restart.vector))
    residue <- as.numeric(dist(rbind(prox.vector, old.prox.vector), method="euclidean"))
    iter <- iter + 1; 
}

if (iter == maxit)
{
    stop(paste("Seed prots:", paste(seed.prots, collapse=","), "did not converge in 100 iterations..."))
}

names(prox.vector) <- rownames(graph.sp.mat)

ord.prox.vector <- sort(prox.vector, decreasing=TRUE)

write.table(ord.prox.vector, file="temp_ranking_results.txt", quote=F, col.names=F)
    
    
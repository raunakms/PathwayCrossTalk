# -------------------------------------------------------------------- #
# THE ANALYSIS IS BASED ON THE FOLLOWING PAPER
#
# CITATION: 
# Yong Li, Pankaj Agarwal, Dilip Rajagopalan, 
# A global pathway crosstalk network, 
# Bioinformatics, Volume 24, Issue 12, 15 June 2008, Pages 1442–1447, 
# https://doi.org/10.1093/bioinformatics/btn200
# -------------------------------------------------------------------- #

### DEFINE LIBRARIES ----
library("stringr")
library("igraph")
library("foreach")
library("doParallel")
library("reshape2")

### DEFINE PARAMETERS ---
n.rmd <- 100 #NO. OF RANDOMIZATION STEPS
batch <- "KEGG"

### SET PATHS ---
dir.wrk <- get.wd()
dir.data <- file.path(dir.wrk, "data")
dir.output <- file.path(dir.wrk, "output")
dir.script <- file.path(dir.wrk, "script")

### CREATE ENVIRONMENT ---
dir.batch <- file.path(dir.output, batch)
dir.rmd <- file.path(dir.batch, "random_set")
dir.plot <- file.path(dir.batch, "plot")

cat(paste(Sys.time()), "DIRECTORY STRUCTURE : CREATING ...","\n", sep="\t")
	dir.create(dir.batch, showWarnings=FALSE)
	dir.create(dir.rmd, showWarnings=FALSE)
	dir.create(dir.plot, showWarnings=FALSE)
cat(paste(Sys.time()), "DIRECTORY STRUCTURE : CREATED ...","\n", sep="\t")

### DEFINE FILES ---
file.network <- file.path(dir.data, "string10_network.tsv")
file.genesets <- file.path(dir.data, "c2.cp.kegg.v5.0.symbols.txt")

### LOAD NETWORK GRAPH ---
dat.network <- read.delim(file.network, header=F, stringsAsFactors=F, na.strings = "null")

### CREATE iGRAPH OBJECT ---
g <- graph.data.frame(dat.network, directed=FALSE)
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE)

### LOAD PATHWAY DATA ---
dat.genesets <- read.delim(file.genesets, header=T, stringsAsFactors=F)

### 1: TRIM PATHWAY DATA -------------------------------------------------------------------------------
# Generate a set of pathways for crosstalk analysis by removing
# 	pathways containing less than 5 genes or more than 200 genes.
# 	Pathways with too many genes might be too generic
# 	and pathways with too few genes may not have sufficient
# 	biological content. These size cutoffs were set up arbitrarily.

dat.genesets$lgenes <- unlist(lapply(str_split(dat.genesets$Genesets, ","), function(x) length(x)))
del.index <- c(which(dat.genesets$lgenes < 5), which(dat.genesets$lgenes > 200))
dat.genesets <- dat.genesets[-del.index, ]

### 2: EVALUATE GENE-OVERLAP BETWEEN PATHWAYS ------------------------------------------------------------
# Evaluate gene overlap between any given pair of pathways by
# 	performing Fisher Exact test (Al Shahrour et al., 2004). Raw
# 	P-values were adjusted by false discovery rate (FDR) Benjamini-
# 	Hochberg (BH) procedure (Benjamini and Hochberg, 1995) to
# 	account for multiple hypothesis testing. Any pathway pairs with
# 	adjusted P-value 0.05 were considered to overlap significantly
# 	and were removed during the pruning step (see below).

genes.U <- unique(unlist(str_split(dat.genesets$Genesets, ",")))

### Declate Cluster 
no_cores <- 60
cl <- makeCluster(no_cores)
registerDoParallel(cl)

lpar <- foreach(i = 1:nrow(dat.genesets), .combine='cbind') %:%
			foreach(j = 1:nrow(dat.genesets), .combine='c', .packages="stringr") %dopar% {
				if(i > j){
					genes.A <- str_split(dat.genesets$Genesets[i], ",")[[1]]
					genes.B <- str_split(dat.genesets$Genesets[j], ",")[[1]]
				
					contingency.table <- matrix(0, nrow=2, ncol=2, dimnames=list(c("YES","NO"), c("YES","NO")))
					contingency.table[1,1] <- length(intersect(genes.A, genes.B))
					contingency.table[1,2] <- length(setdiff(genes.A, genes.B))
					contingency.table[2,1] <- length(setdiff(genes.B, genes.A))
					contingency.table[2,2] <- length(setdiff(genes.U, unique(c(genes.A, genes.B))))
		
					ftest <- fisher.test(contingency.table)
					lpar <- ftest$p.value
				} else{
					lpar <- NA
				}
			}

stopCluster(cl)
colnames(lpar) <- rownames(lpar) <- dat.genesets$Category

dat.overlap <- melt(lpar)
colnames(dat.overlap) <- c("PathwayA","PathwayB", "pvalue")
dat.overlap <- dat.overlap[-which(is.na(dat.overlap$pvalue)),]
dat.overlap$fdr <- p.adjust(dat.overlap$pvalue, method="BH")
dat.overlap$overlap.sig <- ifelse(dat.overlap$fdr < 0.05, 1, 0)     ## we are interested in '0's i.e. significant non-overlaps

dat.overlap$PathwayA <- as.character(dat.overlap$PathwayA)
dat.overlap$PathwayB <- as.character(dat.overlap$PathwayB)

write.table(dat.overlap, file.path(dir.batch, paste("01_dat_pathway_overlap_", batch, ".tsv", sep="")), sep="\t", row.names=F, col.names=T, quote=F)


### 3: EVALUATE PROTEIN INTERACTIONS BETWEEN PATHWAYS ------------------------------------------------------
# Count number of protein interactions between any two pathways.
# 	For each pathway pair, first remove genes common to both
# 	pathways, then count all protein interactions between these two
# 	pathways.

adj <- as.matrix(as_adjacency_matrix(g, type="both", names=T))
lpar <- compute_protein_interactions(dat=dat.genesets, adj)
lpar[is.na(lpar)] <- 0
write.table(lpar, file.path(dir.batch, paste("02_dat_pathway_interaction_matrix_", batch, ".tsv", sep="")), sep="\t", row.names=T, col.names=NA, quote=F)

### 4: ESTIMATE BACKGROUND DISTRIBUTION OF PROTEIN INTERACTIONS BETWEEN PATHWAYS ---------------------------
# Estimate background distribution of protein interaction count of
# 	each pathway pair. Each pathway was randomized as follows. Go
# 	through all genes in a given pathway. If a gene does not have any
# 	interactions, skip it. If a gene has interactions, first count the
# 	number of genes it interacts with, then randomly draw a gene
# 	from the protein interaction dataset which interacts with the
# 	same or similar number of genes, and replace the original
# 	pathway gene with this newly selected gene. Once both pathways
# 	were randomized, Step 3 is performed to count the number of
# 	interactions between them. This randomization step is repeated
# 	1000 times.

deg <- degree(g, c(1:length(V(g))), loops=F, normalized=F)
deg <- data.frame(Gene=names(deg), Degree=as.numeric(deg))
deg$Gene <- as.character(deg$Gene)
deg <- deg[order(deg$Degree, decreasing=T),]
deg$Rank <- c(1:nrow(deg))


##### RANDOMIZE GENESETS -------------
### Declate Cluster 
no_cores <- 60
cl <- makeCluster(no_cores)
registerDoParallel(cl)

lpar <- foreach(k = 1:n.rmd, .combine='rbind') %:%
			foreach(i = 1:nrow(dat.genesets), .combine='rbind', .packages="stringr") %dopar% {
				genes.A <- str_split(dat.genesets$Genesets[i], ",")[[1]]
				genes.A <- intersect(genes.A, deg$Gene)
				genes.rmd <- rep("", length(genes.A))
				
				for(j in 1:length(genes.A)){
					gene0 <- genes.A[j]
					grank <- deg$Rank[which(deg$Gene == gene0)]
		
					if(grank < 10){
						genes.set <- setdiff(deg$Gene[which(deg$Rank %in%  c(1:(grank+10)))], gene0)
					}else if(grank > nrow(deg)-10){
						genes.set <- setdiff(deg$Gene[which(deg$Rank %in%  c((grank-10):nrow(deg)))], gene0)
					}else{
						genes.set <- setdiff(deg$Gene[which(deg$Rank %in% c((grank-10):(grank+10)))], gene0)
					}
		
					genes.rmd[j] <- sample(genes.set, 1)
				}
				lpar <- data.frame(Itr=k, Category=dat.genesets$Category[i], Genesets=paste(genes.rmd, collapse=","))
			}
stopCluster(cl)

dat.rmd.genesets <- lpar
dat.rmd.genesets$Category <- as.character(dat.rmd.genesets$Category)
dat.rmd.genesets$Genesets <- as.character(dat.rmd.genesets$Genesets)

write.table(dat.rmd.genesets, file.path(dir.batch, paste("03_dat_pathway_interaction_matrix_random_", batch, ".tsv", sep="")), sep="\t", row.names=F, col.names=T, quote=F)


adj <- as.matrix(as_adjacency_matrix(g, type="both", names=T))

for(k in 1:n.rmd){
	cat(paste(Sys.time()), "ITERATION START:",k, "\n", sep="\t")
	
	dat.g <- dat.rmd.genesets[which(dat.rmd.genesets$Itr == k),]
	mat <- compute_protein_interactions(dat=dat.g, adj)
	
	file.output <- file.path(dir.rmd, paste("random_genesets_interaction_count_",k,".tsv", sep=""))
	write.table(mat, file.output, sep="\t", row.names=T, col.names=NA, quote=F)
	
	cat(paste(Sys.time()), "ITERATION END:",k, "\n", sep="\t")
}

##### COMPILE RESULTS FORM RANDOMIZED GENESETS -------------	
files.rmd <- list.files(dir.rmd, full.names=T)

.list <- list()
for(i in 1:length(files.rmd)){
	mat <- read.delim(files.rmd[i], header=T, row.names=1, stringsAsFactors=F)
	mat[is.na(mat)] <- 0
	.list[[i]] <- mat
}
mat.sum <- Reduce("+", .list)
mat.avg <- (Reduce("+", .list))/length(.list)

write.table(mat.sum, file.path(dir.batch, paste("04_dat_pathway_interaction_matrix_random_sum_", batch, ".tsv", sep="")), sep="\t", row.names=T, col.names=NA, quote=F)
write.table(mat.avg, file.path(dir.batch, paste("05_dat_pathway_interaction_matrix_random_average_", batch, ".tsv", sep="")), sep="\t", row.names=T, col.names=NA, quote=F)

### 5: ESTIMATE CROSS-TALKS BETWEEN PATHWAYS ------------------------------------------------------------------
# Perform one-sided Fisher Exact test on all pathway pairs using
# 	the 2 X 2 contingency table that include the following numbers:
# 	n, N-n, r, R-r where 
# 		n = the interaction count between original pathways, 
# 		N = the number of total interaction counts of all pathway pairs, 
# 		r = the average of interaction counts between the pair of corresponding randomized pathways after 1000 rounds of randomizations and 
# 		R = the average of total interaction counts of all randomized pathway pairs after 1000 rounds of randomizations. 
# 	The null hypothesis is that the ratio of true interactions between two pathways to all interactions
# 	(n/N) is the same as the ratio of random interactions to all
# 	random interactions (r/R). In our analysis, we only focused on
# 	pathway pairs where n/N is significantly higher than r/R. Fisher
# 	exact test P-values were adjusted using FDR BH procedure
# 	(Benjamini and Hochberg, 1995) to account for multiple
# 	hypothesis testing

### True Interactions
file.dat.int <- file.path(dir.batch, paste("02_dat_pathway_interaction_matrix_", batch, ".tsv", sep=""))
dat.int <- read.delim(file.dat.int, header=T, row.names=1, stringsAsFactors=F)
N <- sum(rowSums(dat.int))

df.int <- melt(as.matrix(dat.int))
colnames(df.int) <- c("PathwayA","PathwayB", "Interactions")
df.int$PathwayA <- as.character(df.int$PathwayA)
df.int$PathwayB <- as.character(df.int$PathwayB)

### Random Interactions
file.mat.sum <- file.path(dir.batch, paste("04_dat_pathway_interaction_matrix_random_sum_", batch, ".tsv", sep=""))
mat.sum <- read.delim(file.mat.sum, header=T, row.names=1, stringsAsFactors=F)
R <- sum(rowSums(mat.sum))/100

file.mat.avg <- file.path(dir.batch, paste("05_dat_pathway_interaction_matrix_random_average_", batch, ".tsv", sep=""))
mat.avg <- read.delim(file.mat.avg, header=T, row.names=1, stringsAsFactors=F)

df.rmd <- melt(as.matrix(mat.avg))
colnames(df.rmd) <- c("PathwayA","PathwayB", "Interactions")
df.rmd$PathwayA <- as.character(df.rmd$PathwayA)
df.rmd$PathwayB <- as.character(df.rmd$PathwayB)

# TP: n
# TN: N-n
# FP: r
# FN: R-r

### Combine Data
df <- df.int
colnames(df)[3] <- "TP"
df$TN <- N - df$TP
df$FP <- df.rmd$Interactions
df$FN <- R - df$FP

df <- subset(df, df$TP != 0)
df$n_N <- df$TP/N
df$r_R <- df$FP/R

df$ratio.status=ifelse(df$n_N > df$r_R, 1, 0)

### Fisher's Exact Test
df$pvalue <- apply(df, 1, function(x) fisher.test(matrix(c(as.numeric(x[3]), as.numeric(x[4]), as.numeric(x[5]), as.numeric(x[6])), 2, 2, byrow=T), alternative="greater")$p.value)
df$fdr <- p.adjust(df$pvalue, method="BH")

### 6: SIGNIFICANT PATHWAY CROSSTAKS -----------------------------------------------------------------
# All pathway pairs with adjusted Fisher exact test P-value50.05
# 	were pulled together to construct a network in which a node is a
# 	pathway and an edge represents crosstalk between two pathways.

df$sig <- ifelse(df$fdr < 0.05, 1, 0)     

### NETWOEK PURNING ----------------------------------------------------------------------------------
#Two types of redundant edges caused by gene redundancy
#	among some pathways were removed during the following
#	pruning steps to clean up the network.
#(a) All edges between two pathways with significant gene overlap
#	identified in Step 2 were considered as not informative and thus
#	removed from the network. Note that it is our intent to discover
#	crosstalk among different biological activities. Pathway pairs
#	where both pathways significantly overlap with each other in
#	terms of gene members represent similar biology and thus were
#	excluded from the network.
#(b) Two overlapping pathways may both interact with the same
#	pathway. In this case, the two edges were considered redundant
#	and one of them was removed. For example, pathway A has four
#	neighbors: pathway B, C, D and E. If at least 75% of genes in B
#	are also in C and B has less number of genes than C, then the
#	edge between A and B was removed. Repeat this process until no
#	more edges between pathway A and its neighbors could be
#	removed. This step is especially necessary for pathways derived
#	from GO where all children gene sets are complete subsets of
#	their parents. If both a child and a parent GO pathway interact
#	with pathway A, the edge between pathway A and the child GO
#	pathway will be removed during this step of pruning.


### (a) Load Pathway Overlap data (from step-2)
file.overlap <- file.path(dir.output, "dat_pathway_overlap_kegg.tsv")
dat.overlap <- read.delim(file.overlap, header=T, stringsAsFactors=F)

df$items <- apply(df, 1, function(x) paste(as.character(x[1]),as.character(x[2]), sep=":"))
dat.overlap$items <- apply(dat.overlap, 1, function(x) paste(as.character(x[1]),as.character(x[2]), sep=":"))
dat.overlap <- subset(dat.overlap, dat.overlap$items %in% df$items)
dat.overlap <- dat.overlap[match(df$items, dat.overlap$items),]
df <- df[,-which(colnames(df) == "items")]

df$overlap.sig <- dat.overlap$overlap.sig
df$crosstalk.status <- ifelse((df$ratio.status == 1) & (df$sig == 1) & (df$overlap.sig == 0), 1, 0)
write.table(df, file.path(dir.batch, paste("06_dat_pathway_crosstalk_ftest_summary_", batch, ".tsv", sep="")), sep="\t", row.names=F, col.names=T, quote=F)

### (b) Remove redundant edges
df.ct <- subset(df, df$crosstalk.status == 1)

### create igraph object for the resulting pathways
gpath <- graph.data.frame(df.ct, directed=FALSE)

for(i in 1:length(V(gpath))){
	neigh.path <- names(neighbors(gpath, V(gpath)[i]))
	if(length(neigh.path) == 1) next
	d <- data.frame(t(combn(x=neigh.path, m=2, simplify=T)))
	colnames(d) <- c("PathA","PathB")
	d$PathA <- as.character(d$PathA)
	d$PathB <- as.character(d$PathB)
	
	genes.pathA <- str_split(dat.genesets$Genesets[apply(d, 1, function(x) which(dat.genesets$Category == x[1]))], ",")
	genes.pathB <- str_split(dat.genesets$Genesets[apply(d, 1, function(x) which(dat.genesets$Category == x[2]))], ",")
	
	d$nA <- unlist(lapply(genes.pathA, length))
	d$nB <- unlist(lapply(genes.pathB, length))
	
	if(nrow(d) == 1){
		d$AnB <- length(intersect(genes.pathA[[1]], genes.pathB[[1]]))
	} else{
		d$AnB <- unlist(lapply(mapply(intersect, genes.pathA, genes.pathB), length))
	}
	
	d$redundant <- ifelse(apply(d, 1, function(x) as.numeric(x[5])/min(as.numeric(x[3]), as.numeric(x[4]))) >= 0.75, 1, 0)
	
	d <- subset(d, d$redundant == 1)
	
	if(nrow(d) == 0) next
	
	d$index <- ifelse(d$nA == apply(d, 1, function(x) min(as.numeric(x[3]), as.numeric(x[4]))), 1, 2)
	d$del.path <- apply(d, 1, function(x) x[as.numeric(x[7])])
	gpath <- gpath %>% delete_edges(paste(V(gpath)$name[i], unique(d$del.path), sep="|"))
	#cat("DONE:", i, "\n", sep="\t")
}

# remove isolated nodes
gpath <- delete_vertices(gpath, which(degree(gpath) == 0))

### Pathway Cross-talk Network 
dat.deg <- data.frame(Pathway=names(degree(gpath, v = V(gpath), mode = "all")), Degree=as.numeric(degree(gpath, v = V(gpath), mode = "all")))
dat.deg$Pathway <- as.character(dat.deg$Pathway)  
dat.deg <- dat.deg[order(dat.deg$Degree, decreasing=T),]



##### FUNCTION: COUNT PROTEIN INTERACTION FOR RANDOMIZED GENESETS -------------
compute_protein_interactions <- function(dat, adj){	
	require("stringr", "reshape2","foreach","doParallel")

	### Declate Cluster 
	no_cores <- detectCores() - 10
	cl <- makeCluster(no_cores)
	registerDoParallel(cl)

	lpar <- foreach(i = 1:nrow(dat), .combine='cbind') %:%
			foreach(j = 1:nrow(dat), .combine='c', .packages="stringr") %dopar% {
				if(i > j){
					genes.A <- str_split(dat$Genesets[i], ",")[[1]]
					genes.B <- str_split(dat$Genesets[j], ",")[[1]]
					genes.AB <- intersect(genes.A, genes.B)
					
					genes.A <- setdiff(genes.A, genes.AB)
					genes.B <- setdiff(genes.B, genes.AB)
					
					genes.A <- intersect(genes.A, rownames(adj))
					genes.B <- intersect(genes.B, rownames(adj))
		
					if((length(genes.A) == 0) | (length(genes.B) == 0)){
						lpar <- NA
					} else{
						dat.combn <- expand.grid(GeneA=genes.A, GeneB=genes.B)
						dat.combn$GeneA <- as.character(dat.combn$GeneA)
						dat.combn$GeneB <- as.character(dat.combn$GeneB)
						dat.combn$Edge <- 0
						for(k in 1:nrow(dat.combn)){
							dat.combn$Edge[k] <- as.numeric(adj[dat.combn$GeneA[k], dat.combn$GeneB[k]])
						}
						lpar <- length(which(dat.combn$Edge == 1))
					}					
				} else{
					lpar <- NA
				}
			}
	stopCluster(cl)
	colnames(lpar) <- rownames(lpar) <- dat$Category

	return(lpar)
}

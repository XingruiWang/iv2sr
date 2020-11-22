load("mice/mice_processed_liver.rda")
load("mice/mice_result_lasso.rda")
load("mice/mice_result_scad.rda")
load("mice/mice_result_mcp.rda")
load("mice/mice_stab.rda")
gene.map <- read.delim("mice/gene_map.txt")
mark.map <- read.delim("mice/vanNas_etal_crossI_genotypeData.txt")[, 1:3]
n <- length(y); p <- ncol(x); q <- ncol(z)
cM <- 10; thr <- 0.5

(sum(bet.pls.lasso != 0))
(sum(bet.pls.scad != 0))
(sum(bet.pls.mcp != 0))
(sum(bet.2sr.lasso != 0))
(sum(bet.2sr.scad != 0))
(sum(bet.2sr.mcp != 0))

gene.chr <- gene.map$Chr[match(colnames(x), gene.map$ID)]
gene.pos <- gene.map$cM.Position[match(colnames(x), gene.map$ID)]
mark.chr <- mark.map$chromosome[match(colnames(z), mark.map$Marker_name)]
mark.pos <- mark.map$marker_pos_cm[match(colnames(z), mark.map$Marker_name)]
is.cis <- matrix(, q, p)
is.cis <- matrix(mark.chr[row(is.cis)] == gene.chr[col(is.cis)] &
				 abs(mark.pos[row(is.cis)] - gene.pos[col(is.cis)]) <= cM, q, p)
ind.snp <- (gam.2sr.lasso != 0 | gam.2sr.scad != 0 | gam.2sr.mcp != 0) & is.cis
cis.snp <- character(p)
for (j in 1:p) cis.snp[j] <- paste(colnames(z)[ind.snp[, j]], collapse=",")

ind <- (1:p)[pmax(prob.2sr.lasso, prob.2sr.scad, prob.2sr.mcp) >= thr]
gene.table <- gene.map[match(colnames(x)[ind], gene.map$ID), c(2, 4:6)]
rownames(gene.table) <- ind
gene.table <- cbind(gene.table, Prob.PLS.Lasso=prob.pls.lasso[ind],
					Prob.PLS.SCAD=prob.pls.scad[ind], Prob.PLS.MCP=prob.pls.mcp[ind],
					Prob.2SR.Lasso=prob.2sr.lasso[ind], Prob.2SR.SCAD=prob.2sr.scad[ind],
					Prob.2SR.MCP=prob.2sr.mcp[ind])
gene.table <- cbind(gene.table, cis.SNPs=cis.snp[ind])
gene.table <- gene.table[order(gene.table$cM.Position), ]
gene.table <- gene.table[order(match(gene.table$Chr, c(1:19, "X"))), ]
write.table(gene.table, "mice/gene_table.txt", quote=FALSE, sep="\t", row=FALSE)

r2.pls.lasso <- cor(y, x %*% bet.pls.lasso)
(r2adj.pls.lasso <- 1 - (1 - r2.pls.lasso)*(n - 1)/(n - sum(bet.pls.lasso != 0) - 1))
r2.pls.scad <- cor(y, x %*% bet.pls.scad)
(r2adj.pls.scad <- 1 - (1 - r2.pls.scad)*(n - 1)/(n - sum(bet.pls.scad != 0) - 1))
r2.pls.mcp <- cor(y, x %*% bet.pls.mcp)
(r2adj.pls.mcp <- 1 - (1 - r2.pls.mcp)*(n - 1)/(n - sum(bet.pls.mcp != 0) - 1))

r2.2sr.lasso <- cor(y, x %*% bet.2sr.lasso)
(r2adj.2sr.lasso <- 1 - (1 - r2.2sr.lasso)*(n - 1)/(n - sum(bet.2sr.lasso != 0) - 1))
r2.2sr.scad <- cor(y, x %*% bet.2sr.scad)
(r2adj.2sr.scad <- 1 - (1 - r2.2sr.scad)*(n - 1)/(n - sum(bet.2sr.scad != 0) - 1))
r2.2sr.mcp <- cor(y, x %*% bet.2sr.mcp)
(r2adj.2sr.mcp <- 1 - (1 - r2.2sr.mcp)*(n - 1)/(n - sum(bet.2sr.mcp != 0) - 1))

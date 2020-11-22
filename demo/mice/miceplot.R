load("mice/mice_processed_liver.rda")
load("mice/mice_result_lasso.rda")
load("mice/mice_result_scad.rda")
load("mice/mice_result_mcp.rda")
load("mice/mice_stab.rda")
mark.map <- read.delim("mice/vanNas_etal_crossI_genotypeData.txt")[, 1:3]

p <- length(prob.pls.lasso); thr <- 0.4
set1 <- (1:p)[prob.pls.lasso >= thr]
set2 <- (1:p)[prob.2sr.lasso >= thr]
set.pls.lasso <- setdiff(set1, set2)
set.2sr.lasso <- setdiff(set2, set1)
set.int.lasso <- intersect(set1, set2)
set.nei.lasso <- setdiff(1:p, union(set1, set2))

set1 <- (1:p)[prob.pls.scad >= thr]
set2 <- (1:p)[prob.2sr.scad >= thr]
set.pls.scad <- setdiff(set1, set2)
set.2sr.scad <- setdiff(set2, set1)
set.int.scad <- intersect(set1, set2)
set.nei.scad <- setdiff(1:p, union(set1, set2))

set1 <- (1:p)[prob.pls.mcp >= thr]
set2 <- (1:p)[prob.2sr.mcp >= thr]
set.pls.mcp <- setdiff(set1, set2)
set.2sr.mcp <- setdiff(set2, set1)
set.int.mcp <- intersect(set1, set2)
set.nei.mcp <- setdiff(1:p, union(set1, set2))

pdf.options(width=5.8, height=4.7)
setEPS(width=5.8, height=4.7)
pdf.out <- TRUE

if (pdf.out) pdf(file="fig/path_pls_lasso.pdf") else postscript(file="fig/path_pls_lasso.eps")
par(mai=c(0.7, 0.8, 0.4, 0.2), mgp=c(2.1, 0.7, 0))
lam.seq <- seq(0, by=log(50)/49, along=lam.pls.lasso)
plot(lam.seq, path.pls.lasso[1, ], type="n", ylim=c(0, 1), main="PLS with Lasso", cex.main=1.2,
	 xlab=expression(-log(mu/mu[max])), ylab="Selection probability", cex.lab=1.2)
for (i in set.nei.lasso) lines(lam.seq, path.pls.lasso[i, ], lty=2)
for (i in set.pls.lasso) lines(lam.seq, path.pls.lasso[i, ], col="blue")
for (i in set.int.lasso) lines(lam.seq, path.pls.lasso[i, ], col="red")
dev.off()

if (pdf.out) pdf(file="fig/path_2sr_lasso.pdf") else postscript(file="fig/path_2sr_lasso.eps")
par(mai=c(0.7, 0.8, 0.4, 0.2), mgp=c(2.1, 0.7, 0))
lam.seq <- seq(0, by=log(50)/49, along=lam.2sr.lasso)
plot(lam.seq, path.2sr.lasso[1, ], type="n", ylim=c(0, 1), main="2SR with Lasso", cex.main=1.2,
	 xlab=expression(-log(mu/mu[max])), ylab="Selection probability", cex.lab=1.2)
for (i in set.nei.lasso) lines(lam.seq, path.2sr.lasso[i, ], lty=2)
for (i in set.2sr.lasso) lines(lam.seq, path.2sr.lasso[i, ], col="blue")
for (i in set.int.lasso) lines(lam.seq, path.2sr.lasso[i, ], col="red")
dev.off()

if (pdf.out) pdf(file="fig/path_pls_scad.pdf") else postscript(file="fig/path_pls_scad.eps")
par(mai=c(0.7, 0.8, 0.4, 0.2), mgp=c(2.1, 0.7, 0))
lam.seq <- seq(0, by=log(50)/49, along=lam.pls.scad)
plot(lam.seq, path.pls.scad[1, ], type="n", ylim=c(0, 1), main="PLS with SCAD", cex.main=1.2,
	 xlab=expression(-log(mu/mu[max])), ylab="Selection probability", cex.lab=1.2)
for (i in set.nei.scad) lines(lam.seq, path.pls.scad[i, ], lty=2)
for (i in set.pls.scad) lines(lam.seq, path.pls.scad[i, ], col="blue")
for (i in set.int.scad) lines(lam.seq, path.pls.scad[i, ], col="red")
dev.off()

if (pdf.out) pdf(file="fig/path_2sr_scad.pdf") else postscript(file="fig/path_2sr_scad.eps")
par(mai=c(0.7, 0.8, 0.4, 0.2), mgp=c(2.1, 0.7, 0))
lam.seq <- seq(0, by=log(50)/49, along=lam.2sr.scad)
plot(lam.seq, path.2sr.scad[1, ], type="n", ylim=c(0, 1), main="2SR with SCAD", cex.main=1.2,
	 xlab=expression(-log(mu/mu[max])), ylab="Selection probability", cex.lab=1.2)
for (i in set.nei.scad) lines(lam.seq, path.2sr.scad[i, ], lty=2)
for (i in set.2sr.scad) lines(lam.seq, path.2sr.scad[i, ], col="blue")
for (i in set.int.scad) lines(lam.seq, path.2sr.scad[i, ], col="red")
dev.off()

if (pdf.out) pdf(file="fig/path_pls_mcp.pdf") else postscript(file="fig/path_pls_mcp.eps")
par(mai=c(0.7, 0.8, 0.4, 0.2), mgp=c(2.1, 0.7, 0))
lam.seq <- seq(0, by=log(50)/49, along=lam.pls.mcp)
plot(lam.seq, path.pls.mcp[1, ], type="n", ylim=c(0, 1), main="PLS with MCP", cex.main=1.2,
	 xlab=expression(-log(mu/mu[max])), ylab="Selection probability", cex.lab=1.2)
for (i in set.nei.mcp) lines(lam.seq, path.pls.mcp[i, ], lty=2)
for (i in set.pls.mcp) lines(lam.seq, path.pls.mcp[i, ], col="blue")
for (i in set.int.mcp) lines(lam.seq, path.pls.mcp[i, ], col="red")
dev.off()

if (pdf.out) pdf(file="fig/path_2sr_mcp.pdf") else postscript(file="fig/path_2sr_mcp.eps")
par(mai=c(0.7, 0.8, 0.4, 0.2), mgp=c(2.1, 0.7, 0))
lam.seq <- seq(0, by=log(50)/49, along=lam.2sr.mcp)
plot(lam.seq, path.2sr.mcp[1, ], type="n", ylim=c(0, 1), main="2SR with MCP", cex.main=1.2,
	 xlab=expression(-log(mu/mu[max])), ylab="Selection probability", cex.lab=1.2)
for (i in set.nei.mcp) lines(lam.seq, path.2sr.mcp[i, ], lty=2)
for (i in set.2sr.mcp) lines(lam.seq, path.2sr.mcp[i, ], col="blue")
for (i in set.int.mcp) lines(lam.seq, path.2sr.mcp[i, ], col="red")
dev.off()

chr.size <- c(98.5, 103.9, 82.7, 88.6, 90.2, 79, 89.1, 76.2, 75.1, 77.9, 88.0, 63.9, 67.3, 66.4,
			  59, 57.8, 61.3, 59.4, 56.9)
chr.start <- c(0, cumsum(chr.size))
mark.chr <- mark.map$chromosome[match(colnames(z), mark.map$Marker_name)]
mark.pos <- mark.map$marker_pos_cm[match(colnames(z), mark.map$Marker_name)]
geno.pos <- (chr.start[mark.chr] + mark.pos)

eff.bgr.lasso <- rowSums(abs(gam.2sr.lasso))
eff.bgr.scad <- rowSums(abs(gam.2sr.scad))
eff.bgr.mcp <- rowSums(abs(gam.2sr.mcp))
ind <- bet.2sr.lasso != 0 | bet.2sr.scad != 0 | bet.2sr.mcp != 0
eff.sel.lasso <- rowSums(abs(gam.2sr.lasso[, ind]))
eff.sel.scad <- rowSums(abs(gam.2sr.scad[, ind]))
eff.sel.mcp <- rowSums(abs(gam.2sr.mcp[, ind]))

snp.bgr.ind <- order(pmax(eff.bgr.lasso, eff.bgr.scad, eff.bgr.mcp), decreasing=TRUE)[1:5]
snp.bgr.table <- mark.map[match(colnames(z)[snp.bgr.ind], mark.map$Marker_name), 1:3]
snp.bgr.table <- cbind(snp.bgr.table, geno.pos=geno.pos[snp.bgr.ind])
snp.sel.ind <- order(pmax(eff.sel.lasso, eff.sel.scad, eff.sel.mcp), decreasing=TRUE)[1:5]
snp.sel.table <- mark.map[match(colnames(z)[snp.sel.ind], mark.map$Marker_name), 1:3]
snp.sel.table <- cbind(snp.sel.table, geno.pos=geno.pos[snp.sel.ind])

pdf.options(width=11.8, height=4.7)
setEPS(width=11.8, height=4.7)
pdf.out <- TRUE

if (pdf.out) pdf(file="fig/snp_bgr.pdf") else postscript(file="fig/snp_bgr.eps")
par(mai=c(0.7, 0.8, 0.4, 0.2), mgp=c(2.1, 0.7, 0))
plot(geno.pos, eff.bgr.lasso, type="n", cex.lab=1.2, xlab="Genome position (cM)",
	 ylab=expression(paste(italic(L)[1], " effect size")), lab=c(15, 5, 7))
lines(geno.pos, eff.bgr.lasso)
lines(geno.pos, eff.bgr.scad, lty=2)
lines(geno.pos, eff.bgr.mcp, lty=3)
legend("topleft", c("Lasso", "SCAD", "MCP"), lty=1:4, inset=c(0.02, 0.06))
dev.off()

if (pdf.out) pdf(file="fig/snp_sel.pdf") else postscript(file="fig/snp_sel.eps")
par(mai=c(0.7, 0.8, 0.4, 0.2), mgp=c(2.1, 0.7, 0))
plot(geno.pos, eff.sel.scad, type="n", cex.lab=1.2, xlab="Genome position (cM)",
	 ylab=expression(paste(italic(L)[1], " effect size")), lab=c(15, 4, 7))
lines(geno.pos, eff.sel.lasso)
lines(geno.pos, eff.sel.scad, lty=2)
lines(geno.pos, eff.sel.mcp, lty=3)
legend("topleft", c("Lasso", "SCAD", "MCP"), lty=1:4, inset=c(0.02, 0.06))
dev.off()

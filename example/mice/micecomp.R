set.seed(43)
load("mice/mice_processed_liver.rda")
n <- length(y); p <- ncol(x); q <- ncol(z); nfold <- 10
ind <- split(1:n, sample(rep(1:nfold, length=n)))
#lasso <- TRUE; scad <- FALSE; mcp <- FALSE

x <- scale(x, scale=FALSE)
y <- scale(y, scale=FALSE)
z <- scale(z, scale=FALSE)

if (lasso) {
	ans <- cd(y, x, "Lasso")
	sol.pls.lasso <- ans$sol; lam <- ans$lam
	i.lam <- cv(y, x, "Lasso", lam, ind)$i
	bet.pls.lasso <- sol.pls.lasso[, i.lam]
	
	cat("2SR-Lasso\n")
	gam.2sr.lasso <- matrix(, q, p)
	for (j in 1:p) {
		if (j %% 100 == 0) cat(round(j/p*100), "% complete.\n", sep="")
		ans <- cd(x[, j], z, "Lasso", upto=0.25*n)
		sol <- ans$sol; lam <- ans$lam
		i.lam <- cv(x[, j], z, "Lasso", lam, ind)$i
		gam.2sr.lasso[, j] <- sol[, i.lam]
	}
	x.hat <- z %*% gam.2sr.lasso
	ind.x <- colSums(gam.2sr.lasso != 0) != 0
	ans <- cd(y, x.hat[, ind.x], "Lasso", upto=0.25*n)
	lam <- ans$lam; sol.2sr.lasso <- matrix(0, p, length(lam))
	sol.2sr.lasso[ind.x, ] <- ans$sol
	i.lam <- cv(y, x.hat[, ind.x], "Lasso", lam, ind)$i
	bet.2sr.lasso <- sol.2sr.lasso[, i.lam]
	
	save(bet.pls.lasso, bet.2sr.lasso, gam.2sr.lasso, sol.pls.lasso, sol.2sr.lasso,
		 file="mice/mice_result_lasso.rda")
}

if (scad) {
	ans <- cd(y, x, "SCAD")
	sol.pls.scad <- ans$sol; lam <- ans$lam
	i.lam <- cv(y, x, "SCAD", lam, ind)$i
	bet.pls.scad <- sol.pls.scad[, i.lam]

	cat("2SR-SCAD\n")
	gam.2sr.scad <- matrix(, q, p)
	for (j in 1:p) {
		if (j %% 100 == 0) cat(round(j/p*100), "% complete.\n", sep="")
		ans <- cd(x[, j], z, "SCAD", upto=0.2*n)
		sol <- ans$sol; lam <- ans$lam
		i.lam <- cv(x[, j], z, "SCAD", lam, ind)$i
		gam.2sr.scad[, j] <- sol[, i.lam]
	}
	x.hat <- z %*% gam.2sr.scad
	ind.x <- colSums(gam.2sr.scad != 0) != 0
	ans <- cd(y, x.hat[, ind.x], "SCAD", upto=0.2*n)
	lam <- ans$lam; sol.2sr.scad <- matrix(0, p, length(lam))
	sol.2sr.scad[ind.x, ] <- ans$sol
	i.lam <- cv(y, x.hat[, ind.x], "SCAD", lam, ind)$i
	bet.2sr.scad <- sol.2sr.scad[, i.lam]
	
	save(bet.pls.scad, bet.2sr.scad, gam.2sr.scad, sol.pls.scad, sol.2sr.scad,
		 file="mice/mice_result_scad.rda")
}

if (mcp) {
	ans <- cd(y, x, "MCP")
	sol.pls.mcp <- ans$sol; lam <- ans$lam
	i.lam <- cv(y, x, "MCP", lam, ind)$i
	bet.pls.mcp <- sol.pls.mcp[, i.lam]

	cat("2SR-MCP\n")
	gam.2sr.mcp <- matrix(, q, p)
	for (j in 1:p) {
		if (j %% 100 == 0) cat(round(j/p*100), "% complete.\n", sep="")
		ans <- cd(x[, j], z, "MCP", upto=0.15*n)
		sol <- ans$sol; lam <- ans$lam
		i.lam <- cv(x[, j], z, "MCP", lam, ind)$i
		gam.2sr.mcp[, j] <- sol[, i.lam]
	}
	x.hat <- z %*% gam.2sr.mcp
	ind.x <- colSums(gam.2sr.mcp != 0) != 0
	ans <- cd(y, x.hat[, ind.x], "MCP", upto=0.15*n)
	lam <- ans$lam; sol.2sr.mcp <- matrix(0, p, length(lam))
	sol.2sr.mcp[ind.x, ] <- ans$sol
	i.lam <- cv(y, x.hat[, ind.x], "MCP", lam, ind)$i
	bet.2sr.mcp <- sol.2sr.mcp[, i.lam]
	
	save(bet.pls.mcp, bet.2sr.mcp, gam.2sr.mcp, sol.pls.mcp, sol.2sr.mcp,
		 file="mice/mice_result_mcp.rda")
}

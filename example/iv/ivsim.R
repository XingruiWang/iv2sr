gen.data <- function(n, p, q, s=5, s0=10, r=5, r0=50, weak=FALSE, mixed=FALSE) {
	require(MASS)
	gam <- matrix(0, q, p)
	for (j in 1:p)
		if (mixed) {
			ind <- sample(1:q, r)
			gam[ind, j] <- (2*rbinom(r, 1, 0.5) - 1)*runif(r, 0.5, 1)
			gam[sample(setdiff(1:q, ind), r0 - r), j] <- (2*rbinom(r0 - r, 1, 0.5) - 1)*runif(r0 - r, 0.05, 0.1)
		} else if (weak)
			gam[sample(1:q, r), j] <- (2*rbinom(r, 1, 0.5) - 1)*runif(r, 0.5, 0.75)
		else
			gam[sample(1:q, r), j] <- (2*rbinom(r, 1, 0.5) - 1)*runif(r, 0.75, 1)
	bet <- rep(0, p)
	ind <- sample(1:p, s)
	bet[ind] <- (2*rbinom(s, 1, 0.5) - 1)*runif(s, 0.5, 1)
	
	sig <- matrix(, p + 1, p + 1)
	sig <- 0.2^abs(row(sig) - col(sig))
	sig[1:p, p + 1] <- 0
	sig[c(ind, sample(setdiff(1:p, ind), s0 - s)), p + 1] <- 0.3
	sig[p + 1, 1:p] <- sig[1:p, p + 1]
	
	e <- mvrnorm(n, rep(0, p + 1), sig)
	if (mixed)
		z <- matrix(rbinom(n*q, 1, runif(n*q, 0, 0.5)), n, q)
	else
		z <- matrix(rbinom(n*q, 1, 0.5), n, q)
	x <- z %*% gam + e[, 1:p]
	y <- drop(x %*% bet) + e[, p + 1]
	
	list(x=x, y=y, z=z, bet=bet, gam=gam)
}

ivgen <- function(n, p, q, nsim=50, seed=19, weak=FALSE, mixed=FALSE) {
	set.seed(seed)
	data <- vector("list", nsim)
	for (i in 1:nsim) {
		data[[i]] <- gen.data(n, p, q, weak=weak, mixed=mixed)
		cat("Sim.", i, "generated.\n")
	}
	fn <- paste(ifelse(weak, "w", ifelse(mixed, "m", "s")), paste(n, p, q, sep="_"), ".rda", sep="")
	save(data, file=paste("data/", fn, sep=""))
}

ivsim <- function(n, p, q, nsim=50, nfold=10, seed=42, weak=FALSE, mixed=FALSE, batch) {
	set.seed(seed)
	fn <- paste(ifelse(weak, "w", ifelse(mixed, "m", "s")), paste(n, p, q, sep="_"), sep="")
	load(paste("data/", fn, ".rda", sep=""))
	if (missing(batch))	result <- vector("list", nsim)
	
	for (i in 1:nsim) {
		ind <- split(1:n, sample(rep(1:nfold, length=n)))
		if (!missing(batch) && batch != i) next
		
		x <- data[[i]]$x; y <- data[[i]]$y; z <- data[[i]]$z
		bet <- data[[i]]$bet; gam <- data[[i]]$gam
		x <- scale(x, scale=FALSE)
		y <- scale(y, scale=FALSE)
		z <- scale(z, scale=FALSE)
		
		# PLS-oracle
		supp <- bet != 0
		y1 <- scale(y, scale=FALSE); x1 <- x[, supp]
		bet.pls.oracle <- rep(0, p)
		bet.pls.oracle[supp] <- solve(crossprod(x1), crossprod(x1, y1))
		
		# PLS-Lasso
		ans <- cd(y, x, "Lasso")
		sol.pls.lasso <- ans$sol; lam <- ans$lam
		i.lam <- cv(y, x, "Lasso", lam, ind)$i
		bet.pls.lasso <- sol.pls.lasso[, i.lam]
		
		# PLS-SCAD
		ans <- cd(y, x, "SCAD")
		sol.pls.scad <- ans$sol; lam <- ans$lam
		i.lam <- cv(y, x, "SCAD", lam, ind)$i
		bet.pls.scad <- sol.pls.scad[, i.lam]
		
		# PLS-MCP
		ans <- cd(y, x, "MCP")
		sol.pls.mcp <- ans$sol; lam <- ans$lam
		i.lam <- cv(y, x, "MCP", lam, ind)$i
		bet.pls.mcp <- sol.pls.mcp[, i.lam]
		
		# 2SR-oracle
		x1 <- scale(x, scale=FALSE)
		gam.2sr.oracle <- matrix(0, q, p)
		for (j in 1:p) {
			supp <- gam[, j] != 0
			z1 <- z[, supp]
			gam.2sr.oracle[supp, j] <- solve(crossprod(z1), crossprod(z1, x1[, j]))
		}
		x.hat <- z %*% gam.2sr.oracle
		supp <- bet != 0
		x1.hat <- x.hat[, supp]
		y1 <- scale(y, scale=FALSE)
		bet.2sr.oracle <- rep(0, p)
		bet.2sr.oracle[supp] <- solve(crossprod(x1.hat), crossprod(x1.hat, y1))
		
		# 2SR-Lasso
		cat("2SR-Lasso\n")
		gam.2sr.lasso <- matrix(, q, p)
		for (j in 1:p) {
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
		
		# 2SR-SCAD
		cat("2SR-SCAD\n")
		gam.2sr.scad <- matrix(, q, p)
		for (j in 1:p) {
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
		
		# 2SR-MCP
		cat("2SR-MCP\n")
		gam.2sr.mcp <- matrix(, q, p)
		for (j in 1:p) {
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
		
		if (missing(batch))
			result[[i]] <- list(bet.pls.oracle=bet.pls.oracle, bet.pls.lasso=bet.pls.lasso,
				bet.pls.scad=bet.pls.scad, bet.pls.mcp=bet.pls.mcp,
				bet.2sr.oracle=bet.2sr.oracle, bet.2sr.lasso=bet.2sr.lasso,
 				bet.2sr.scad=bet.2sr.scad, bet.2sr.mcp=bet.2sr.mcp,
				gam.2sr.oracle=gam.2sr.oracle, gam.2sr.lasso=gam.2sr.lasso,
 				gam.2sr.scad=gam.2sr.scad, gam.2sr.mcp=gam.2sr.mcp,
				sol.pls.lasso=sol.pls.lasso, sol.pls.scad=sol.pls.scad,
				sol.pls.mcp=sol.pls.mcp, sol.2sr.lasso=sol.2sr.lasso,
				sol.2sr.scad=sol.2sr.scad, sol.2sr.mcp=sol.2sr.mcp)
		else
			result1 <- list(bet.pls.oracle=bet.pls.oracle, bet.pls.lasso=bet.pls.lasso,
				bet.pls.scad=bet.pls.scad, bet.pls.mcp=bet.pls.mcp,
				bet.2sr.oracle=bet.2sr.oracle, bet.2sr.lasso=bet.2sr.lasso,
 				bet.2sr.scad=bet.2sr.scad, bet.2sr.mcp=bet.2sr.mcp,
				gam.2sr.oracle=gam.2sr.oracle, gam.2sr.lasso=gam.2sr.lasso,
 				gam.2sr.scad=gam.2sr.scad, gam.2sr.mcp=gam.2sr.mcp,
				sol.pls.lasso=sol.pls.lasso, sol.pls.scad=sol.pls.scad,
				sol.pls.mcp=sol.pls.mcp, sol.2sr.lasso=sol.2sr.lasso,
				sol.2sr.scad=sol.2sr.scad, sol.2sr.mcp=sol.2sr.mcp)
		cat("Sim.", i, "done.\n")
	}
	if (missing(batch))
		save(result, file=paste("result/", fn, ".rda", sep=""))
	else
		save(result1, file=paste("result/", fn, "b", batch, ".rda", sep=""))
}

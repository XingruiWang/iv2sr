load("mice/mice_processed_liver.rda")
load("mice/mice_result_lasso.rda")
load("mice/mice_result_scad.rda")
load("mice/mice_result_mcp.rda")

x <- scale(x, scale=FALSE)
y <- scale(y, scale=FALSE)
z <- scale(z, scale=FALSE)
p <- ncol(x)

lam.pls.lasso <- cd(y, x, "Lasso", upto=Inf)$lam
ans <- stab(y, x, "Lasso", lam.pls.lasso)
path.pls.lasso <- ans$path; prob.pls.lasso <- ans$prob

x.hat <- z %*% gam.2sr.lasso
ind.x <- colSums(gam.2sr.lasso != 0) != 0
lam.2sr.lasso <- cd(y, x.hat[, ind.x], "Lasso", upto=Inf)$lam
ans <- stab(y, x.hat[, ind.x], "Lasso", lam.2sr.lasso)
path.2sr.lasso <- matrix(0, p, length(lam.2sr.lasso))
prob.2sr.lasso <- numeric(p)
path.2sr.lasso[ind.x, ] <- ans$path
prob.2sr.lasso[ind.x] <- ans$prob

lam.pls.scad <- cd(y, x, "SCAD", upto=Inf)$lam
ans <- stab(y, x, "SCAD", lam.pls.scad)
path.pls.scad <- ans$path; prob.pls.scad <- ans$prob

x.hat <- z %*% gam.2sr.scad
ind.x <- colSums(gam.2sr.scad != 0) != 0
lam.2sr.scad <- cd(y, x.hat[, ind.x], "SCAD", upto=Inf)$lam
ans <- stab(y, x.hat[, ind.x], "SCAD", lam.2sr.scad)
path.2sr.scad <- matrix(0, p, length(lam.2sr.scad))
prob.2sr.scad <- numeric(p)
path.2sr.scad[ind.x, ] <- ans$path
prob.2sr.scad[ind.x] <- ans$prob

lam.pls.mcp <- cd(y, x, "MCP", upto=Inf)$lam
ans <- stab(y, x, "MCP", lam.pls.mcp)
path.pls.mcp <- ans$path; prob.pls.mcp <- ans$prob

x.hat <- z %*% gam.2sr.mcp
ind.x <- colSums(gam.2sr.mcp != 0) != 0
lam.2sr.mcp <- cd(y, x.hat[, ind.x], "MCP", upto=Inf)$lam
ans <- stab(y, x.hat[, ind.x], "MCP", lam.2sr.mcp)
path.2sr.mcp <- matrix(0, p, length(lam.2sr.mcp))
prob.2sr.mcp <- numeric(p)
path.2sr.mcp[ind.x, ] <- ans$path
prob.2sr.mcp[ind.x] <- ans$prob

save(path.pls.lasso, path.pls.scad, path.pls.mcp, path.2sr.lasso, path.2sr.scad, path.2sr.mcp,
	 prob.pls.lasso, prob.pls.scad, prob.pls.mcp, prob.2sr.lasso, prob.2sr.scad, prob.2sr.mcp,
	 lam.pls.lasso, lam.pls.scad, lam.pls.mcp, lam.2sr.lasso, lam.2sr.scad, lam.2sr.mcp,
	 file="mice/mice_stab.rda")

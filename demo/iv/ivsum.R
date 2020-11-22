#fn <- "s200_100_100"
load(paste("data/", fn, ".rda", sep = ""))
load(paste("result/", fn, ".rda", sep = ""))
nsim <- length(data)
n <- length(data[[1]]$y)

est1.pls.oracle <- numeric(nsim)
est2.pls.oracle <- numeric(nsim)
pred.pls.oracle <- numeric(nsim)

est1.pls.lasso <- numeric(nsim)
est2.pls.lasso <- numeric(nsim)
pred.pls.lasso <- numeric(nsim)
tp.pls.lasso <- numeric(nsim)
size.pls.lasso <- numeric(nsim)
mcc.pls.lasso <- numeric(nsim)

est1.pls.scad <- numeric(nsim)
est2.pls.scad <- numeric(nsim)
pred.pls.scad <- numeric(nsim)
tp.pls.scad <- numeric(nsim)
size.pls.scad <- numeric(nsim)
mcc.pls.scad <- numeric(nsim)

est1.pls.mcp <- numeric(nsim)
est2.pls.mcp <- numeric(nsim)
pred.pls.mcp <- numeric(nsim)
tp.pls.mcp <- numeric(nsim)
size.pls.mcp <- numeric(nsim)
mcc.pls.mcp <- numeric(nsim)

est1.2sr.oracle <- numeric(nsim)
est2.2sr.oracle <- numeric(nsim)
pred.2sr.oracle <- numeric(nsim)

est1.2sr.lasso <- numeric(nsim)
est2.2sr.lasso <- numeric(nsim)
pred.2sr.lasso <- numeric(nsim)
tp.2sr.lasso <- numeric(nsim)
size.2sr.lasso <- numeric(nsim)
mcc.2sr.lasso <- numeric(nsim)

est1.2sr.scad <- numeric(nsim)
est2.2sr.scad <- numeric(nsim)
pred.2sr.scad <- numeric(nsim)
tp.2sr.scad <- numeric(nsim)
size.2sr.scad <- numeric(nsim)
mcc.2sr.scad <- numeric(nsim)

est1.2sr.mcp <- numeric(nsim)
est2.2sr.mcp <- numeric(nsim)
pred.2sr.mcp <- numeric(nsim)
tp.2sr.mcp <- numeric(nsim)
size.2sr.mcp <- numeric(nsim)
mcc.2sr.mcp <- numeric(nsim)

for (i in 1:nsim) {
	x <- data[[i]]$x; bet <- data[[i]]$bet
	
	bet.pls.oracle <- result[[i]]$bet.pls.oracle
	bet.pls.lasso <- result[[i]]$bet.pls.lasso
	bet.pls.scad <- result[[i]]$bet.pls.scad
	bet.pls.mcp <- result[[i]]$bet.pls.mcp
	
	bet.2sr.oracle <- result[[i]]$bet.2sr.oracle
	bet.2sr.lasso <- result[[i]]$bet.2sr.lasso
	bet.2sr.scad <- result[[i]]$bet.2sr.scad
	bet.2sr.mcp <- result[[i]]$bet.2sr.mcp
	
	gam.2sr.oracle <- result[[i]]$gam.2sr.oracle
	gam.2sr.lasso <- result[[i]]$gam.2sr.lasso
	gam.2sr.scad <- result[[i]]$gam.2sr.scad
	gam.2sr.mcp <- result[[i]]$gam.2sr.mcp
	
	sol.pls.oracle <- result[[i]]$sol.pls.oracle
	sol.pls.lasso <- result[[i]]$sol.pls.lasso
	sol.pls.scad <- result[[i]]$sol.pls.scad
	sol.pls.mcp <- result[[i]]$sol.pls.mcp
	
	sol.2sr.oracle <- result[[i]]$sol.2sr.oracle
	sol.2sr.lasso <- result[[i]]$sol.2sr.lasso
	sol.2sr.scad <- result[[i]]$sol.2sr.scad
	sol.2sr.mcp <- result[[i]]$sol.2sr.mcp
	
	est1.pls.oracle[i] <- sum(abs(bet.pls.oracle - bet))
	est2.pls.oracle[i] <- sqrt(sum((bet.pls.oracle - bet)^2))
	pred.pls.oracle[i] <- sqrt(sum((x %*% (bet.pls.oracle - bet))^2)/n)
	
	est1.pls.lasso[i] <- sum(abs(bet.pls.lasso - bet))
	est2.pls.lasso[i] <- sqrt(sum((bet.pls.lasso - bet)^2))
	pred.pls.lasso[i] <- sqrt(sum((x %*% (bet.pls.lasso - bet))^2)/n)
	tp.pls.lasso[i] <- sum(bet.pls.lasso != 0 & bet != 0)
	size.pls.lasso[i] <- sum(bet.pls.lasso != 0)
	tp <- sum(bet.pls.lasso != 0 & bet != 0); tn <- sum(bet.pls.lasso == 0 & bet == 0)
	fp <- sum(bet.pls.lasso != 0 & bet == 0); fn <- sum(bet.pls.lasso == 0 & bet != 0)
	mcc.pls.lasso[i] <- (tp*tn - fp*fn)/sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))
	
	est1.pls.scad[i] <- sum(abs(bet.pls.scad - bet))
	est2.pls.scad[i] <- sqrt(sum((bet.pls.scad - bet)^2))
	pred.pls.scad[i] <- sqrt(sum((x %*% (bet.pls.scad - bet))^2)/n)
	tp.pls.scad[i] <- sum(bet.pls.scad != 0 & bet != 0)
	size.pls.scad[i] <- sum(bet.pls.scad != 0)
	tp <- sum(bet.pls.scad != 0 & bet != 0); tn <- sum(bet.pls.scad == 0 & bet == 0)
	fp <- sum(bet.pls.scad != 0 & bet == 0); fn <- sum(bet.pls.scad == 0 & bet != 0)
	mcc.pls.scad[i] <- (tp*tn - fp*fn)/sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))
	
	est1.pls.mcp[i] <- sum(abs(bet.pls.mcp - bet))
	est2.pls.mcp[i] <- sqrt(sum((bet.pls.mcp - bet)^2))
	pred.pls.mcp[i] <- sqrt(sum((x %*% (bet.pls.mcp - bet))^2)/n)
	tp.pls.mcp[i] <- sum(bet.pls.mcp != 0 & bet != 0)
	size.pls.mcp[i] <- sum(bet.pls.mcp != 0)
	tp <- sum(bet.pls.mcp != 0 & bet != 0); tn <- sum(bet.pls.mcp == 0 & bet == 0)
	fp <- sum(bet.pls.mcp != 0 & bet == 0); fn <- sum(bet.pls.mcp == 0 & bet != 0)
	mcc.pls.mcp[i] <- (tp*tn - fp*fn)/sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))
	
	est1.2sr.oracle[i] <- sum(abs(bet.2sr.oracle - bet))
	est2.2sr.oracle[i] <- sqrt(sum((bet.2sr.oracle - bet)^2))
	pred.2sr.oracle[i] <- sqrt(sum((x %*% (bet.2sr.oracle - bet))^2)/n)
	
	est1.2sr.lasso[i] <- sum(abs(bet.2sr.lasso - bet))
	est2.2sr.lasso[i] <- sqrt(sum((bet.2sr.lasso - bet)^2))
	pred.2sr.lasso[i] <- sqrt(sum((x %*% (bet.2sr.lasso - bet))^2)/n)
	tp.2sr.lasso[i] <- sum(bet.2sr.lasso != 0 & bet != 0)
	size.2sr.lasso[i] <- sum(bet.2sr.lasso != 0)
	tp <- sum(bet.2sr.lasso != 0 & bet != 0); tn <- sum(bet.2sr.lasso == 0 & bet == 0)
	fp <- sum(bet.2sr.lasso != 0 & bet == 0); fn <- sum(bet.2sr.lasso == 0 & bet != 0)
	mcc.2sr.lasso[i] <- (tp*tn - fp*fn)/sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))
	
	est1.2sr.scad[i] <- sum(abs(bet.2sr.scad - bet))
	est2.2sr.scad[i] <- sqrt(sum((bet.2sr.scad - bet)^2))
	pred.2sr.scad[i] <- sqrt(sum((x %*% (bet.2sr.scad - bet))^2)/n)
	tp.2sr.scad[i] <- sum(bet.2sr.scad != 0 & bet != 0)
	size.2sr.scad[i] <- sum(bet.2sr.scad != 0)
	tp <- sum(bet.2sr.scad != 0 & bet != 0); tn <- sum(bet.2sr.scad == 0 & bet == 0)
	fp <- sum(bet.2sr.scad != 0 & bet == 0); fn <- sum(bet.2sr.scad == 0 & bet != 0)
	mcc.2sr.scad[i] <- (tp*tn - fp*fn)/sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))
	
	est1.2sr.mcp[i] <- sum(abs(bet.2sr.mcp - bet))
	est2.2sr.mcp[i] <- sqrt(sum((bet.2sr.mcp - bet)^2))
	pred.2sr.mcp[i] <- sqrt(sum((x %*% (bet.2sr.mcp - bet))^2)/n)
	tp.2sr.mcp[i] <- sum(bet.2sr.mcp != 0 & bet != 0)
	size.2sr.mcp[i] <- sum(bet.2sr.mcp != 0)
	tp <- sum(bet.2sr.mcp != 0 & bet != 0); tn <- sum(bet.2sr.mcp == 0 & bet == 0)
	fp <- sum(bet.2sr.mcp != 0 & bet == 0); fn <- sum(bet.2sr.mcp == 0 & bet != 0)
	mcc.2sr.mcp[i] <- (tp*tn - fp*fn)/sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))
}

#sink(paste("result/", fn, ".out", sep = ""))
cat("PLS-Lasso\n")
cat("Est1\t", mean(est1.pls.lasso), " (", sd(est1.pls.lasso), ")\n", sep = "")
cat("Est2\t", mean(est2.pls.lasso), " (", sd(est2.pls.lasso), ")\n", sep = "")
cat("Pred\t", mean(pred.pls.lasso), " (", sd(pred.pls.lasso), ")\n", sep = "")
cat("TP  \t", mean(tp.pls.lasso), " (", sd(tp.pls.lasso), ")\n", sep = "")
cat("Size\t", mean(size.pls.lasso), " (", sd(size.pls.lasso), ")\n", sep = "")
cat("MCC \t", mean(mcc.pls.lasso), " (", sd(mcc.pls.lasso), ")\n", sep = "")

cat("PLS-SCAD\n")
cat("Est1\t", mean(est1.pls.scad), " (", sd(est1.pls.scad), ")\n", sep = "")
cat("Est2\t", mean(est2.pls.scad), " (", sd(est2.pls.scad), ")\n", sep = "")
cat("Pred\t", mean(pred.pls.scad), " (", sd(pred.pls.scad), ")\n", sep = "")
cat("TP  \t", mean(tp.pls.scad), " (", sd(tp.pls.scad), ")\n", sep = "")
cat("Size\t", mean(size.pls.scad), " (", sd(size.pls.scad), ")\n", sep = "")
cat("MCC \t", mean(mcc.pls.scad), " (", sd(mcc.pls.scad), ")\n", sep = "")

cat("PLS-MCP\n")
cat("Est1\t", mean(est1.pls.mcp), " (", sd(est1.pls.mcp), ")\n", sep = "")
cat("Est2\t", mean(est2.pls.mcp), " (", sd(est2.pls.mcp), ")\n", sep = "")
cat("Pred\t", mean(pred.pls.mcp), " (", sd(pred.pls.mcp), ")\n", sep = "")
cat("TP  \t", mean(tp.pls.mcp), " (", sd(tp.pls.mcp), ")\n", sep = "")
cat("Size\t", mean(size.pls.mcp), " (", sd(size.pls.mcp), ")\n", sep = "")
cat("MCC \t", mean(mcc.pls.mcp), " (", sd(mcc.pls.mcp), ")\n", sep = "")

cat("PLS-oracle\n")
cat("Est1\t", mean(est1.pls.oracle), " (", sd(est1.pls.oracle), ")\n", sep = "")
cat("Est2\t", mean(est2.pls.oracle), " (", sd(est2.pls.oracle), ")\n", sep = "")
cat("Pred\t", mean(pred.pls.oracle), " (", sd(pred.pls.oracle), ")\n", sep = "")

cat("2SR-Lasso\n")
cat("Est1\t", mean(est1.2sr.lasso), " (", sd(est1.2sr.lasso), ")\n", sep = "")
cat("Est2\t", mean(est2.2sr.lasso), " (", sd(est2.2sr.lasso), ")\n", sep = "")
cat("Pred\t", mean(pred.2sr.lasso), " (", sd(pred.2sr.lasso), ")\n", sep = "")
cat("TP  \t", mean(tp.2sr.lasso), " (", sd(tp.2sr.lasso), ")\n", sep = "")
cat("Size\t", mean(size.2sr.lasso), " (", sd(size.2sr.lasso), ")\n", sep = "")
cat("MCC \t", mean(mcc.2sr.lasso), " (", sd(mcc.2sr.lasso), ")\n", sep = "")

cat("2SR-SCAD\n")
cat("Est1\t", mean(est1.2sr.scad), " (", sd(est1.2sr.scad), ")\n", sep = "")
cat("Est2\t", mean(est2.2sr.scad), " (", sd(est2.2sr.scad), ")\n", sep = "")
cat("Pred\t", mean(pred.2sr.scad), " (", sd(pred.2sr.scad), ")\n", sep = "")
cat("TP  \t", mean(tp.2sr.scad), " (", sd(tp.2sr.scad), ")\n", sep = "")
cat("Size\t", mean(size.2sr.scad), " (", sd(size.2sr.scad), ")\n", sep = "")
cat("MCC \t", mean(mcc.2sr.scad), " (", sd(mcc.2sr.scad), ")\n", sep = "")

cat("2SR-MCP\n")
cat("Est1\t", mean(est1.2sr.mcp), " (", sd(est1.2sr.mcp), ")\n", sep = "")
cat("Est2\t", mean(est2.2sr.mcp), " (", sd(est2.2sr.mcp), ")\n", sep = "")
cat("Pred\t", mean(pred.2sr.mcp), " (", sd(pred.2sr.mcp), ")\n", sep = "")
cat("TP  \t", mean(tp.2sr.mcp), " (", sd(tp.2sr.mcp), ")\n", sep = "")
cat("Size\t", mean(size.2sr.mcp), " (", sd(size.2sr.mcp), ")\n", sep = "")
cat("MCC \t", mean(mcc.2sr.mcp), " (", sd(mcc.2sr.mcp), ")\n", sep = "")

cat("2SR-oracle\n")
cat("Est1\t", mean(est1.2sr.oracle), " (", sd(est1.2sr.oracle), ")\n", sep = "")
cat("Est2\t", mean(est2.2sr.oracle), " (", sd(est2.2sr.oracle), ")\n", sep = "")
cat("Pred\t", mean(pred.2sr.oracle), " (", sd(pred.2sr.oracle), ")\n", sep = "")
#sink()

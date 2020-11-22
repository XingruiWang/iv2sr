nseq <- seq(200, 1500, 100)
p <- 100; q <- 100
nn <- length(nseq)

est1s.pls.oracle <- numeric(nn)
preds.pls.oracle <- numeric(nn)

est1s.pls.lasso <- numeric(nn)
preds.pls.lasso <- numeric(nn)
mccs.pls.lasso <- numeric(nn)

est1s.pls.scad <- numeric(nn)
preds.pls.scad <- numeric(nn)
mccs.pls.scad <- numeric(nn)

est1s.pls.mcp <- numeric(nn)
preds.pls.mcp <- numeric(nn)
mccs.pls.mcp <- numeric(nn)

est1s.2sr.oracle <- numeric(nn)
preds.2sr.oracle <- numeric(nn)

est1s.2sr.lasso <- numeric(nn)
preds.2sr.lasso <- numeric(nn)
mccs.2sr.lasso <- numeric(nn)

est1s.2sr.scad <- numeric(nn)
preds.2sr.scad <- numeric(nn)
mccs.2sr.scad <- numeric(nn)

est1s.2sr.mcp <- numeric(nn)
preds.2sr.mcp <- numeric(nn)
mccs.2sr.mcp <- numeric(nn)

for (i in 1:nn) {
	load(paste("data/s", paste(nseq[i], p, q, sep="_"), ".rda", sep=""))
	load(paste("result/s", paste(nseq[i], p, q, sep="_"), ".rda", sep=""))
	n <- nseq[i]; nsim <- length(data)
	
	est1.pls.oracle <- numeric(nsim)
	pred.pls.oracle <- numeric(nsim)
	
	est1.pls.lasso <- numeric(nsim)
	pred.pls.lasso <- numeric(nsim)
	mcc.pls.lasso <- numeric(nsim)
	
	est1.pls.scad <- numeric(nsim)
	pred.pls.scad <- numeric(nsim)
	mcc.pls.scad <- numeric(nsim)
	
	est1.pls.mcp <- numeric(nsim)
	pred.pls.mcp <- numeric(nsim)
	mcc.pls.mcp <- numeric(nsim)
	
	est1.2sr.oracle <- numeric(nsim)
	pred.2sr.oracle <- numeric(nsim)
	
	est1.2sr.lasso <- numeric(nsim)
	pred.2sr.lasso <- numeric(nsim)
	mcc.2sr.lasso <- numeric(nsim)
	
	est1.2sr.scad <- numeric(nsim)
	pred.2sr.scad <- numeric(nsim)
	mcc.2sr.scad <- numeric(nsim)
	
	est1.2sr.mcp <- numeric(nsim)
	pred.2sr.mcp <- numeric(nsim)
	mcc.2sr.mcp <- numeric(nsim)
	
	for (j in 1:nsim) {
		x <- data[[j]]$x; bet <- data[[j]]$bet
		
		bet.pls.oracle <- result[[j]]$bet.pls.oracle
		bet.pls.lasso <- result[[j]]$bet.pls.lasso
		bet.pls.scad <- result[[j]]$bet.pls.scad
		bet.pls.mcp <- result[[j]]$bet.pls.mcp
		
		bet.2sr.oracle <- result[[j]]$bet.2sr.oracle
		bet.2sr.lasso <- result[[j]]$bet.2sr.lasso
		bet.2sr.scad <- result[[j]]$bet.2sr.scad
		bet.2sr.mcp <- result[[j]]$bet.2sr.mcp
		
		est1.pls.oracle[j] <- sum(abs(bet.pls.oracle - bet))
		pred.pls.oracle[j] <- sqrt(sum((x %*% (bet.pls.oracle - bet))^2)/n)
		
		est1.pls.lasso[j] <- sum(abs(bet.pls.lasso - bet))
		pred.pls.lasso[j] <- sqrt(sum((x %*% (bet.pls.lasso - bet))^2)/n)
		tp <- sum(bet.pls.lasso != 0 & bet != 0); tn <- sum(bet.pls.lasso == 0 & bet == 0)
		fp <- sum(bet.pls.lasso != 0 & bet == 0); fn <- sum(bet.pls.lasso == 0 & bet != 0)
		mcc.pls.lasso[j] <- (tp*tn - fp*fn)/sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))
		
		est1.pls.scad[j] <- sum(abs(bet.pls.scad - bet))
		pred.pls.scad[j] <- sqrt(sum((x %*% (bet.pls.scad - bet))^2)/n)
		tp <- sum(bet.pls.scad != 0 & bet != 0); tn <- sum(bet.pls.scad == 0 & bet == 0)
		fp <- sum(bet.pls.scad != 0 & bet == 0); fn <- sum(bet.pls.scad == 0 & bet != 0)
		mcc.pls.scad[j] <- (tp*tn - fp*fn)/sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))
		
		est1.pls.mcp[j] <- sum(abs(bet.pls.mcp - bet))
		pred.pls.mcp[j] <- sqrt(sum((x %*% (bet.pls.mcp - bet))^2)/n)
		tp <- sum(bet.pls.mcp != 0 & bet != 0); tn <- sum(bet.pls.mcp == 0 & bet == 0)
		fp <- sum(bet.pls.mcp != 0 & bet == 0); fn <- sum(bet.pls.mcp == 0 & bet != 0)
		mcc.pls.mcp[j] <- (tp*tn - fp*fn)/sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))
		
		est1.2sr.oracle[j] <- sum(abs(bet.2sr.oracle - bet))
		pred.2sr.oracle[j] <- sqrt(sum((x %*% (bet.2sr.oracle - bet))^2)/n)
		
		est1.2sr.lasso[j] <- sum(abs(bet.2sr.lasso - bet))
		pred.2sr.lasso[j] <- sqrt(sum((x %*% (bet.2sr.lasso - bet))^2)/n)
		tp <- sum(bet.2sr.lasso != 0 & bet != 0); tn <- sum(bet.2sr.lasso == 0 & bet == 0)
		fp <- sum(bet.2sr.lasso != 0 & bet == 0); fn <- sum(bet.2sr.lasso == 0 & bet != 0)
		mcc.2sr.lasso[j] <- (tp*tn - fp*fn)/sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))
		
		est1.2sr.scad[j] <- sum(abs(bet.2sr.scad - bet))
		pred.2sr.scad[j] <- sqrt(sum((x %*% (bet.2sr.scad - bet))^2)/n)
		tp <- sum(bet.2sr.scad != 0 & bet != 0); tn <- sum(bet.2sr.scad == 0 & bet == 0)
		fp <- sum(bet.2sr.scad != 0 & bet == 0); fn <- sum(bet.2sr.scad == 0 & bet != 0)
		mcc.2sr.scad[j] <- (tp*tn - fp*fn)/sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))
		
		est1.2sr.mcp[j] <- sum(abs(bet.2sr.mcp - bet))
		pred.2sr.mcp[j] <- sqrt(sum((x %*% (bet.2sr.mcp - bet))^2)/n)
		tp <- sum(bet.2sr.mcp != 0 & bet != 0); tn <- sum(bet.2sr.mcp == 0 & bet == 0)
		fp <- sum(bet.2sr.mcp != 0 & bet == 0); fn <- sum(bet.2sr.mcp == 0 & bet != 0)
		mcc.2sr.mcp[j] <- (tp*tn - fp*fn)/sqrt((tp + fp)*(tp + fn)*(tn + fp)*(tn + fn))
		if (is.nan(mcc.2sr.mcp[j])) mcc.2sr.mcp[j] <- 0
	}
	
	est1s.pls.oracle[i] <- mean(est1.pls.oracle)
	preds.pls.oracle[i] <- mean(pred.pls.oracle)
	
	est1s.pls.lasso[i] <- mean(est1.pls.lasso)
	preds.pls.lasso[i] <- mean(pred.pls.lasso)
	mccs.pls.lasso[i] <- mean(mcc.pls.lasso)
	
	est1s.pls.scad[i] <- mean(est1.pls.scad)
	preds.pls.scad[i] <- mean(pred.pls.scad)
	mccs.pls.scad[i] <- mean(mcc.pls.scad)
	
	est1s.pls.mcp[i] <- mean(est1.pls.mcp)
	preds.pls.mcp[i] <- mean(pred.pls.mcp)
	mccs.pls.mcp[i] <- mean(mcc.pls.mcp)
	
	est1s.2sr.oracle[i] <- mean(est1.2sr.oracle)
	preds.2sr.oracle[i] <- mean(pred.2sr.oracle)
	
	est1s.2sr.lasso[i] <- mean(est1.2sr.lasso)
	preds.2sr.lasso[i] <- mean(pred.2sr.lasso)
	mccs.2sr.lasso[i] <- mean(mcc.2sr.lasso)
	
	est1s.2sr.scad[i] <- mean(est1.2sr.scad)
	preds.2sr.scad[i] <- mean(pred.2sr.scad)
	mccs.2sr.scad[i] <- mean(mcc.2sr.scad)
	
	est1s.2sr.mcp[i] <- mean(est1.2sr.mcp)
	preds.2sr.mcp[i] <- mean(pred.2sr.mcp)
	mccs.2sr.mcp[i] <- mean(mcc.2sr.mcp)
}

save(est1s.pls.oracle, est1s.pls.lasso, est1s.pls.scad, est1s.pls.mcp,
	 est1s.2sr.oracle, est1s.2sr.lasso, est1s.2sr.scad, est1s.2sr.mcp,
	 preds.pls.oracle, preds.pls.lasso, preds.pls.scad, preds.pls.mcp,
	 preds.2sr.oracle, preds.2sr.lasso, preds.2sr.scad, preds.2sr.mcp,
	 mccs.pls.lasso, mccs.pls.scad, mccs.pls.mcp, mccs.2sr.lasso, mccs.2sr.scad, mccs.2sr.mcp,
	 file="result/plot_low.rda")

pdf.options(width=5.8, height=4.7)
setEPS(width=5.8, height=4.7)
pdf.out <- TRUE

if (pdf.out) pdf(file="fig/low_est1_pls.pdf") else postscript(file="fig/low_est1_pls.eps")
par(mai=c(0.7, 0.8, 0.4, 0.2), mgp=c(2.1, 0.7, 0))
plot(nseq, est1s.pls.lasso, type="n", ylim=c(0, 4), main="PLS", cex.main=1.2, xlab="Sample size",
	 ylab=expression(paste(italic(L)[1], " estimation loss")), cex.lab=1.2)
lines(nseq, est1s.pls.lasso, lwd=1.6)
lines(nseq, est1s.pls.scad, lty=2, lwd=1.6)
lines(nseq, est1s.pls.mcp, lty=3, lwd=1.6)
lines(nseq, est1s.pls.oracle, col="red", lty=4, lwd=1.6)
legend("topright", c("Lasso", "SCAD", "MCP", "Oracle"), col=c("black", "black", "black", "red"),
	   lty=1:4, lwd=1.6)
dev.off()

if (pdf.out) pdf(file="fig/low_est1_2sr.pdf") else postscript(file="fig/low_est1_2sr.eps")
par(mai=c(0.7, 0.8, 0.4, 0.2), mgp=c(2.1, 0.7, 0))
plot(nseq, est1s.2sr.lasso, type="n", ylim=c(0, 4), main="2SR", cex.main=1.2, xlab="Sample size",
	 ylab=expression(paste(italic(L)[1], " estimation loss")), cex.lab=1.2)
lines(nseq, est1s.2sr.lasso, lwd=1.6)
lines(nseq, est1s.2sr.scad, lty=2, lwd=1.6)
lines(nseq, est1s.2sr.mcp, lty=3, lwd=1.6)
lines(nseq, est1s.2sr.oracle, col="red", lty=4, lwd=1.6)
legend("topright", c("Lasso", "SCAD", "MCP", "Oracle"), col=c("black", "black", "black", "red"),
	   lty=1:4, lwd=1.6)
dev.off()

if (pdf.out) pdf(file="fig/low_pred_pls.pdf") else postscript(file="fig/low_pred_pls.eps")
par(mai=c(0.7, 0.8, 0.4, 0.2), mgp=c(2.1, 0.7, 0))
plot(nseq, preds.pls.lasso, type="n", ylim=c(0, 1.5), main="PLS", cex.main=1.2, xlab="Sample size",
	 ylab="Prediction loss", cex.lab=1.2)
lines(nseq, preds.pls.lasso, lwd=1.6)
lines(nseq, preds.pls.scad, lty=2, lwd=1.6)
lines(nseq, preds.pls.mcp, lty=3, lwd=1.6)
lines(nseq, preds.pls.oracle, col="red", lty=4, lwd=1.6)
legend("topright", c("Lasso", "SCAD", "MCP", "Oracle"), col=c("black", "black", "black", "red"),
	   lty=1:4, lwd=1.6)
dev.off()

if (pdf.out) pdf(file="fig/low_pred_2sr.pdf") else postscript(file="fig/low_pred_2sr.eps")
par(mai=c(0.7, 0.8, 0.4, 0.2), mgp=c(2.1, 0.7, 0))
plot(nseq, preds.2sr.lasso, type="n", ylim=c(0, 1.5), main="2SR", cex.main=1.2, xlab="Sample size",
	 ylab="Prediction loss", cex.lab=1.2)
lines(nseq, preds.2sr.lasso, lwd=1.6)
lines(nseq, preds.2sr.scad, lty=2, lwd=1.6)
lines(nseq, preds.2sr.mcp, lty=3, lwd=1.6)
lines(nseq, preds.2sr.oracle, col="red", lty=4, lwd=1.6)
legend("topright", c("Lasso", "SCAD", "MCP", "Oracle"), col=c("black", "black", "black", "red"),
	   lty=1:4, lwd=1.6)
dev.off()

if (pdf.out) pdf(file="fig/low_mcc_pls.pdf") else postscript(file="fig/low_mcc_pls.eps")
par(mai=c(0.7, 0.8, 0.4, 0.2), mgp=c(2.1, 0.7, 0))
plot(nseq, mccs.pls.lasso, type="n", ylim=c(0, 1), main="PLS", cex.main=1.2, xlab="Sample size",
	 ylab="MCC", cex.lab=1.2)
lines(nseq, mccs.pls.lasso, lwd=1.6)
lines(nseq, mccs.pls.scad, lty=2, lwd=1.6)
lines(nseq, mccs.pls.mcp, lty=3, lwd=1.6)
legend("topright", c("Lasso", "SCAD", "MCP"), lty=1:3, lwd=1.6)
dev.off()

if (pdf.out) pdf(file="fig/low_mcc_2sr.pdf") else postscript(file="fig/low_mcc_2sr.eps")
par(mai=c(0.7, 0.8, 0.4, 0.2), mgp=c(2.1, 0.7, 0))
plot(nseq, mccs.2sr.lasso, type="n", ylim=c(0, 1), main="2SR", cex.main=1.2, xlab="Sample size",
	 ylab="MCC", cex.lab=1.2)
lines(nseq, mccs.2sr.lasso, lwd=1.6)
lines(nseq, mccs.2sr.scad, lty=2, lwd=1.6)
lines(nseq, mccs.2sr.mcp, lty=3, lwd=1.6)
legend("bottomright", c("Lasso", "SCAD", "MCP"), lty=1:3, lwd=1.6)
dev.off()

#' Description of the file
#'
#' @param x x
#' @param y y
#' @return return
#' @examples
#' cd()

#' @export
cd <- function(y, x, method, lam, nlam=50, a=3.7, upto=0.4*length(y), maxit=50, tol=1e-4) {
	pen <- as.integer(switch(method, Lasso=0, SCAD=1, MCP=2))
	y <- scale(y, scale=FALSE)
	x <- scale(x, scale=apply(x, 2, sd)*sqrt(nrow(x)-1))
	if (missing(lam)) {
		lam.max <- max(abs(crossprod(x, y)))
		lam <- lam.max*exp(seq(0, -log(nlam), length=nlam))
	}
	#' @useDynLib iv2sr coordinate_descent
	ans <- .Call("coordinate_descent", y, x, pen, lam, a, upto, as.integer(maxit), tol)
	names(ans) <- c("sol", "lam")
	ans$sol <- ans$sol/attr(x, "scaled:scale")
	ans
}

#' @export
cv <- function(y, x, method, lam, ind, nfold=10) {
	n <- length(y); nlam <- length(lam)
	if (missing(ind)) ind <- split(1:n, sample(rep(1:nfold, length=n)))
	pred <- matrix(, nfold, nlam)
	for (i in 1:nfold) {
		ind2 <- ind[[i]]; ind1 <- setdiff(1:n, ind2)
		y1 <- y[ind1]; y2 <- y[ind2]
		x1 <- x[ind1, ]; x2 <- x[ind2, ]
		sol <- cd(y1, x1, method, lam, upto=Inf)$sol
		pred[i, ] <- colSums((y2 - x2 %*% sol)^2)
	}
	cv.err <- colMeans(pred)
	se <- apply(pred, 2, sd)
	i <- which.min(cv.err)
	list(cv.err=cv.err, se=se, i=i)
}

#' @export
stab <- function(y, x, method, lam, nsample=100, nsub=floor(0.5*length(y)), seed=41) {
	set.seed(seed)
	n <- length(y); p <- ncol(x); nlam <- length(lam)
	path <- matrix(0, p, nlam)
	for (i in 1:nsample) {
		ind <- sample(1:n, nsub)
		y1 <- y[ind]; x1 <- x[ind, ]
		sol <- cd(y1, x1, method, lam, upto=Inf)$sol
		path <- path + (sol != 0)
	}
	path <- path/nsample
	prob <- apply(path, 1, max)
	list(path=path, prob=prob)
}

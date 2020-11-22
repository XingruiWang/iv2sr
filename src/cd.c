#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

SEXP coordinate_descent(SEXP y, SEXP x, SEXP method, SEXP lam_seq, SEXP shape, SEXP upto,
						SEXP max_it, SEXP tol) {
	R_len_t	n = length(y), p = ncols(x), n_lam = length(lam_seq), i, j, k;
	int		pen = INTEGER(method)[0], max_iter = INTEGER(max_it)[0], iter, brk, s;
	double	*r_y = REAL(y), *r_x = REAL(x), *r_lam_seq = REAL(lam_seq), a = REAL(shape)[0],
			a1 = NA_REAL, a2 = NA_REAL, max_s = REAL(upto)[0], eps = REAL(tol)[0], *sol, *bet,
			*res, t0, a_t0, lam, del, old, dif, norm, norm_old;
	SEXP	ans;
	
	sol = (double *) R_alloc(p*n_lam, sizeof(double));
	bet = (double *) R_alloc(p, sizeof(double));
	res = (double *) R_alloc(n, sizeof(double));
	
	for (i = 0; i < p; i++) bet[i] = 0.0;
	for (i = 0; i < n; i++) res[i] = r_y[i];
	switch(pen) {
		case 1:		/* SCAD */
			a1 = a/(a - 1.0); a2 = (a - 1.0)/(a - 2.0); break;
		case 2:		/* MCP */
			a1 = a/(a - 1.0); break;
	}
	brk = 0; norm = 0.0;
	for (j = 0; j < n_lam && !brk; j++) {
		lam = r_lam_seq[j];
		for (iter = 0; iter < max_iter; iter++) {
			del = 0.0; norm_old = norm; norm = 0.0;
			for (i = 0; i < p; i++) {
				old = bet[i]; t0 = bet[i];
				for (k = 0; k < n; k++) t0 += res[k]*r_x[k + i*n];
				a_t0 = fabs(t0);
				switch (pen) {
					case 0:		/* Lasso */
						bet[i] = copysign(fdim(a_t0, lam), t0); break;
					case 1:		/* SCAD */
						if (a_t0 > a*lam)
							bet[i] = t0;
						else if (a_t0 > 2.0*lam)
							bet[i] = copysign((a_t0 - a1*lam)*a2, t0);
						else
							bet[i] = copysign(fdim(a_t0, lam), t0);
						break;
					case 2:		/* MCP */
						if (a_t0 > a*lam)
							bet[i] = t0;
						else
							bet[i] = copysign(fdim(a_t0, lam)*a1, t0);
						break;
				}
				norm += fabs(bet[i]);
				if (bet[i] != old) {
					dif = bet[i] - old;
					for (k = 0; k < n; k++) res[k] -= r_x[k + i*n]*dif;
					del += fabs(dif);
				}
			}
			if (del < eps*norm_old) break;
		}
		s = 0;
		for (i = 0; i < p; i++) {
			sol[i + j*p] = bet[i];
			if (bet[i] != 0.0) s++;
		}
		if (R_FINITE(max_s) && s > max_s) brk = 1;
	}
	
	PROTECT(ans = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(ans, 0, allocMatrix(REALSXP, p, j));
	for (i = 0; i < p*j; i++) REAL(VECTOR_ELT(ans, 0))[i] = sol[i];
	SET_VECTOR_ELT(ans, 1, allocVector(REALSXP, j));
	for (i = 0; i < j; i++) REAL(VECTOR_ELT(ans, 1))[i] = r_lam_seq[i];
	
	UNPROTECT(1);
	return ans;
}

library(qtl2pleio)

context("testing parametric bootstrap code")



probs_pre <- rbinom(n = 100 * 10, size = 1, prob = 1 / 2)
probs <- array(data = probs_pre, dim = c(100, 1, 10))

X <- as.matrix(probs[ , , 5])

pvl_boot(X = X, B = c(1, 2), Vg_initial = diag(2), Ve_initial = diag(2), kinship = diag(100), probs = probs, start_snp = 1, n_snp = 10, nboot = 5)

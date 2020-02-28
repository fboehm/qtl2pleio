## ----setup2-------------------------------------------------------------------
is_inst <- function(pkg) {
    nzchar(system.file(package = pkg))
}
qtl2_indic <- is_inst("qtl2")

## ----load-qtl2pleio, eval = qtl2_indic----------------------------------------
library(qtl2pleio)
library(qtl2)

## ----download-DOex, eval = qtl2_indic-----------------------------------------
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/DOex/DOex.zip")
DOex <- read_cross2(file)

## ----calc_pr, eval = qtl2_indic-----------------------------------------------
probs <- calc_genoprob(DOex)
pr <- genoprob_to_alleleprob(probs)

## ----calc-K, eval = qtl2_indic------------------------------------------------
calc_kinship(probs = pr, type = "loco") -> kinship

## ----pp-def, eval = qtl2_indic------------------------------------------------
# set up the design matrix, X
pp <- pr[[2]] #we'll work with Chr 3's genotype data
dim(pp)

## ----X-def, eval = qtl2_indic-------------------------------------------------
#Next, we prepare a design matrix X
X <- gemma2::stagger_mats(pp[ , , 50], pp[ , , 50], pp[ , , 50])
dim(X)

## ----B-def, eval = qtl2_indic-------------------------------------------------
# assemble B matrix of allele effects
B <- matrix(data = rep(rep(c(-1, 1), each = 4), times = 3), nrow = 8, ncol = 3, byrow = FALSE)
# verify that B is what we want:
B
# set.seed to ensure reproducibility
set.seed(2018-01-30)
# call to sim1
Ypre <- sim1(X = X, B = B, Vg = diag(3), Ve = diag(3), kinship = kinship[[2]])
Y <- matrix(Ypre, nrow = 261, ncol = 3, byrow = FALSE)
rownames(Y) <- rownames(pp)
colnames(Y) <- c("tr1", "tr2", "tr3")

## ----1d-scans, eval = qtl2_indic----------------------------------------------
s1 <- scan1(genoprobs = pr, pheno = Y, kinship = kinship)

## ----1d-lod-plots, eval = qtl2_indic------------------------------------------
plot(s1, DOex$pmap$`3`)
plot(s1, DOex$pmap$`3`, lod = 2, col ="violetred", add=TRUE)
plot(s1, DOex$pmap$`3`, lod = 3, col ="green", add=TRUE)
legend("topleft", colnames(s1), lwd = 2, col=c("darkslateblue", "violetred", "green"), bg="gray92")


## ----find-peaks, eval = qtl2_indic--------------------------------------------
find_peaks(s1, map = DOex$pmap, threshold = 8)

## ----scan, eval = qtl2_indic--------------------------------------------------
start_index <- 43
out <- scan_pvl(probs = pp,
                pheno = Y,
                kinship = kinship$`3`,
                start_snp = start_index,
                n_snp = 15, n_cores = 1
                )

## ----detect-cores, eval = FALSE-----------------------------------------------
#  parallel::detectCores()

## ----check-out, eval = qtl2_indic---------------------------------------------
out

## ----profile-plot, eval = qtl2_indic------------------------------------------
library(dplyr)
library(ggplot2)
out_lods <- out %>%
  calc_profile_lods() %>%
  add_pmap(pmap = DOex$pmap$`3`)
out_lods %>% 
  ggplot() + geom_line(aes(x = marker_position, y = profile_lod, colour = trait))

## ----lrt-calc, eval = qtl2_indic----------------------------------------------
(lrt <- max(out_lods$profile_lod))

## ----get-pleio-index, eval = qtl2_indic---------------------------------------
(pleio_index <- find_pleio_peak_tib(out, start_snp = start_index))

## ----boot, eval = qtl2_indic--------------------------------------------------
set.seed(2018 - 11 - 25)
suppressMessages(b_out <- boot_pvl(probs = pp,
                                   pheno = Y,
                                   pleio_peak_index = pleio_index,
                                   kinship = kinship$`3`, 
                                   nboot_per_job = 2,
                                   start_snp = start_index,
                                   n_snp = 15
                                   )
                 )


## ----pval, eval = qtl2_indic--------------------------------------------------
b_out
(pvalue <- mean(b_out >= lrt))

## -----------------------------------------------------------------------------
devtools::session_info()


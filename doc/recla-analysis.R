## ----setup2-------------------------------------------------------------------
is_inst <- function(pkg) {
    nzchar(system.file(package = pkg))
}
qtl2_indic <- is_inst("qtl2")
knitr::opts_chunk$set(eval = qtl2_indic)

## ----pkgs---------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(qtl2pleio)
library(qtl2)

## ----load-recla---------------------------------------------------------------
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/DO_Recla/recla.zip")
recla <- read_cross2(file)
# make sex a covariate for use in qtl2pleio::scan_pvl
recla[[6]][ , 1, drop = FALSE] -> sex
# insert pseudomarkers
insert_pseudomarkers(recla, step = 0.10) -> pseudomap
gm <- pseudomap$`8`

## ----calc-genoprobs-code------------------------------------------------------
probs <- calc_genoprob(recla, map = pseudomap, cores = 1)

## ----calc-aprobs--------------------------------------------------------------
aprobs <- genoprob_to_alleleprob(probs)

## ----calc-kinship-------------------------------------------------------------
kinship <- calc_kinship(aprobs, "loco")

## ----log-phenos---------------------------------------------------------------
recla$pheno -> ph
log(ph) -> lph
apply(FUN = broman::winsorize, X = lph, MARGIN = 2) -> wlph
as_tibble(wlph) -> wlph_tib

## ----scan1--------------------------------------------------------------------
sex2 <- matrix(as.numeric(sex == "female"), ncol = 1)
colnames(sex2) <- "female"
rownames(sex2) <- rownames(aprobs[[1]])
out <- scan1(genoprobs = aprobs, 
             pheno = wlph, 
             kinship = kinship, 
             addcovar = sex2, 
             reml = TRUE
             )

## ----get-peaks----------------------------------------------------------------
(peaks <- find_peaks(out, pseudomap, threshold = 5) %>%
  dplyr::arrange(chr, pos) %>%
   dplyr::select(- lodindex))
peaks8 <- peaks %>%
  dplyr::filter(chr == 8, pos > 50, pos < 60)
pos_LD_light_pct <- peaks8 %>%
  dplyr::filter(lodcolumn == "LD_light_pct") %>%
  dplyr::select(pos)
pos_HP_latency <- peaks8 %>%
  dplyr::filter(lodcolumn == "HP_latency") %>%
  dplyr::select(pos)

## ----cors---------------------------------------------------------------------
cor(wlph[ , 7], wlph[ , 10], use = "complete.obs")
cor(wlph[ , 22], wlph[ , 10], use = "complete.obs")
cor(wlph[ , 7], wlph[ , 22], use = "complete.obs")

## ----scatter------------------------------------------------------------------
ggplot() + 
  geom_point(data = wlph_tib, aes(y = HP_latency, x = LD_light_pct)) +
  labs(x = "Percent time in light", y = "Hot plate latency")

## ----lod10-plot---------------------------------------------------------------
plot(out, map = pseudomap, 
     lodcolumn = 10, 
     main = "percent time in light"
     )

## ----lod22-plot---------------------------------------------------------------
plot(out, map = pseudomap, 
     lodcolumn = 22, 
     main = "hot plate latency"
     )

## ----coefs-calc---------------------------------------------------------------
scan1coef(aprobs[ , 8], pheno = wlph[ , 10], kinship = kinship$`8`, 
          reml = TRUE,
          addcovar = sex2) -> s1c_10
scan1coef(aprobs[ , 8], pheno = wlph[ , 22], kinship = kinship$`8`, 
          reml = TRUE,
          addcovar = sex2) -> s1c_22

## ----coefs-subset-------------------------------------------------------------
# subset scan1output objects
s1c_10s <- s1c_10[650:999, ] 
# 650:999 is the same as the interval for the two-dimensional scan.
s1c_22s <- s1c_22[650:999, ]

## ----plot-coefs---------------------------------------------------------------
plot_coefCC(s1c_10s, scan1_output = out[ , 10, drop = FALSE], map = pseudomap, main = "percent time in light")
plot_coefCC(s1c_22s, scan1_output = out[ , 22, drop = FALSE], map = pseudomap, main = "hot plate latency")

## ----2d-scan, eval = FALSE----------------------------------------------------
#  scan_pvl(probs = aprobs$`8`,
#           pheno = wlph[, c(10, 22)],
#           addcovar = sex2,
#           kinship = kinship$`8`,
#           start_snp = 650,
#           n_snp = 350) -> pvl1022

## ----2d-scan-results, echo = TRUE---------------------------------------------
as_tibble(read.table("https://zenodo.org/record/3210710/files/recla-10-22.txt?download=1", stringsAsFactors = FALSE)) -> pvl1022

## ----lrt-calc-----------------------------------------------------------------
(mylrt <- calc_lrt_tib(pvl1022))

## ----profile-plot-------------------------------------------------------------
colnames(recla$pheno)[c(10, 22)] <- c("Percent time in light", "Hot plate latency")
pvl1022 %>%
  mutate(log10lik = loglik / log(10)) %>%
  dplyr::select(- loglik) %>%
  calc_profile_lods() %>%
  add_pmap(pmap = recla$pmap$`8`) %>%
  ggplot() + geom_line(aes(x = marker_position, y = profile_lod, colour = trait))

## ----get-pleio-peak-----------------------------------------------------------
find_pleio_peak_tib(pvl1022, start_snp = 650)

## ----read-boot----------------------------------------------------------------
## read boot lrt files
boot_lrt <- list()
for (i in 1:1000){
  n <- i - 1
  fn <- paste0("https://raw.githubusercontent.com/fboehm/qtl2pleio-manuscript-chtc/master/Recla-bootstrap/results/recla-boot-run561_", n, ".txt")
  boot_lrt[i] <- read.table(fn)
}
# convert list to numeric vector
boot_lrt <- unlist(boot_lrt)

## ----pval---------------------------------------------------------------------
sum(boot_lrt >= mylrt) / length(boot_lrt)

## ----session-info, eval = TRUE------------------------------------------------
devtools::session_info()


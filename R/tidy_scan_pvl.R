#' Assemble tibble from matrix of log-likelihood values
#'
#' @family profile log-likelihood tibble functions
#' @param loglik_mat a square matrix of log-likelihood values with rownames and colnames matching
#' @return tibble with 3 columns: marker1, marker2 and log-likelihood
#' @examples
#' llmat <- matrix(nrow = 3, ncol = 3, data = rgamma(9, 5))
#' rownames(llmat) <- paste0('m', 1:3)
#' colnames(llmat) <- paste0('m', 1:3)
#' transform_loglik_mat(llmat)
#' @export
transform_loglik_mat <- function(loglik_mat) {
    marker1 <- rep(rownames(loglik_mat), times = ncol(loglik_mat))
    marker2 <- rep(colnames(loglik_mat), each = nrow(loglik_mat))
    ll <- as.vector(loglik_mat)
    dat <- tibble::tibble(marker1, marker2, ll)
    return(dat)
}

#' Add physical map contents to tibble
#'
#' @family profile log-likelihood tibble functions
#' @param tib a tibble with 3 columns: marker1 name, marker2 name and log-likelihood values
#' @param pmap a physical map for a single chromosome
#' @return a tibble with 5 columns: marker1, marker2, log-likelihood, marker1_position, marker2_position
#' @examples
#' pm <- 1:3
#' names(pm) <- as.character(paste0('m', 1:3))
#' expand.grid(paste0('m', 1:3), paste0('m', 1:3)) -> foo
#' tib <- tibble::tibble(marker1 = as.character(foo[ , 1]),
#'   marker2 = as.character(foo[ , 2]))
#' tib$ll <- rgamma(9, 5)
#' add_pmap(tib, pm)
#' @export
#' @importFrom rlang .data

add_pmap <- function(tib, pmap) {
    pmap_tib <- tibble::tibble(marker = as.character(names(pmap)), marker_position = pmap)
    out <- pmap_tib %>%
        dplyr::right_join(tib, by = c("marker" = "marker1")) %>%
        dplyr::rename(marker1_position = .data$marker_position, marker1 = .data$marker) %>%
        dplyr::left_join(pmap_tib, by = c("marker2" = "marker")) %>%
        dplyr::rename(marker2_position = .data$marker_position)
    return(out)
}

#' Assemble a profile log-likelihood tibble
#'
#' @family profile log-likelihood tibble functions
#' @param tib a tibble, derived from loglik_mat and with pmap info
#' @param trace character vector to identify which profile log likelihood to calculate. Takes values 'profile1' or 'profile2'
#' @export
#' @return a tibble with 3 columns: marker position, profile log-likelihood, and trace (profile1 or profile2)
#'
#' @importFrom rlang .data
#' @examples
#' pm <- 1:3
#' names(pm) <- as.character(paste0('m', 1:3))
#' expand.grid(paste0('m', 1:3), paste0('m', 1:3)) -> foo
#' tib <- tibble::tibble(marker1 = as.character(foo[ , 1]),
#'   marker2 = as.character(foo[ , 2]))
#' tib$ll <- rgamma(9, 5)
#' add_pmap(tib, pm) -> tib_pm
#' assemble_profile_tib(tib = tib_pm, trace = "profile1")

assemble_profile_tib <- function(tib, trace = "profile1") {
    if (trace == "profile1")
        tib2 <- dplyr::group_by(tib, .data$marker1)
    if (trace == "profile2")
        tib2 <- dplyr::group_by(tib, .data$marker2)
    tib3 <- dplyr::filter(tib2, .data$ll == max(.data$ll))
    tib4 <- dplyr::ungroup(tib3)
    if (trace == "profile1") {
        tib4$trace <- "profile1"
        tib4$marker_position <- tib4$marker1_position
    }
    if (trace == "profile2") {
        tib4$trace <- "profile2"
        tib4$marker_position <- tib4$marker2_position
    }
    tib5 <- dplyr::select(tib4, .data$marker_position, .data$ll, .data$trace)
    return(tib5)
}

#' Tidy the data frame outputted by scan_pvl for further analysis & plotting
#'
#' @family profile log-likelihood tibble functions
#' @param mytib outputted tibble from scan_pvl
#' @param pmap physical map (in Mb) or genetic map (in cM) for exactly one chromosome
#' @export
#' @importFrom rlang .data
#' @examples
#' # read data
#' iron <- qtl2::read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' # insert pseudomarkers into map
#' map <- qtl2::insert_pseudomarkers(iron$gmap, step=1)
#' # calculate genotype probabilities
#' probs <- qtl2::calc_genoprob(iron, map, error_prob=0.002)
#' # grab phenotypes and covariates; ensure that covariates have names attribute
#' pheno <- iron$pheno
#' # leave-one-chromosome-out kinship matrices
#' kinship <- qtl2::calc_kinship(probs, "loco")$`1`
#' # get founder allele probabilites
#' probs <- qtl2::genoprob_to_alleleprob(probs)$`1`
#' ac <- matrix(as.numeric(iron$covar$sex == "m", ncol = 1))
#' colnames(ac) <- "sex"
#' rownames(ac) <- rownames(probs)
#' ss <- scan_pvl(probs = probs, pheno = pheno, kinship = kinship, addcovar = ac,
#' start_snp = 1, n_snp = 80, n_cores = 1)
#' tidy_scan_pvl(ss, pm = map)


tidy_scan_pvl <- function(mytib, pmap) {
    dat <- mytib %>%
        dplyr::rename(marker1 = .data$Var1, marker2 = .data$Var2, ll = .data$loglik) %>%
        add_pmap(pmap)
    pl2 <- dat %>%
        dplyr::filter(.data$marker1 == .data$marker2) %>%
        dplyr::rename(marker_position = .data$marker1_position)
    pl2$trace <- "pleio"
    pleio_ll <- dplyr::select(pl2, .data$marker_position, .data$ll, .data$trace)
    # assemble pro1_ll
    pro1 <- assemble_profile_tib(dat, "profile1")
    pro2 <- assemble_profile_tib(dat, "profile2")
    # bind 3 tibbles
    foo <- dplyr::bind_rows(pleio_ll, pro1, pro2)
    foo$lod <- (foo$ll - max(pleio_ll$ll)) / log(10)  # convert from base e to base 10
    dat <- dplyr::select(foo, .data$marker_position, .data$lod, .data$trace)
    return(dat)
}

#' Add intercepts to tidied log-likelihood tibble
#'
#' @family profile log-likelihood tibble functions
#' @param tib a tibble that results from tidy_scan_pvl acting on a log likelihood matrix
#' @param intercepts_univariate a vector of length 2 that contains the x coordinate values for the ordered univariate peaks
#' @export
#' @importFrom rlang .data
#' @examples
#' # read data
#' iron <- qtl2::read_cross2(system.file("extdata", "iron.zip", package="qtl2"))
#' # insert pseudomarkers into map
#' map <- qtl2::insert_pseudomarkers(iron$gmap, step=1)
#' # calculate genotype probabilities
#' probs <- qtl2::calc_genoprob(iron, map, error_prob=0.002)
#' # grab phenotypes and covariates; ensure that covariates have names attribute
#' pheno <- iron$pheno
#' # leave-one-chromosome-out kinship matrices
#' kinship <- qtl2::calc_kinship(probs, "loco")$`1`
#' # get founder allele probabilites
#' probs <- qtl2::genoprob_to_alleleprob(probs)$`1`
#' ac <- matrix(as.numeric(iron$covar$sex == "m", ncol = 1))
#' colnames(ac) <- "sex"
#' rownames(ac) <- rownames(probs)
#' ss <- scan_pvl(probs = probs, pheno = pheno, kinship = kinship, addcovar = ac,
#' start_snp = 1, n_snp = 80, n_cores = 1)
#' tidy_scan_pvl(ss, pm = map) -> tidy_ss
#' # in practice, use the univariate peaks for each trait.
#' # here, we add arbitrary univariate peaks to illustrate function use.
#' add_intercepts(tib = tidy_ss, intercepts_univariate = c(8, 10))


add_intercepts <- function(tib, intercepts_univariate) {
    tib$intercept_uni <- NA
    tib$intercept_uni[tib$trace == "profile1"] <- intercepts_univariate[1]
    tib$intercept_uni[tib$trace == "profile2"] <- intercepts_univariate[2]
    # determine pleiotropy peak location
    pleio_lod <- dplyr::filter(tib, .data$trace == "pleio")
    pleio_peak <- pleio_lod$marker_position[which.max(pleio_lod$lod)]
    tib$intercept_pleio <- NA
    tib$intercept_pleio[tib$trace == "pleio"] <- pleio_peak
    # determine bivariate peak locations
    tib$intercept_biv <- NA
    profile1_lod <- dplyr::filter(tib, .data$trace == "profile1")
    profile1_peak <- profile1_lod$marker_position[which.max(profile1_lod$lod)]
    profile2_lod <- dplyr::filter(tib, .data$trace == "profile2")
    profile2_peak <- profile2_lod$marker_position[which.max(profile2_lod$lod)]
    tib$intercept_biv[tib$trace == "profile1"] <- profile1_peak
    tib$intercept_biv[tib$trace == "profile2"] <- profile2_peak
    return(tib)
}


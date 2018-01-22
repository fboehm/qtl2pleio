library(qtl2pleio)
context("Testing family of tibble-building helper functions")

marker1 <- rep(c("m1", "m2", "m3"), times = 3)
marker2 <- rep(c("m1", "m2", "m3"), each = 3)
set.seed(2018-01-22)
ll <- rgamma(n = 9, shape = 5)
mytib <- tibble::tibble(marker1, marker2, ll)
# mymat assembly
mymat <- matrix(ll, nrow = 3, ncol = 3, byrow = FALSE)
rownames(mymat) <- c("m1", "m2", "m3")
colnames(mymat) <- c("m1", "m2", "m3")
# define pm
pm <- 1:3
names(pm) <- paste0("m", 1:3)

test_that("transform_loglik_mat works on a small matrix of ll values", {
  expect_equal(mytib, transform_loglik_mat(mymat))
})

mytib -> mytib2
mytib2$marker1_position <- rep(1:3, times = 3)
mytib2$marker2_position <- rep(1:3, each = 3)
test_that("add_pmap works on the small tibble", {
  expect_equal(add_pmap(tib = mytib, pmap = pm), mytib2)
})

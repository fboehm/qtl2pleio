library(qtl2pleio)

context("testing add_pmap")

## example code
pm <- 1:3
names(pm) <- as.character(paste0('m', 1:3))
expand.grid(paste0('m', 1:3), paste0('m', 1:3)) -> foo
tib <- tibble::tibble(marker1 = as.character(foo[ , 1]),
   marker2 = as.character(foo[ , 2]))
tib$ll <- rgamma(9, 5)
pm2 <- 1:4
names(pm2) <- as.character(paste0('m', 1:4))
# tests

test_that("add_pmap returns a tibble with 5 columns", {
  expect_s3_class(add_pmap(tib, pm), class = "tbl_df")
  expect_true(tibble::is_tibble(add_pmap(tib, pm)))
  expect_equal(ncol(add_pmap(tib, pm)), 5)
})

test_that("add_pmap returns object with correct number of rows when inputs differ in length or number of rows",
          {
            expect_equal(nrow(add_pmap(tib, pm2)), 3 ^ 2)
            expect_equal(nrow(add_pmap(tib, pm2[1:2])), 3 ^ 2)
            expect
          })


library("falconnr")

test_that("parameter changes effective", {
    n <- 1000
    d <- 10
    p <- LshParameterSetter$new(n, d)
    pv0 <- p$asList()

    expect_equal(pv0$points, n)
    expect_equal(pv0$dimension, d)
    expect_equal(pv0$lshFamily, "cross_polytope")
    expect_equal(pv0$distance, "euclidean_squared")
    expect_equal(pv0$storage, "bit_packed_flat_hash_table")

    pv1 <- p$distance("negative_inner_product")$family("hyperplane")$asList()
    expect_equal(pv1$lshFamily, "hyperplane")
    expect_equal(pv1$distance, "negative_inner_product")
})

test_that("data closest to itself", {
    n <- 1000
    d <- 10
    X <- matrix(rnorm(n * d), n, d)
    L <- LshTable(X)

    expect_equal(similar(L, as.vector(X[1,])), 1)
    expect_equal(similar(L, as.vector(X[2,])), 2)
    expect_equal(similar(L, as.vector(X[3,])), 3)
})


# example test
test_that("multiplication works", {
				  expect_equal(2 * 2, 4)
})

# load data for testing
data(L)

test_that("method='foo' fails", {
				  expect_error(conquerDE(L, method="foo"))
})

test_that("method='edgeRQLF' works", {
				  expect_s3_class(conquerDE(L, method="edgeRQLF"), "data.frame")
})

test_that("method='DESeq2' works", {
				  expect_s3_class(conquerDE(L, method="DESeq2"), "data.frame")
})

test_that("method='glmGamPoi' works", {
				  expect_s3_class(conquerDE(L, method="glmGamPoi"), "data.frame")
})

test_that("method='limma' works", {
				  expect_s3_class(conquerDE(L, method="limma"), "data.frame")
})

test_that("method='limmatrend' works", {
				  expect_s3_class(conquerDE(L, method="limmatrend"), "data.frame")
})

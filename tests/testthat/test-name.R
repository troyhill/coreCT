context("test-name")
library(coreCT)

test_that("voxDims gives correct output", {
  expect_equal(voxDims(directory = "core_426")$thickness.mm, 0.625)
  expect_error(voxDims(directory = "RDATA"))
  expect_error(voxDims(directory = "RDATA.dcm"))
  expect_error(voxDims(directory = "RDATA.dcm")) 
  expect_error(voxDims(directory = names(core_426))) # tough error message to trigger!
})

test_that("coreHist gives correct output", {
  expect_equal(length(coreHist(directory = "core_426", units = "percent")), 3)
  expect_equal(length(coreHist(directory = "core_426", units = "monkeys")), 3)
  expect_error(coreHist(directory = "RDATA"))
  expect_error(coreHist(directory = "RDATA.dcm"))
  expect_error(coreHist(directory = names(core_426))) # tough error message to trigger!
})

test_that("rootSize test", {
  #expect_equal(sum(rootSize(core_426$img, pixelA = 0.0596, waterHU = 663.912)[, 2]), 43)
  expect_equal(nrow(rootSize(mat.list = core_426$img, pixelA = 0.0596)), 3)
})


test_that("rootSizeDir test", {
  expect_equal(nrow(rootSizeDir("core_426", diameter.classes = c(2.5, 10))), 3)
  expect_error(rootSizeDir(directory = "RDATA"))
  expect_error(rootSizeDir(directory = "RDATA.dcm"))
  expect_error(rootSizeDir(directory = tempdir())) 
  expect_error(rootSizeDir(directory = names(core_426))) 
})



test_that("conv test", {
  expect_equal(nrow(conv(core_426$img, pixelA = 0.0596)), 3)
})

test_that("convDir test", {
  expect_equal(nrow(convDir("core_426", diameter.classes = c(2.5, 10))), 3)
  expect_equal(nrow(convDir("core_426", diameter.classes = c(2.5, 10), rootData = FALSE)), 3)
  expect_error(convDir(directory = "RDATA"))
  expect_error(convDir(directory = "RDATA.dcm"))
  expect_error(convDir(directory = tempdir())) 
  expect_error(convDir(directory = names(core_426))) 
})



test_that("getSurface test", {
  expect_equal(nrow(getSurface(a <- conv(core_426$img, pixelA = 0.0596), material = "peat", threshold = 1)), 3)
  expect_equal(nrow(getSurface(a <- conv(core_426$img, pixelA = 0.0596), material = "peat", threshold = 1, start = "both")), 3)
 # expect_equal(nrow(getSurface(a <- conv(core_426$img, pixelA = 0.0596), material = "rock_shell", threshold = 0.0011, start = "bottom")), 2)
})


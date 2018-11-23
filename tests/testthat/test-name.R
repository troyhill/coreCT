context("test-name")
library(coreCT)

test_that("voxDims gives correct output", {
  expect_equal(voxDims(directory = "core_426")$thickness.mm, 0.625)
  expect_error(voxDims(directory = "RDATA"))
  expect_error(voxDims(directory = "RDATA.dcm"))
  expect_error(voxDims(directory = "RDATA.dcm")) 
  expect_error(voxDims(directory = names(core_426))) # tough error message to trigger!
})

# test_that("fld.frq gives correct output", {
#   expect_equal(fld.frq(2, 1:10, units = "percent"), 0.8)
#   expect_equal(fld.frq(2, 1:10, units = "tides"), 9)
#   expect_error(fld.frq("string", 1:10)) 
#   expect_error(fld.frq(2, "string"))
#   expect_error(fld.frq(2, 1:10, units = "furlongs")) 
# })

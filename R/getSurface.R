#' @title Remove artificial surface layers from processed CT data
#'
#' @description Identifies and removes artificial surface layers from processed CT data
#'
#' @details Identifies and removes artificial surface layers from processed CT data. Areas can be removed from one or both ends of the core (set by \code{start}), based on exceeding a \code{threshold} proportion of material (e.g., 75% particles, sand, etc.)
#' 
#' @usage getSurface(x, material = "particles", threshold = 0.40, start = "top", thickness = 0.625)
#' 
#' @param x dataframe created by \code{conv}
#' @param material material used for determining where the surface begins
#' @param threshold decimal fraction of total area, used to determine the surface layer. Surface slices where \code{material} exceeds threshold value are removed.
#' @param start should core be processed from the top, bottom, or both?
#' @param thickness CT image thickness (mm)
#' 
#' @return value \code{getSurface} shortens the output of \code{conv} to remove artificial surface layers. The output is thus a subset of the input, and identical in structure to the /code{conv} output. 
#' 
#' @seealso \code{\link{conv}}
#' 
#' @examples
#' ### Not run:
#' \dontrun{data(core_426)
#' ct.slope <- unique(extractHeader(core_426$hdr, "RescaleSlope"))
#' ct.int   <- unique(extractHeader(core_426$hdr, "RescaleIntercept")) 
#' # convert raw units to Hounsfield units
#' HU_426 <- lapply(core_426$img, function(x) x*ct.slope + ct.int)
#' 
#' materials <- conv(HU_426)
#' head(materials[, 1:6], 20)
#' 
#' materials2 <- getSurface(materials)
#' head(materials2[, 1:6])
#' }
#' @export


getSurface <- function (x, material = "particles", threshold = 0.40, start = "top", thickness = 0.625) {
  # function identifies/removes areas at one end of the core (set by "start") exceeding a
  # threshold proportion of material (e.g., 75% particles, sand, etc.)
  x2 <- x
  if (start %in% c("top", "both")) {
    temp <- (x[, paste0(material, ".cm2")] / x$tot.cm2) > threshold
    if (temp[1]) { # if first value exceeds threshold
      rle.temp <- rle(temp)
      x2 <- x[-c(1:rle.temp$lengths[1]), ]
      x2$depth <- 1:nrow(x2) * thickness / 10
      rownames(x2) <-  1:nrow(x2)
    } else x2 <- x
  }
  if (start %in% c("bottom", "both")) {
    x2 <- x2[-order(x2[, "depth"]), ]
    temp <- (x2[, paste0(material, ".cm2")] / x$tot.cm2) > threshold
    if (temp[1]) { # if first value exceeds threshold
      rle.temp <- rle(temp)
      x2 <- x[-c(1:rle.temp$lengths[1]), ]
      x2$depth <- 1:nrow(x2) * thickness / 10
      rownames(x2) <-  1:nrow(x2)
    } else x2 <- x
  }
  x2
}
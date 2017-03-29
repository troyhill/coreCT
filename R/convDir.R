#' @title Convert a directory of raw DICOM images to material classes
#'
#' @description Calculates the area and volume of material classes for each CT slice in a directory
#'
#' @details Calculates the area and volume of material classes for each CT slice in a directory. Unlike \code{\link{conv}}, \code{\link{convDir}} accepts a folder of raw values and makes the conversion to Hounsfield Units using the metadata associated with the DICOM images.
#' 
#' @usage convDir(directory, upperLim = 3045, lowerLim = -1024, 
#' airHU = -850.3233, airSD = 77.6953, 
#' SiHU = 271.7827, SiSD = 39.2814,
#' glassHU = 1345.0696, glassSD = 45.4129,
#' waterHU = 63.912, waterSD = 14.1728,
#' densities = c(0.0012, 1, 1.23, 2.2))
#' 
#' @param directory folder of DICOM images(raw values)
#' @param upperLim upper bound cutoff for pixels (Hounsfield Units)
#' @param lowerLim lower bound cutoff for pixels (Hounsfield Units)
#' @param airHU mean value for air-filled calibration rod (Hounsfield Units)
#' @param airSD standard deviation for air-filled calibration rod
#' @param SiHU mean value for colloidal silica calibration rod 
#' @param SiSD standard deviation for colloidal Si calibration rod
#' @param glassHU mean value for glass calibration rod
#' @param glassSD standard deviation for glass calibration rod
#' @param waterHU mean value for water filled calibration rod
#' @param waterSD standard deviation for water filled calibration rod
#' @param densities numeric vector of known cal rod densities. Format must be c(air, water, Si, glass)
#' 
#' @return value \code{convDir} returns a dataframe with one row per CT slice. Values returned are the area and volume of 7 material classes: gas, peat, roots and rhizomes, rock and shell, fine mineral particles, sand, and water.
#' 
#' @seealso \code{\link{convDir}} is a wrapper for \code{\link{conv}}. \code{\link{rootSizeDir}} operates similarly.
#' 
#' @examples
#' \dontrun{ 
#' CTdir <- "C:/RDATA/CT/426/"
#' materials <- convDir(CTdir)
#' 
#' # plot using "ggplot" package after transforming with "reshape2" package
#' mass.long <- reshape2::melt(materials, id.vars = c("depth"), 
#'    measure.vars = grep(".g", names(materials)))
#' ggplot2::ggplot(data = mass.long, ggplot2::aes(y = -depth, x = value, 
#'    color = variable)) + ggplot2::geom_point() + ggplot2::theme_classic() + 
#'    ggplot2::xlab("mass per section (g)") 
#' }
#' 
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @importFrom stats aggregate
#' @importFrom stats coef
#' @importFrom stats lm
#' @importFrom oro.dicom extractHeader
#' @importFrom oro.dicom readDICOM
#' 
#' @export

# TODO: add command like 
# temp$bin <- cut(temp[(temp > lowerLim) & (temp < upperLim)], breaks = c(), labels = splits$material)
# and summarize by category
convDir <- function(directory, upperLim = 3045, lowerLim = -1024,
                 airHU = -850.3233, airSD = 77.6953, # all cal rod arguments are in Hounsfield Units
                 SiHU = 271.7827, SiSD = 39.2814,
                 glassHU = 1345.0696, glassSD = 45.4129,
                 waterHU = 63.912, waterSD = 14.1728,
                 densities = c(0.0012, 1, 1.23, 2.2) # format = air, water, Si, glass
) {
  # load DICOMs, takes a couple minutes
  fname   <- readDICOM(directory, verbose = TRUE) 
  # scrape some metadata
  pixelArea <- as.numeric(strsplit(fname$hdr[[1]]$value[fname$hdr[[1]]$name %in% "PixelSpacing"], " ")[[1]][1])^2
  ct.slope <- unique(extractHeader(fname$hdr, "RescaleSlope"))
  ct.int   <- unique(extractHeader(fname$hdr, "RescaleIntercept")) 
  # convert raw units to Hounsfield units
  HU <- lapply(fname$img, function(x) x*ct.slope + ct.int)
  thick <- unique(extractHeader(fname$hdr, "SliceThickness"))
  # pass data to conv()
  returnDat <- conv(mat.list = HU, pixelA = pixelArea, thickness = thick,
                         airHU, airSD,
                         SiHU, SiSD,
                         glassHU, glassSD,
                         waterHU, waterSD)
  returnDat
}
#' @title Convert a directory of raw DICOM images to material classes
#'
#' @description Calculates the area and volume of material classes for each CT slice in a directory
#'
#' @details Calculates the area and volume of material classes for each CT slice in a directory. Unlike \code{\link{conv}}, \code{\link{convDir}} accepts a folder of raw values and makes the conversion to Hounsfield Units using the metadata associated with the DICOM images.
#' 
#' @usage convDir(directory = file.choose(), upperLim = 3045, lowerLim = -1025, 
#' airHU = -850.3233, airSD = 77.6953, 
#' SiHU = 271.7827, SiSD = 39.2814,
#' glassHU = 1345.0696, glassSD = 45.4129,
#' waterHU = 63.912, waterSD = 14.1728,
#' densities = c(0.0012, 1, 1.23, 2.2),
#' rootData = TRUE, 
#' diameter.classes = c(1, 2, 2.5, 10), 
#' class.names = diameter.classes,
#' pixel.minimum = 4)
#' 
#' @param directory a character string that can be a matrix of DICOM images or the address of an individual DICOM file in a folder of DICOM images. The default action is <code>file.choose()</code>; a browser menu appears so the user can select the the desired directory by identifying a single DICOM file in the folder of images.
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
#' @param rootData if TRUE, \code{rootSize} is also called on the matrix
#' @param diameter.classes if rootData is TRUE, this argument provides an integer vector of diameter cut points used by \code{rootSize}. Units are mm (zero is added in automatically).
#' @param class.names placeholder, not used presently
#' @param pixel.minimum minimum number of pixels needed for a clump to be identified as a root
#' 
#' @return value \code{convDir} returns a dataframe with one row per CT slice. Values returned are the area and volume of seven material classes: gas, peat, roots and rhizomes, rock and shell, fine mineral particles, sand, and water. If \code{rootData = TRUE}, the output will also contain data on the abundance (number of particles), volume (cm3), and external surface area (cm2) of the root size classes specified in the \code{diameter.classes} argument.
#' 
#' @seealso \code{\link{convDir}} is a wrapper for \code{\link{conv}}. \code{\link{rootSizeDir}} operates similarly.
#' 
#' @examples
#' materials <- convDir("core_426", rootData = FALSE)
#' 
#' 
#' \dontrun{ 
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
#' @importFrom plyr join_all
#' @importFrom oro.dicom extractHeader
#' @importFrom oro.dicom readDICOM
#' 
#' @export

# TODO: add command like 
# temp$bin <- cut(temp[(temp > lowerLim) & (temp < upperLim)], breaks = c(), labels = splits$material)
# and summarize by category
convDir <- function(directory = file.choose(), 
                 upperLim = 3045, lowerLim = -1025,
                 airHU = -850.3233, airSD = 77.6953, # all cal rod arguments are in Hounsfield Units
                 SiHU = 271.7827, SiSD = 39.2814,
                 glassHU = 1345.0696, glassSD = 45.4129,
                 waterHU = 63.912, waterSD = 14.1728,
                 densities = c(0.0012, 1, 1.23, 2.2), # format = air, water, Si, glass
                 rootData = TRUE, diameter.classes = c(1, 2, 2.5, 10), 
                 class.names = diameter.classes,
                 pixel.minimum = 4
) {
  if (!exists(directory)) { # is "directory" an existing object (user-loaded DICOM matrix)
    
    if (substr(directory, nchar(directory) - 3, nchar(directory)) %in% ".dcm") {
      directory <- dirname(directory)
    } else if (grep("/", directory) == 1) { # dangerously assumes that if there's a forward slash, it's a valid address
      # directory <- directory # do nothing
    } else stop("Incorrect directory name: directory specified in a character string must end with a '/'; if 'file.choose()' is used, the selected file must be a dicom image")
    # load DICOMs, takes a couple minutes
    fname   <- readDICOM(directory, verbose = TRUE) 
    
  } else if (exists(directory) & (sum(names(get(directory)) %in% c("hdr", "img")) == 2)){ # could have better error checking here
    fname <- get(directory)
  } else stop("Invalid input: 'directory' object or file location is incorrectly specified.")
  # scrape some metadata
  # pixelArea <- voxDims(directory)$pixelArea.mm2
  # thick     <- voxDims(directory)$thickness.mm
  pixelArea <- as.numeric(strsplit(fname$hdr[[1]]$value[fname$hdr[[1]]$name %in% "PixelSpacing"], " ")[[1]][1])^2
  thick <- unique(extractHeader(fname$hdr, "SliceThickness"))
  
  # convert raw units to Hounsfield units
  ct.slope <- unique(extractHeader(fname$hdr, "RescaleSlope"))
  ct.int   <- unique(extractHeader(fname$hdr, "RescaleIntercept")) 
  HU <- lapply(fname$img, function(x) x*ct.slope + ct.int)
  
  # pass data to conv()
  returnDat <- conv(mat.list = HU, pixelA = pixelArea, thickness = thick,
                         upperLim = upperLim, lowerLim = lowerLim,
                         airHU = airHU, airSD = airSD, 
                         waterHU = waterHU, waterSD = waterSD,
                         SiHU = SiHU, SiSD = SiSD,
                         glassHU = glassHU, glassSD = glassSD,
                         densities = densities)
  if (rootData == TRUE) {
    rootsDat <- rootSize(mat.list = HU, pixelA = pixelArea, thickness = thick, 
                         diameter.classes = diameter.classes, class.names = diameter.classes, 
                         airHU = airHU, airSD = airSD, 
                         waterHU = waterHU, waterSD = waterSD,
                         pixel.minimum = pixel.minimum)
    # returnDat <- cbind(returnDat, rootsDat)
    returnDat <- plyr::join_all(list(returnDat, rootsDat), by = "depth")
  }
  
  returnDat
}
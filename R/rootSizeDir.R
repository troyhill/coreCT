#' @title Convert a directory of raw DICOM images to  root particle counts and surface areas
#'
#' @description Calculates the number of root/rhizome particles and surface areas, for different size classes
#'
#' @details Calculates the number of root/rhizome particles and surface areas, for different size classes. Unlike \code{\link{rootSize}}, \code{\link{rootSizeDir}} accepts a folder of raw values and makes the conversion to Hounsfield Units using the metadata associated with the DICOM images.
#' 
#' @usage rootSizeDir(directory = file.choose(), diameter.classes = c(1, 2, 5, 10, 20),
#' class.names = diameter.classes,
#' airHU = -850.3233,
#' airSD = 77.6953,
#' waterHU = 63.912,
#' waterSD = 14.1728,
#' pixel.minimum = 1)
#' 
#' @param directory a character string that can be a matrix of DICOM images or the address of an individual DICOM file in a folder of DICOM images. The default action is <code>file.choose()</code>; a browser menu appears so the user can select the the desired directory by identifying a single DICOM file in the folder of images.
#' @param diameter.classes an integer vector of diameter cut points. Units are mm (zero is added in automatically).
#' @param class.names not used presently
#' @param airHU mean value for air-filled calibration rod (all rod arguments are in Hounsfield Units)
#' @param airSD standard deviation for air-filled calibration rod
#' @param waterHU mean value for water-filled calibration rod
#' @param waterSD standard deviation for water-filled calibration rod
#' @param pixel.minimum minimum number of pixels needed for a clump to be identified as a root
#' 
#' @return value \code{rootSize} returns a dataframe with one row per CT slice. Values returned are the number, volume (cm3), and surface area (cm2) of particles in each size class with an upper bound defined in \code{diameter.classes}.
#' 
#' @seealso \code{\link{rootSizeDir}} is a wrapper for \code{\link{rootSize}}. \code{\link{rootSizeDir}} operates similarly.
#' 
#' @examples
#' rootChars <- rootSizeDir("core_426", diameter.classes = c(2.5, 10))
#' 
#' \dontrun{
#' # plot using "ggplot" package after transforming with "reshape2" package
#' area.long <- reshape2::melt(rootChars, id.vars = c("depth"), 
#'    measure.vars = grep("Area", names(rootChars)))
#' ggplot2::ggplot(data = area.long, ggplot2::aes(y = -depth, x = value, 
#'    color = variable)) + ggplot2::geom_point() + ggplot2::theme_classic() + 
#'    ggplot2::xlab("root external surface area per slice (cm2)")
#' }
#' 
#' @importFrom raster freq
#' @importFrom raster match
#' @importFrom raster boundaries
#' @importFrom raster clump
#' @importFrom igraph union
#' @importFrom igraph decompose
#' @importFrom igraph spectrum
#' @importFrom stats aggregate
#' @importFrom stats coef
#' @importFrom stats lm
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @importFrom oro.dicom extractHeader
#' @importFrom oro.dicom readDICOM
#' 
#' @export

rootSizeDir <- function (directory = file.choose(), 
                        diameter.classes = c(1, 2, 5, 10, 20), # mm, targets clumps less than or equal to "diameter" argument
                        class.names = diameter.classes,
                        airHU = -850.3233, # Hounsfield Units
                        airSD = 77.6953,
                        waterHU = 63.912,
                        waterSD = 14.1728,
                        pixel.minimum = 1) {
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

  # pass data to rootSize()
  returnDat <- rootSize(mat.list = HU, pixelA = pixelArea, thickness = thick, 
                        diameter.classes, class.names, airHU, airSD, waterHU, waterSD,
                        pixel.minimum = pixel.minimum)
  
  returnDat
}

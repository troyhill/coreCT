#' @title Extract voxel dimensions from DICOM image
#'
#' @description Extract pixel area and slice thickness from DICOM header to characterize voxel (3D pixel) dimensions.
#'
#' @usage voxDims(directory = file.choose())
#' 
#' @param directory a character string that can be a matrix of DICOM images or the address of an individual DICOM file in a folder of DICOM images. The default action is <code>file.choose()</code>; a browser menu appears so the user can select the the desired directory by identifying a single DICOM file in the folder of images.
#' 
#' @return value \code{voxDims} returns a two-column dataframe showing the pixel area and slice thickness. Values in the DICOM headers are assumed to be millimeters; pixel area and slice thickness columns are labeled based on this assumption.
#' 
#' @examples
#' # data(core_426)
#' voxDims("core_426")
#' 
#' @importFrom oro.dicom extractHeader
#' @importFrom oro.dicom readDICOM
#' 
#' @export

voxDims <- function(directory = file.choose()) {
  if (!exists(directory)) { # is "directory" an existing object (user-loaded DICOM matrix)
    
    if (substr(directory, nchar(directory) - 3, nchar(directory)) %in% ".dcm") {
      #directory <- dirname(directory)
    } else if (grep("/", directory) == 1) { # dangerously assumes that if there's a forward slash, it's a valid address
      # directory <- directory
    } else stop("Incorrect directory name: directory specified in a character string must end with a '/'; if 'file.choose()' is used, the selected file must be a dicom image")
    # load DICOMs, takes a couple minutes
    fname   <- oro.dicom::readDICOM(directory, verbose = TRUE) 
    
  } else if (exists(directory) & (sum(names(get(directory)) %in% c("hdr", "img")) == 2)){ # could have better error checking here
    fname <- get(directory)
  } else stop("Invalid input: 'directory' object or file location is incorrectly specified.")
  
  # scrape some metadata
  pixelArea <- as.numeric(strsplit(fname$hdr[[1]]$value[fname$hdr[[1]]$name %in% "PixelSpacing"], " ")[[1]][1])^2
  thick     <- unique(oro.dicom::extractHeader(fname$hdr, "SliceThickness"))
  
  returnDat <- data.frame(pixelArea.mm2 = pixelArea, thickness.mm = thick)
  returnDat
}
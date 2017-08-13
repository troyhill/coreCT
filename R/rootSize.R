#' @title Convert a matrix of semi-processed DICOM images to root particle counts, volumes, and surface areas
#'
#' @description Calculates the number of root/rhizome particles, volumes,  and surface areas, for different size classes
#'
#' @details Calculates the number of root/rhizome particles, volumes, and surface areas, for different size classes. This function requires that values be Hounsfield Units (i.e., data must be semi-processed from the raw DICOM imagery).
#' 
#' @usage rootSize(mat.list, pixelA, diameter.classes = c(1, 2, 2.5, 10),
#' class.names = diameter.classes,
#' thickness = 0.625,
#' airHU = -850.3233,
#' airSD = 77.6953,
#' waterHU = 63.912,
#' waterSD = 14.1728,
#' pixel.minimum = 4)
#' 
#' @param mat.list list of DICOM images for a sediment core (values in Hounsfield Units)
#' @param pixelA pixel area (mm2)
#' @param diameter.classes an integer vector of diameter cut points. Units are mm (zero is added in automatically).
#' @param class.names not used presently
#' @param thickness CT image thickness (mm)
#' @param airHU mean value for air-filled calibration rod (all rod arguments are in Hounsfield Units)
#' @param airSD standard deviation for air-filled calibration rod
#' @param waterHU mean value for water-filled calibration rod
#' @param waterSD standard deviation for water-filled calibration rod
#' @param pixel.minimum minimum number of pixels needed for a clump to be identified as a root
#' 
#' @return value \code{rootSize} returns a dataframe with one row per CT slice. Values returned are the number, volume (cm3), and surface area (cm2) of particles in each size class with an upper bound defined in \code{diameter.classes}.
#' 
#' @seealso \code{\link{conv}}
#' 
#' @examples
#' ct.slope <- unique(extractHeader(core_426$hdr, "RescaleSlope"))
#' ct.int   <- unique(extractHeader(core_426$hdr, "RescaleIntercept")) 
#' # convert raw units to Hounsfield units
#' HU_426 <- lapply(core_426$img, function(x) x*ct.slope + ct.int)
#' 
#' rootChars <- rootSize(HU_426, pixelA = 0.0596,
#' diameter.classes = c(2.5, 10))
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
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @importFrom oro.dicom extractHeader
#' @importFrom oro.dicom readDICOM
#' 
#' @export

rootSize <- function (mat.list, pixelA, 
                      diameter.classes = c(1, 2, 2.5, 10), # mm, targets clumps less than or equal to "diameter" argument
                      class.names = diameter.classes,
                      thickness = 0.625, # mm
                      airHU = -850.3233, # Hounsfield Units
                      airSD = 77.6953,
                      waterHU = 63.912,
                      waterSD = 14.1728,
                      pixel.minimum = 4) {
  # function creates dataframe with root/rhizome numbers and perimeter(!) for each depth interval
  pb <- utils::txtProgressBar(min = 0, max = length(mat.list), initial = 0, style = 3)
  voxelVol <- pixelA * thickness / 1e3 # cm3
  air.UB <- airHU + airSD # lower bound on R&R, exclusive. earl adds 1 to this quantity to start next category
  water.LB <- waterHU - waterSD # upper bound on R&R, inclusive
  # clump IDs extracted on basis of area
  # with rule: clumps must be greater than 1 pixel
  # this rule is based on comparison of ImageJ and R results
  diams <- c(0, diameter.classes) # revise diameter classes to include zero
  thresh.A <- pi*(diams/2)^2 # convert diameter to area
  thresh.pixels <- round(thresh.A / pixelA)  # convert area to no of pixels (irregular clumps are included; should maybe use floor() instead of round())
  thresh.pixels[1] <- pixel.minimum
  
  for (i in 1:length(mat.list)) {
    depth <- thickness * i / 10 # cm
    temp <- mat.list[[i]]
    
    ### adjust here to get size class mass estimates
    temp[(temp > air.UB) & (temp <= water.LB)] <- 1 # isolate pixels with R&R density, change to 1
    temp[!temp == 1] <- NA # make sure only two values exist: 1, NA
    if(length(temp) == 0) next 
    
    if (length(!is.na(temp)) > 0) {
      rmat <- raster::raster(temp)
      rmat.int <- raster::clump(rmat, directions = 8, gaps = FALSE)
      clump.sub1 <- data.frame(freq(rmat.int))
      
      for (j in 2:length(diams)) {
        # select clumps with less than x pixels
        clump.sub <- clump.sub1[ (clump.sub1$count <= thresh.pixels[j]) & 
                                   (clump.sub1$count > thresh.pixels[j-1]) & 
                                   !is.na(clump.sub1$value), ] #clump.sub$A <- clump.sub$count * pixelArea
        includeID <- as.vector(clump.sub$value) # record IDs from clumps which met the criteria in previous step
        rootVol  <- sum(as.vector(clump.sub$count) * voxelVol) # sum clump volumes in cm3
        # get root perimeter (then multiply by thickness to calculate external surface area)
        # make a new raster to be sieved
        maskSieve <- rmat.int
        # assign NA to all clumps whose IDs are NOT found in excludeID
        # maskSieve[!maskSieve %in% includeID] <- NA # 
        maskSieve <- raster::match(maskSieve, includeID)
        
        # calculate perimeter based on clump results # plot(boundaries(maskSieve, classes = FALSE, directions = 8, asNA = TRUE))
        a2 <- raster::freq(raster::boundaries(maskSieve, classes = FALSE, directions = 8, asNA = TRUE)) # number of "1"s x pixelSide = edge length
        if (nrow(clump.sub) > 0) {
          a3 <- a2[[3]] * sqrt(pixelA) # mm of edge length; sqrt() reflects assumption that one side of pixel contributes to perimeter (neglects corners; lower-bound estimate) # (sqrt(2*sqrt(pixelArea)^2) - sqrt(pixelArea)) / sqrt(pixelArea)
        } else if (nrow(clump.sub) == 0) {
          a3 <- 0
        }
        numberOfClumps <- nrow(clump.sub)
        rootSurfaceArea <- a3 * thickness / 100 # edge length (mm) * thickness (mm) = mm2 /100 = cm2 of external root surface area in slice
        # rootSurfaceVol <- rootSurfaceArea * (thickness / 10) # cm3 # tdh: probably not a meaningful parameter
        outDatInt <- data.frame(particles = numberOfClumps, surfArea = rootSurfaceArea, rootVolume = rootVol) #,surfaceVol = rootSurfaceVol)
        names(outDatInt) <- paste0(names(outDatInt), ".", diams[j - 1], "_", diams[j], "mm")
        
        if (j == 2) {
          outDatInt2 <- outDatInt 
        } else {
          outDatInt2 <- cbind(outDatInt2, outDatInt)
        }
      }
    } else if (length(!is.na(temp)) == 0) { # fill in zeroes if there aren't any clumps
      for (j in 2:length(diams)) {
        outDatInt <- data.frame(particles = 0, surfArea = 0) #,surfaceVol = 0)
        names(outDatInt) <- paste0(names(outDatInt2), ".", diams[j - 1], "_", diams[j], "mm")
        if (j == 2) {
          outDatInt2 <- outDatInt 
        } else {
          outDatInt2 <- cbind(outDatInt2, outDatInt)
        }
      }
    }
    # combine dfs from individual layers
    # convert NAs to zeroes
    outDatInt2$depth <- depth
    
    if ((i > 1)  & exists("outDat")) {
      outDat <- rbind(outDat, outDatInt2)
    } else {
      outDat <- outDatInt2
    }
    utils::setTxtProgressBar(pb, i)
  }
  outDat <- outDat[ , order(names(outDat))] # sort columns alphabetically (by material), rather than size class
  ### TODO: aggregate parameters
  outDat$structures <- base::rowSums(outDat[, grep("particles", names(outDat))])
  outDat$totArea    <- base::rowSums(outDat[, grep("Area", names(outDat))])
  # outDat$totVol <- base::rowSums(outDat[, grep("Vol", names(outDat))])
  outDat
}

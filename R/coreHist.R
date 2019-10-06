#' @title Whole-core frequency distribution of Hounsfield units
#'
#' @description Provides the raw data and plots a frequency distibution for Hounsfield Units in the entire core, also delineating material classes. As of coreCT version 1.3.0, this code accommodates calibration curves with >4 calibrants, and uses density thresholds converted to Hounsfield Units using the calibration curve (rather than direct calibration rod values) to partition sediment components. 
#'
#' @usage coreHist(directory = file.choose(), 
#' units = "percent",
#' upperLim = 3045, lowerLim = -1025,
#' means     = c(-850.3233, 63.912, 271.7827, 1345.0696),
#' sds       = c(77.6953, 14.1728, 39.2814, 45.4129),
#' densities = c(0.0012, 1, 1.23, 2.2),
#' returnData = TRUE, pngName = NULL)
#' 
#' @param directory a character string that can be (1) a matrix of DICOM images that exists in the global environment, or (2) the address of an individual DICOM file in a folder of DICOM images. The default action is <code>file.choose()</code>; a browser menu appears so the user can select the the desired directory by identifying a single DICOM file in the folder of images.
#' @param units units to be used for plotting purposes: either "percent" (the default) or "absolute"
#' @param upperLim upper bound cutoff for pixels (Hounsfield Units); upper bound is inclusive
#' @param lowerLim lower bound cutoff for pixels (Hounsfield Units); lower bound is exclusive
#' @param means mean values (units = Hounsfield Units) for calibration rods used.
#' @param sds standard deviations (units = Hounsfield Units) for calibration rods used. Must be in the same order as \code{means}.
#' @param densities numeric vector of known cal rod densities. Must be in the same order as \code{means} and \code{sds}.
#' @param returnData if \code{TRUE}, voxel counts for each Hounsfield unit from \code{lowerLim} to \code{upperLim} are returned, as are material class definitions. These are the data needed to re-create and modify the frequency plot.
#' @param pngName if this is not \code{NULL}, the frequency plot is saved to disk. In that case, \code{pngName} should be a character string containing the name and address of the file. 
#' 
#' @return list if \code{returnData = TRUE}, a list is returned containing (1) the frequencies for each Hounsfield unit value from \code{lowerLim} to \code{upperLim}, (2) the boundaries for material classes, and (3) a summary of the calibration curve applied. Lower boundaries for a component class are exclusive, while upper bounds are inclusive. These materials allow the frequency distribution to be plotted by the user. If \code{returnData = FALSE} the data are plotted in the graphics window, but nothing is preserved.
#' 
#' @examples
#' # data(core_426)
#' coreHist("core_426", returnData = FALSE)
#' 
#' @importFrom oro.dicom extractHeader
#' @importFrom oro.dicom readDICOM
#' @importFrom grDevices dev.off
#' @importFrom grDevices png
#' @importFrom graphics plot
#' @importFrom graphics abline
#' @importFrom graphics text
#' @importFrom graphics par
#' @importFrom stats predict
#' 
#' @export

coreHist <- function(directory = file.choose(), 
                     units = "percent",
                     upperLim = 3045, lowerLim = -1025,
                     means     = c(-850.3233, 63.912, 271.7827, 1345.0696),  # all cal rod units are in Hounsfield Units
                     sds       = c(77.6953, 14.1728, 39.2814, 45.4129),  # in same order as means. units are in Hounsfield Units
                     densities = c(0.0012, 1, 1.23, 2.2), # must be in the same order as the means and SDs. units are g/cm3
                     returnData = TRUE, pngName = NULL) {
  if (!exists(directory)) { # is "directory" an existing object (user-loaded DICOM matrix)
    
    if (substr(directory, nchar(directory) - 3, nchar(directory)) %in% ".dcm") {
      directory <- dirname(directory)
    } else if (grepl("/", directory) == 1) { # dangerously assumes that if there's a forward slash, it's a valid address
      # directory <- directory # do nothing
    } else stop("Incorrect directory name: directory specified in a character string must end with a '/'; if 'file.choose()' is used, the selected file must be a dicom image")
    # load DICOMs, takes a couple minutes
    fname   <- oro.dicom::readDICOM(directory, verbose = TRUE) 
    
  } else if (exists(directory) & (sum(names(get(directory)) %in% c("hdr", "img")) == 2)){ # could have better error checking here
    fname <- get(directory)
  } else stop("Invalid input: 'directory' object or file location is incorrectly specified.")
  
  ##### section added 20190922 to separate calibration curve from component partitioning
  densitydf <- data.frame(HU = means, density = densities)
  calCurve <- summary(lm1 <- stats::lm(density ~ HU, data = densitydf)) # density in g/cm3
  if (summary(lm1)$r.squared < 0.95) { # print message to console if r2 < 0.95 
    message(cat("\n Note: Calibration curve has limited explanatory power; r2 = ", round(summary(lm1)$r.squared, 3), "\n"))
  }
  densToHU <- stats::lm(HU ~ density, data = densitydf)
  
  partition.means <- stats::predict(densToHU, newdata = data.frame(density = c(0.0012, 1, 1.23, 2.2))) # this is hard-coded to reflect partitioning in Davey et al. 2011
  
  air.proxy   <- which.min(abs(means - partition.means[1])) # which element to use for air SD?
  water.proxy <- which.min(abs(means - partition.means[2])) # which element to use for water SD?
  Si.proxy    <- which.min(abs(means - partition.means[3])) # which element to use for Si SD?
  glass.proxy <- which.min(abs(means - partition.means[4])) # which element to use for glass SD?
  
  partition.sds   <- c(sds[c(air.proxy, water.proxy, Si.proxy, glass.proxy)])
  
  ### now, set means and SDs used for partitioning
  airHU      <- partition.means[1]
  airSD      <- partition.sds[1]
  waterHU    <- partition.means[2]
  waterSD    <- partition.sds[2]
  SiHU       <- partition.means[3]
  SiSD       <- partition.sds[3]
  glassHU    <- partition.means[4]
  glassSD    <- partition.sds[4]
  ##### end section added 20190922
  
  
  # divisions between material classes
  splits <- data.frame(material = c("air",             "RR",                "water",         "peat",            "particles",         "sand",                   "rock_shell"),
                       lower = c(round(lowerLim),      round(airHU + airSD), round(waterHU - waterSD), round(waterHU + waterSD),    round(SiHU + SiSD), 750,                      round(glassHU + glassSD)), 
                       #lower = c(round(lowerLim),        round(airHU+airSD) + 1, round(water.LB) + 1, round(water.UB) + 1,    round(SiHU + SiSD) + 1, 750 + 1,                      round(glassHU + glassSD) + 1), 
                       upper = c(round(airHU + airSD), round(waterHU - waterSD),      round(waterHU + waterSD), round(SiHU + SiSD), 750,                round(glassHU + glassSD), round(upperLim)))
  
  
  pixelArea <- as.numeric(strsplit(fname$hdr[[1]]$value[fname$hdr[[1]]$name %in% "PixelSpacing"], " ")[[1]][1])^2
  thick <- unique(oro.dicom::extractHeader(fname$hdr, "SliceThickness"))
  
  # convert raw units to Hounsfield units
  ct.slope <- unique(oro.dicom::extractHeader(fname$hdr, "RescaleSlope"))
  ct.int   <- unique(oro.dicom::extractHeader(fname$hdr, "RescaleIntercept")) 
  HU <- lapply(fname$img, function(x) x*ct.slope + ct.int)
  
  # tempDat <- as.data.frame(table(unlist(HU))) # Takes a long time. started at 7:27pm
  # tempDat <- lapply(HU, table)
  # tempDat2 <- plyr::join_all(lapply(tempDat, as.data.frame, stringsAsFactors = FALSE), by = "Var1")
  # tempDat2.temp <- do.call(rbind, HU)
  # tempDat2 <- data.frame(table(tempDat2.temp))
  # names(tempDat2)[1] <- "Var1"
  # tempDat2$Var1 <- as.numeric(tempDat2$Var1)
  # tempDat2$finalFreq <- rowSums(tempDat2[, 2:ncol(tempDat2)], na.rm = TRUE)
  
  for ( i in 1:length(HU)) {
    if ( i == 1) {
      test <- tabulate((HU[[i]] - lowerLim), nbins = upperLim - lowerLim + 1) # all values need to be positive integers
    } else {
      test <- test + tabulate((HU[[i]]  - lowerLim), nbins = upperLim - lowerLim + 1)
    }
  }
  tempDat2 <- data.frame(Var1 = (lowerLim + 1):(upperLim + 1), finalFreq = test)
  
  
  if (units == "percent") {
    tot <- sum(tempDat2$finalFreq, na.rm = TRUE)
    ylabel <- "Frequency (% of total voxels)"
  } else {
    tot <- 1
    ylabel <- "Frequency (no. of voxels)"
  }
  if (!is.null(pngName)) {
    grDevices::png(file = pngName, width = 4, height = 3.5, units = "in")
  }
  graphics::par(mar = c(4, 4, 0.5, 0.5))
  graphics::plot(finalFreq / tot ~ Var1, tempDat2[(tempDat2$Var1 < upperLim) & (tempDat2$Var1 > lowerLim), ], cex = 0.7, 
       pch = 19, las = 1, ylab = ylabel, xlab = "HU", xlim = c(lowerLim, upperLim))
  # add lines and label material classes
  graphics::abline(v = c(lowerLim + 1, splits$upper))
  verticalTextPosition <- max(tempDat2[(tempDat2$Var1 < upperLim) & (tempDat2$Var1 > lowerLim), ], na.rm = TRUE) * c(0.90, 0.70) / tot
  graphics::text(x = (splits$upper + splits$lower)/2, y = verticalTextPosition, labels = as.character(splits$material)) # all classes
  # graphics::text(x = (splits$upper[-3] + splits$lower[-3])/2, y = verticalTextPosition, labels = as.character(splits$material)[-3]) # without water
  
  if (!is.null(pngName)) {
    grDevices::dev.off()
  }
  
  if (returnData == TRUE) {
    outDat <- tempDat2[(tempDat2$Var1 < upperLim) & (tempDat2$Var1 > lowerLim), c("Var1", "finalFreq")]
    return(list(histData = outDat, splits = splits, calCurve = calCurve))
  }
  
}

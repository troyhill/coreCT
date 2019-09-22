#' @title Convert a matrix of semi-processed DICOM images to mass and volume of material classes
#'
#' @description Converts raw CT units to material classes for each CT slice.
#'
#' @details Calculates average Hounsfield units, cross-sectional areas (cm2), volumes (cm3), and masses (g) of material classes for each CT slice. This function assumes that core walls and all non-sediment material have been removed from the raw DICOM imagery. This function converts data from raw x-ray attenuation values to Hounsfield Units, and then uses user-defined calibration rod inputs to categorize sediment components: air, roots and rhizomes, peat, water, particulates, sand, and rock/shell. The input style for calibration rods ensures sediment components are partitioned following the density divisions in Davey et al. 2011. Calibration rods and are used to develop the calibration curve. Separately, the densities used for partitioning in Davey et al. 2011 (0.0012, 1, 1.23, 2.2 g/cm3) are converted to Hounsfield Units and used for partitioning sediment components. The standard deviation for the calibration rod nearest to the target value is used for the standard deviation for the division between two sediment components.
#' 
#' @usage conv(mat.list, upperLim = 3045, lowerLim = -1025, 
#' pixelA, thickness = 0.625, # all in mm 
#' means     = c(-850.3233, 63.912, 271.7827, 1345.0696),
#' sds       = c(77.6953, 14.1728, 39.2814, 45.4129),
#' densities = c(0.0012, 1, 1.23, 2.2))
#' 
#' @param mat.list list of DICOM images for a sediment core (values in Hounsfield Units)
#' @param upperLim upper bound cutoff for pixels (Hounsfield Units)
#' @param lowerLim lower bound cutoff for pixels (Hounsfield Units)
#' @param pixelA pixel area (mm2)
#' @param thickness slice thickness for computed tomography image series (mm)
#' @param means mean values (units = Hounsfield Units) for calibration rods used.
#' @param sds standard deviations (units = Hounsfield Units) for calibration rods used. Must be in the same order as \code{means}.
#' @param densities numeric vector of known cal rod densities. Must be in the same order as \code{means} and \code{sds}.
#' 
#' @return value \code{conv} returns a dataframe with one row per CT slice. Values returned are the average Hounsfield Unit value, the area (cm2), volume (cm3), and mass (grams) of 7 material classes: gas, peat, roots and rhizomes, particulates, sand, water, and rock/shell. If <code>rootData = TRUE</code>, data for specified root size classes are also returned. See <code>rootSize</code> for more detail on those values.
#' 
#' @seealso \code{\link{rootSize}} operates similarly.
#' 
#' @examples
#' ct.slope <- unique(extractHeader(core_426$hdr, "RescaleSlope"))
#' ct.int   <- unique(extractHeader(core_426$hdr, "RescaleIntercept")) 
#' # convert raw units to Hounsfield units
#' HU_426 <- lapply(core_426$img, function(x) x*ct.slope + ct.int)
#' 
#' materials <- conv(HU_426, pixelA = 0.0596)
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
#' @importFrom stats predict
#' @importFrom oro.dicom extractHeader
#' @importFrom oro.dicom readDICOM
#' 
#' @export

conv <- function(mat.list, upperLim = 3045, lowerLim = -1025,
                  pixelA, thickness = 0.625, # all in mm
                  means     = c(-850.3233, 63.912, 271.7827, 1345.0696),  # all cal rod units are in Hounsfield Units
                  sds       = c(77.6953, 14.1728, 39.2814, 45.4129),  # in same order as means. units are in Hounsfield Units
                  densities = c(0.0012, 1, 1.23, 2.2) # must be in the same order as the means and SDs. units are g/cm3
) {
  voxelVol <- pixelA * thickness / 1e3 # cm3
  
  
  
  ##### section added 20190922 to separate calibration curve from component partitioning
  densitydf <- data.frame(HU = means, density = densities)
  summary(lm1 <- stats::lm(density ~ HU, data = densitydf)) # density in g/cm3
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
  
  
  
  ### now calculate borders between sediment components based on density targets; irrespective of calibrants
  water.LB <- waterHU - waterSD
  water.UB <- waterHU + waterSD
  # note: Earl adds 1 to switch between categories
  splits <- data.frame(material = c("air",             "RR",                "water",         "peat",            "particulates",         "sand",                   "rock_shell"),
                       lower = c(round(lowerLim),      round(airHU + airSD), round(water.LB), round(water.UB),    round(SiHU + SiSD), 750,                      round(glassHU + glassSD)), 
                       #lower = c(round(lowerLim),        round(airHU+airSD) + 1, round(water.LB) + 1, round(water.UB) + 1,    round(SiHU + SiSD) + 1, 750 + 1,                      round(glassHU + glassSD) + 1), 
                       upper = c(round(airHU + airSD), round(water.LB),      round(water.UB), round(SiHU + SiSD), 750,                round(glassHU + glassSD), round(upperLim)))
  
  # start replace -----------------------------------------------------------
  areaDat <- list()
  volDat  <- list()
  HU_Dat  <- list()
  massDat <- list()
  
  pb <- utils::txtProgressBar(min = 0, max = length(mat.list), initial = 0, style = 3)
  
  for ( i in 1:length(mat.list)) {
    test <- tabulate((mat.list[[i]] - lowerLim), nbins = upperLim - lowerLim + 1) #counts from mat.list = lowerLim:mat.list = upperLim
    component.borders <- diff(unique(c(splits$lower, upperLim)))
    area_output.int <- diff(c(0,
                              sum(test[1:component.borders[1]]),
                              sum(test[1:sum(component.borders[1:2])]),
                              sum(test[1:sum(component.borders[1:3])]),
                              sum(test[1:sum(component.borders[1:4])]),
                              sum(test[1:sum(component.borders[1:5])]),
                              sum(test[1:sum(component.borders[1:6])]),
                              sum(test[1:sum(component.borders[1:7])])
    ) * pixelA / 1e2) # area of each category - cm2
    
    areaDat[[i]]    <- area_output.int
    volDat[[i]]     <- area_output.int * (thickness / 10) # cm3
    HU_Dat[[i]]     <- tapply(mat.list[[i]], cut(mat.list[[i]], breaks = c(splits$lower, upperLim), right = TRUE, include.lowest = FALSE, labels = splits$material), mean)
    
    # wetMass     <- (mat.list[[i]] * stats::coef(lm1)[2] + stats::coef(lm1)[1]) * voxelVol  # convert to g/cm3 and then to g (wet) in each pixel
    # massDat[[i]]    <- tapply(wetMass, cut(wetMass, 
    #                       breaks = c(c(splits$lower, upperLim) * stats::coef(lm1)[2] + stats::coef(lm1)[1]) * voxelVol,
    #                       right = TRUE, include.lowest = FALSE, labels = splits$material),
    #                       sum)
    # could skip this by calculating as ((mean_HU * slope + coef) * voxelVol) * area_output.int (i.e., mean HU converted to density * voxel count)
    massDat[[i]]    <- (HU_Dat[[i]] * stats::coef(lm1)[2] + stats::coef(lm1)[1]) * volDat[[i]] # * voxelVol * (area_output.int * thickness / 10)
    utils::setTxtProgressBar(pb, i)
  }
  
  HU_final     <- do.call(rbind, HU_Dat)
  area_final   <- do.call(rbind, areaDat)
  vols_final   <- do.call(rbind, volDat)
  masses_final <- do.call(rbind, massDat)
  
  names(area_final)   <- paste0(names(HU_final), ".cm2")
  names(vols_final)   <- paste0(names(HU_final), ".cm3")
  names(masses_final) <- paste0(names(HU_final), ".g")
  names(HU_Dat)     <- paste0(names(HU_final), ".HU")
  outDat <- as.data.frame(do.call(cbind, list(1:length(mat.list) * thickness / 10, HU_final, area_final, vols_final, masses_final) ))
  names(outDat) <- c("depth", paste0(splits$material, c(".HU")), paste0(splits$material, c(".cm2")),
                     paste0(splits$material, c(".cm3")),
                     paste0(splits$material, c(".g")))
  outDat$tot.cm2    <- base::rowSums(outDat[, grep(".cm2", names(outDat))], na.rm = TRUE)
  outDat$tot.cm3    <- base::rowSums(outDat[, grep(".cm3", names(outDat))], na.rm = TRUE)
  outDat$tot.g      <- base::rowSums(outDat[, grep(".g", names(outDat))], na.rm = TRUE)
  # outDat$tot.meanHU <- base::rowSums(outDat[, grep(".HU", names(outDat))]) * 
  
  return(outDat)
} 
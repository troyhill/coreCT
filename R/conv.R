#' @title Convert Hounsfield Units to material classes
#'
#' @description Calculates the area and volume of material classes for each CT slice
#'
#' @details Calculates the area and volume of material classes for each CT slice
#' 
#' @usage conv(mat.list, upperLim = 3045, lowerLim = -1024, 
#' pixelA, thickness = 0.625, # all in mm 
#' airHU = -850.3233, 
#' SiHU = 271.7827, glassHU = 1345.0696,
#' waterHU = 63.912, waterSD = 14.1728,
#' densities = c(0.0012, 1, 1.23, 2.2))
#' 
#' @param mat.list list of DICOM images for a sediment core
#' @param upperLim upper bound cutoff for pixels (Hounsfield Units)
#' @param lowerLim lower bound cutoff for pixels (Hounsfield Units)
#' @param pixelA pixel area (mm2)
#' @param thickness CT image thickness (mm)
#' @param airHU mean value for air-filled calibration rod (Hounsfield Units)
#' @param SiHU mean value for colloidal silica calibration rod (Hounsfield Units)
#' @param glassHU mean value for glass calibration rod (Hounsfield Units)
#' @param waterHU mean value for water-filled calibration rod (Hounsfield Units)
#' @param waterSD standard deviation for water-filled calibration rod (Hounsfield Units)
#' @param densities numeric vector of known cal rod densities. Format must be c(air, water, Si, glass)
#' 
#' @return value \code{conv} returns a dataframe with one row per CT slice. Values returned are the area and volume of 7 material classes: gas, peat, roots and rhizomes, rock and shell, fine mineral particles, sand, and water.
#' 
#' @seealso \code{\link{rootSize}} operates similarly.
#' 
#' @examples
#' ### Not run:
#' \dontrun{
#' data(core_426)
#' ct.slope <- unique(extractHeader(core_426$hdr, "RescaleSlope"))
#' ct.int   <- unique(extractHeader(core_426$hdr, "RescaleIntercept")) 
#' # convert raw units to Hounsfield units
#' HU_426 <- lapply(core_426$img, function(x) x*ct.slope + ct.int)
#' 
#' materials <- conv(HU_426)
#' 
#' plot using "ggplot" package after transforming with "reshape2" package
#' mass.long <- melt(materials, id.vars = c("depth"), measure.vars = grep(".g", names(rootChars2)))
#' ggplot(data = mass.long, aes(y = -depth, x = value, color = variable)) + 
#' geom_point() + theme_classic() + xlab("mass per section (g)")
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
conv <- function(mat.list, upperLim = 3045, lowerLim = -1024,
                 pixelA, thickness = 0.625, # all in mm
                 airHU = -850.3233,
                 SiHU = 271.7827,
                 glassHU = 1345.0696,
                 waterHU = 63.912,
                 waterSD = 14.1728,
                 densities = c(0.0012, 1, 1.23, 2.2) # format = air, water, Si, glass
) {
  pb <- txtProgressBar(min = 0, max = length(mat.list), initial = 0, style = 3)
  voxelVol <- pixelA * thickness / 1e3 # cm3
  water.LB <- waterHU - waterSD
  water.UB <- waterHU + waterSD
  splits <- data.frame(material = c("air", "R&R", "water", "peat", "particles", "sand", "rock_shell"),
                       lower = c(round(lowerLim), round(airHU), round(water.LB), round(water.UB), round(SiHU), 750, round(glassHU)), 
                       upper = c(round(airHU), round(water.LB), round(water.UB), round(SiHU), 750, round(glassHU), round(upperLim)))
  
  densitydf <- data.frame(HU = c(airHU, waterHU, SiHU, glassHU), density = densities)
  summary(lm1 <- lm(density ~ HU, data = densitydf)) # density in g/cm3
  
  for (i in 1:length(mat.list)) {
    depth <- thickness * i / 10 # cm
    temp <- as.vector(mat.list[[i]])
    # convert from HU to g/cm3
    temp <- temp[(temp > lowerLim) & (temp < upperLim)] # subset based on upper and lower limits
    bin <- cut(temp, breaks = c(splits$lower, upperLim), labels = splits$material, right = FALSE)
    
    if (length(temp) > 0) {
      temp.output.int <- table(bin) * pixelA / 1e2 # number of pixels * pixel area = area in class (cm2)
      temp.output <- data.frame(t(as.vector(temp.output.int))) # cm2
      names(temp.output) <- paste0(names(temp.output.int), ".cm2")
      temp.output$tot.cm2 <- rowSums(temp.output[, 1:7])
      temp.output$depth <- depth
      
      vol.output <- temp.output[, 1:8] * (thickness / 10) # cm3
      names(vol.output) <- gsub("2", "3", names(temp.output)[1:8])
      
      wetMass <- (temp * coef(lm1)[2] + coef(lm1)[1]) * voxelVol  # convert to g/cm3 and then to g (wet) in each pixel
      df1 <- data.frame(bin = bin, wetMass = wetMass)
      test <- aggregate(wetMass ~ bin, data = df1, sum)  # sum mass by category
      test <- merge(data.frame(bin = splits$material), test, all = TRUE)
      mass.output <- data.frame(t(as.vector(test[, 2])))
      mass.output[is.na(mass.output)] <- 0
      names(mass.output) <- paste0(test[, 1], ".g")
      
      outDat.init <- do.call(cbind, list(temp.output, vol.output, mass.output))
    } else {
      outDat.init <- outDat[1, ]
      outDat.init$depth <- depth
    }
    
    if (i > 1) {
      outDat <- rbind(outDat, outDat.init)
    } else if (i == 1) {
      outDat <- outDat.init
    }
    setTxtProgressBar(pb, i)
  }
  outDat <- outDat[, c("depth", names(outDat)[!names(outDat) %in% "depth"])]
  outDat
}
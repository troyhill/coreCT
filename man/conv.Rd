% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/conv.R
\name{conv}
\alias{conv}
\title{Convert a matrix of semi-processed DICOM images to mass and volume of material classes. This is a deprecated version.}
\usage{
conv(mat.list, upperLim = 3045, lowerLim = -1025, 
pixelA, thickness = 0.625, # all in mm 
airHU = -850.3233, airSD = 77.6953, 
SiHU = 271.7827, SiSD = 39.2814,
glassHU = 1345.0696, glassSD = 45.4129,
waterHU = 63.912, waterSD = 14.1728,
densities = c(0.0012, 1, 1.23, 2.2))
}
\arguments{
\item{mat.list}{list of DICOM images for a sediment core (values in Hounsfield Units)}

\item{upperLim}{upper bound cutoff for pixels (Hounsfield Units)}

\item{lowerLim}{lower bound cutoff for pixels (Hounsfield Units)}

\item{pixelA}{pixel area (mm2)}

\item{thickness}{CT image thickness (mm)}

\item{airHU}{mean value for air-filled calibration rod (Hounsfield Units)}

\item{airSD}{standard deviation for air-filled calibration rod}

\item{SiHU}{mean value for colloidal silica calibration rod}

\item{SiSD}{standard deviation for colloidal Si calibration rod}

\item{glassHU}{mean value for glass calibration rod}

\item{glassSD}{standard deviation for glass calibration rod}

\item{waterHU}{mean value for water filled calibration rod}

\item{waterSD}{standard deviation for water filled calibration rod}

\item{densities}{numeric vector of known cal rod densities. Format must be c(air, water, Si, glass)}
}
\value{
value \code{conv} returns a dataframe with one row per CT slice. Values returned are the average Hounsfield Unit value, the area (cm2), volume (cm3), and mass (grams) of 7 material classes: gas, peat, roots and rhizomes, particulates, sand, water, and rock/shell. If <code>rootData = TRUE</code>, data for specified root size classes are also returned. See <code>rootSize</code> for more detail on those values.
}
\description{
Converts raw CT units to material classes for each CT slice, directly replicating Earl Davey's manual classification approach. This method is deprecated as of coreCT version 1.3.0.
}
\details{
Calculates average Hounsfield units, cross-sectional areas (cm2), volumes (cm3), and masses (g) of material classes for each CT slice. This function assumes that core walls and all non-sediment material have been removed from the raw DICOM imagery. This function converts data from raw x-ray attenuation values to Hounsfield Units, and then uses user-defined calibration rod inputs to categorize sediment components: air, roots and rhizomes, peat, water, particulates, sand, and rock/shell.
}
\examples{
ct.slope <- unique(extractHeader(core_426$hdr, "RescaleSlope"))
ct.int   <- unique(extractHeader(core_426$hdr, "RescaleIntercept")) 
# convert raw units to Hounsfield units
HU_426 <- lapply(core_426$img, function(x) x*ct.slope + ct.int)

materials <- conv(HU_426, pixelA = 0.0596)

\dontrun{
# plot using "ggplot" package after transforming with "reshape2" package
mass.long <- reshape2::melt(materials, id.vars = c("depth"), 
   measure.vars = grep(".g", names(materials)))
ggplot2::ggplot(data = mass.long, ggplot2::aes(y = -depth, x = value, 
   color = variable)) + ggplot2::geom_point() + ggplot2::theme_classic() + 
   ggplot2::xlab("mass per section (g)")
}

}
\seealso{
\code{\link{rootSize}} operates similarly.
}

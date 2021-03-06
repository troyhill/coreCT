% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/voxDims.R
\name{voxDims}
\alias{voxDims}
\title{Extract voxel dimensions from DICOM image}
\usage{
voxDims(directory = file.choose())
}
\arguments{
\item{directory}{a character string that can be a matrix of DICOM images or the address of an individual DICOM file in a folder of DICOM images. The default action is <code>file.choose()</code>; a browser menu appears so the user can select the the desired directory by identifying a single DICOM file in the folder of images.}
}
\value{
value \code{voxDims} returns a two-column dataframe showing the pixel area and slice thickness. Values in the DICOM headers are assumed to be millimeters; pixel area and slice thickness columns are labeled based on this assumption.
}
\description{
Extract pixel area and slice thickness from DICOM header to characterize voxel (3D pixel) dimensions.
}
\examples{
# data(core_426)
voxDims("core_426")

}

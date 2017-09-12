[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/coreCT)](https://cran.r-project.org/package=coreCT)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.889651.svg)](https://doi.org/10.5281/zenodo.889651)

[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/coreCT)](https://cran.r-project.org/package=coreCT)




### About **coreCT**

**coreCT** is an R package for programmatic analysis of sediment cores that have been digitized by computed tomography. The package converts Hounsfield Units to material classes (e.g., peat, root/rhizome, sand) and quantifies component masses and volumes. To get started quickly, check out the package vignette:

  - [vignette: 	Using the coreCT package](https://cran.r-project.org/web/packages/coreCT/vignettes/Using_coreCT.html)


### Sediment core characterization

**coreCT** output can be easily plotted using the **reshape2** and **ggplot2** packages (figures below were produced using the same code included in the help file examples).


<img src="https://raw.githubusercontent.com/troyhill/images/master/221_20160607_Vol.png" width="400" height="500" />


**Figure 1.** Volumes of various compartments in a sediment core



**coreCT** also quantifies the number of root/rhizome particles in a user-defined range of size classes, and calculates the external surface area and volume attributable to each size class. This allows the estimation of belowground plant organ contributions to total soil volume and, dividing by core area, elevation.


<img src="https://raw.githubusercontent.com/troyhill/images/master/221_20160607_Particles.png" width="400" height="500" /> <img src="https://raw.githubusercontent.com/troyhill/images/master/221_20160607_rootVol.png" width="400" height="500" />


**Figure 2.** Number (left) and combined volume (right) of root/rhizome particles in each of four size classes.




### Installing **coreCT**

**coreCT** is available on CRAN and can be installed by typing:
    install.packages("coreCT") 
    

To install and load the development version of coreCT, run the following commands in R:

    install.packages("devtools") # if devtools isn't already installed

    devtools::install_github("troyhill/coreCT")

    library(coreCT)


# coreCT

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/coreCT)](https://cran.r-project.org/package=coreCT) [![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/coreCT)](https://cran.r-project.org/package=coreCT) [![](http://cranlogs.r-pkg.org/badges/grand-total/coreCT)](https://cran.r-project.org/package=coreCT) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.889651.svg)](https://doi.org/10.5281/zenodo.889651) 
[![Build Status](https://travis-ci.org/troyhill/coreCT.svg?branch=master)](https://travis-ci.org/troyhill/coreCT) [![codecov.io](https://codecov.io/github/troyhill/coreCT/coverage.svg?branch=master)](https://codecov.io/github/troyhill/coreCT?branch=master)

  
    

## Installing **coreCT**

**coreCT** is available on CRAN and can be installed from the R console:

    install.packages("coreCT") 
    


Once installed, the package can be loaded:
    
    library(coreCT)

  
  

## About **coreCT**

**coreCT** is an [R package](https://cran.r-project.org/package=coreCT) for programmatic analysis of sediment cores that have been digitized by computed tomography. The package converts Hounsfield Units to material classes (e.g., peat, root/rhizome, sand) and quantifies component masses and volumes. To get started quickly, check out the package vignette.


## Sediment core characterization

**coreCT** output can be easily plotted using the **reshape2** and **ggplot2** packages (figures below were produced using the same code included in the help file examples).


<img src="https://raw.githubusercontent.com/troyhill/images/master/221_20160607_Vol.png" width="400" height="500" />


**Figure 1.** Volumes of various compartments in a sediment core



**coreCT** also quantifies the number of root/rhizome particles in a user-defined range of size classes, and calculates the external surface area and volume attributable to each size class. This allows the estimation of belowground plant organ contributions to total soil volume and, dividing by core area, elevation.


<img src="https://raw.githubusercontent.com/troyhill/images/master/221_20160607_Particles.png" width="400" height="500" /> <img src="https://raw.githubusercontent.com/troyhill/images/master/221_20160607_rootVol.png" width="400" height="500" />


**Figure 2.** Number (left) and combined volume (right) of root/rhizome particles in each of four size classes.





## License and disclaimer

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program.  If not, see <https://www.gnu.org/licenses/>.

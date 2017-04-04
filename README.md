# coreCT

'coreCT' is an R package that streamlines analysis of CT-scanned sediment cores. 



To install and load coreCT, run the following commands in R:

    install.packages("devtools") # if devtools isn't already installed

    devtools::install_github("troyhill/coreCT")

    library(coreCT)


### About coreCT

coreCT programmatically ananlyzes CT-scanned sediment cores, to rapidly convert Hounsfield Units to material classes (peat, root/rhizome, sand) and quantify their masses and volumes. 

<img src="https://raw.githubusercontent.com/troyhill/images/master/221_20160607_Vol.png" width="200" height="200" />
Figure 1. Volumes of various compartments in a sediment core

coreCT also quantifies the number of root/rhizome particles in a user-defined range of size classes, and calculates the external surface area and volume attributable to each size class. This allows the estimation of belowground plant organ contributions to total soil volume and, dividing by core area, elevation.
<img src="https://raw.githubusercontent.com/troyhill/images/master/221_20160607_particles.png" width="200" height="200" /> <img src="https://raw.githubusercontent.com/troyhill/images/master/221_20160607_rootVol.png" width="200" height="200" />
Figure 2. Number (left) and combined volume (right) of root/rhizome particles in each of four size classes.

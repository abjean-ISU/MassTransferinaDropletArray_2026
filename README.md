This repository contains the data and code relavent to the manuscript Predicting Mass Transfer in a Vertical Droplet Array, submitted to AIChE Journal. 

Structure 
Folders:   
  Droplet Tower Modeling Files
    Modeling files used to predict mass transfer and liquid holdup in the experimental apparatus. 
    Running NDFApproach_dropletTowerModel_09242024PDF.m requires the use of the AllDropletData.csv available in Experimental Data, which is extracted from droplet images. 
    ODE files are for each of the 4 analyzed kL models
    PD-gqmom_0_1.m was shared by Alberto Passalacqua and should be cited: 
      Fox, R. O., Laurent, F., & Passalacqua, A. (2023). The generalized quadrature method of moments. Journal of Aerosol Science, 167, 106096. https://doi.org/10.1016/j.jaerosci.2022.106096
      
  Experimental Data
    Data on mass transfer, liquid holdup, and droplet array characteristics. The mass transfer and liquid holdup data are raw, and the droplet data are available in 2 forms: 
      Cleaned droplet data for each of the 24 locations in the droplet array where images were collected, and the composite file of all droplet data used to compute average array characteristics. 

  Image Processing Files
    Matlab files used to process droplet images and clean stacks of data to produce AllDropletData.csv

  Model Results
    Droplet tower model results for 18 experimental trials 

  DataProcessingCodeShare.R: code used to compute experimental reults shared in manuscript from both raw data and modeling results. 

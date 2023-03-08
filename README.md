# Rschistox

This is a package for modelling schistosomiasis treatment in low to high prevalence settings. 
The model can be used for both *S.haematobium* and S. mansoni, and also for treatment of:
- Pre school aged children (pre-SAC)
- School aged children (SAC)
- Adult population

The package is developed from the Julia package Schistoxpkg (https://github.com/mattg3004/Schistoxpkg.jl) developed by Matt Graham.

The model allows the user to specify the proportion of population receiving treatment of praziquantel to the three groups mentioned above, and the efficacy of the treatment.
It is also modelled to accomodate vaccination, once available for use.

To install the package you need to have devtools installed. If this isn't installed, then install it with:
```R
> install.packages("devtools")
```
After this, you can run the following to install the package:
```R
> library(devtools) 
> devtools::install_github('mutono/Rschistox', build_vignettes=T)
```

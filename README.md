# metaIntegrator #

This is README space for metaIntegrator_public. To be filled with a README if necessary, or for internal use if required.

### What is this repository for? ###

*  This repository is for building and testing the metaIntegrator_public R package

### Who do I talk to? ###

* Either create an issue, email adityamr [at] stanford.edu, or email hayneswa [at] stanford.edu

### How do I install the MetaIntegrator package? ###

For the current version on CRAN:

`install.packages("MetaIntegrator")`

For the current development version:

Run the following lines in R:

`library(devtools)`

`install_bitbucket("khatrilab/metaIntegrator_public", auth_user=USERNAME, password=PASSWORD)`

where USERNAME and PASSWORD are your bitbucket username and password.  
  
  
If you get "ERROR: dependency ‘deapathways’ is not available for package ‘MetaIntegrator’", follow these steps:  

download the `deapathways_1.0.tar.gz` file (which is stored in metaIntegrator_public/inst).

install deapathways, either by going to the command line, navigating to the folder containing  
the `tar.gz` file in the command line and running `R CMD INSTALL deapathways_1.0.tar.gz`,  
or by running `install.packages(PATH_TO_FILE, repos = NULL, type="source")`  
where `PATH_TO_FILE` is the path to the `deapathways_1.0.tar.gz` file  


### Documentation ###
To read the Vignette go to: http://khatrilab.bitbucket.org/MetaIntegratorVignette/
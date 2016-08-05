# metaIntegrator #

This is README space for metaIntegrator_public. To be filled with a README if necessary, or for internal use if required.

### What is this repository for? ###

*  This repository is for building and testing the metaIntegrator_public R package

### Who do I talk to? ###

* Either create an issue or email me: hayneswa [at] stanford.edu

### How do I install the MetaIntegrator package? ###

Method #1 (preferred):

Run the following lines in R:

library(devtools)

install_bitbucket("khatrilab/metaintegrator_public", auth_user=USERNAME, password=PASSWORD)

where USERNAME and PASSWORD are your bitbucket username and password.

Method #2:

Download the .tar.gz of the current snapshot of MetaIntegrator, available at: 

https://bitbucket.org/khatrilab/metaintegrator_public/downloads/metaintegrator_public.tar.gz

On your computer, navigate to the directory of the downloaded file. 

Run the following lines in R:

install.packages("metaintegrator_public.tar.gz", repos=NULL, type="source")

### Documentation ###
To read the Vignette go to: http://khatrilab.bitbucket.org/MetaIntegratorVignette/
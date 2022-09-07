## ################################
##
## Make Package Codes
##
## ################################

## Remove pkg 
remove.packages("ppfunreg")

## Create/update documentation and (re-)write NAMESPACE
devtools::document("ppfunreg")

## CRAN-check pkg
# devtools::check("ppfunreg")       # check the package

## Install
devtools::install_local("ppfunreg", force = TRUE)
##
## #################################

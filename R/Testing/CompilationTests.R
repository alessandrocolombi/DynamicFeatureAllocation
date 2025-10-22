wd = "C:/Users/colom/DynamicFeatureAllocation/R/Testing/"
wd = "/home/lucia.paci/Lucia/Ale/DynamicFeatureAllocation/R/Testing/"
setwd(wd)
Rcpp::sourceCpp("../../src/RcppFunctions.cpp")
log_dnorm(0,0,1)

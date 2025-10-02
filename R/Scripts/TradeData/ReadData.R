wd = "C:/Users/colom/DynamicFeatureAllocation/R/Scripts/TradeData/"
setwd(wd)
list.files()
library(readxl)
library(tidyverse)
library(tibble)

# Country names, ID and ISO code
CountryDictionary <- read_excel("TradeData.xlsm", 
                        sheet = "Country List", range = "B2:D197", 
                        col_types = c("text", "numeric", "text"))
CountryDictionary[] <- lapply(CountryDictionary, factor)
# save(CountryDictionary,file = "CountryDictionary.rda")
# load("CountryDictionary.rda")
p = nrow(CountryDictionary)

# Read country trades
# Read the needed rows (1:37831) and all useful columns
data_raw <- read_excel("TradeData.xlsm", 
                       sheet = "Data Sheet", 
                       range = cell_limits(c(1, 2), c(37831, 78))) 
# Explanation: (1,2) = row 1 col B, (37831,78) = row 37831 col BZ (BZ is col 78)

# Now select the exact columns you want by Excel letters
data <- data_raw %>% select("IDEX","ISOEX","IDIM","ISOIM",
                            "1995","1996","1997","1998","1999","2000",             
                            "2001","2002","2003","2004","2005",             
                            "2006","2007","2008","2009","2010",             
                            "2011","2012","2013","2014","2015")

head(data)

for (j in 5:ncol(data)) {
  noc_count <- sum(data[[j]] == "NoCty", na.rm = TRUE)
  cat("Column", names(data)[j], "has", noc_count, "NoCty entries\n")
  
  # Replace "NoCty" with 0
  data[[j]][data[[j]] == "NoCty"] <- 0
}
nyears <- ncol(data) - 4
data <- data %>% mutate(across(-c(2, 4), as.integer))  # transform all except 2nd and 4th

# Transform data into a list of matrices
data_list = vector("list",length = nyears)
for(ii in 5:ncol(data)){
  Mat_t = matrix(0,nrow = p, ncol = p)
  for(hh in 1:nrow(data)){
    Mat_t[ data[[1]][hh], data[[3]][hh] ] = data[[ii]][hh]
  }
  data_list[[ii-4]] = Mat_t
}

# Transform data into a 3-dim. array
data_array <- array(data = unlist(data_list), dim = c(p, p, nyears))

# save(data_list, file = "TradeData_Asym_list.rda")
# save(data_array, file = "TradeData_Asym_array.rda")


# Save data in symmetric form
data_list_sym = lapply(data_list, function(x){
  y = x + t(x)
  y[which(y >= 1)] = 1
  y
})
# Transform data into a 3-dim. array and save
data_array_sym <- array(data = unlist(data_list_sym), dim = c(p, p, nyears))
 
# save(data_list_sym, file = "TradeData_Sym_list.rda")
# save(data_array_sym, file = "TradeData_Sym_array.rda")



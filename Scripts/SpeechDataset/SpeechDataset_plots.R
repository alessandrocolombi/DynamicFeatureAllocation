# Required packages
library(quanteda)
library(tidyverse)
library(tidytext)
wd = "C:/Users/colom/DynamicFeatureAllocation/Scripts/SpeechDataset"
setwd(wd)


# Read President names
Presidents_all <- read.csv("C:/Users/colom/DynamicFeatureAllocation/Scripts/data/Presidents_all.csv")
Presidents_all[,2]

save_img = FALSE
save_name_base = "img/SpeechData_mat_"

# Read all data -----------------------------------------------------------

data = read.table(paste0("../data/SpeechData_all.txt"))
colnames(data) = Presidents_all[,2]

# Dimensions
V = nrow(data)
Ttot = ncol(data)
cat("\n Vocabulary size: ",V,"; Number of documents: ",Ttot,"\n")

mycol = c(rep("white",1), # 0
          rep("#A6D59D",1), # 1
          rep("#004616",9), # 2:10
          rep("darkred",100)  # > 10
          )

if(save_img)
  pdf( paste0(save_name_base,"all",".pdf") )
par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
image( 1:Ttot, 1:V, 
       t(data),   
       col = mycol,    
       xlab = "Time", 
       ylab = "Words",
       main = "Data - full",
       axes = FALSE )
axis(1, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
     labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
     cex.axis = 0.7 )
axis(2, at = seq(1, V, length.out = min(V, 10)), 
     labels = round(seq(1, V, length.out = min(V, 10))),
     cex.axis = 0.7)
box()
fields::image.plot(
  1:Ttot, 1:V, t(data),
  col = mycol,
  legend.only = TRUE,
  horizontal = FALSE,
  legend.width = 1.2,            # controls legend thickness
  legend.shrink = 0.8,           # smaller legend
  legend.mar = 8.5,                # margin from image
  legend.args = list(text = " ", side = 3, line = 1, cex = 0.8)
)
if(save_img)
  dev.off()


# Num. words per speech
Nt = colSums(data)
pdf(NULL)
  bp1 <- barplot(height = Nt)
dev.off()
xlabs = sapply(names(Nt), substr,1,4)
ylabs = round(seq(0, max(Nt), by=500) * 1e-3,1)

# Saved in custom size: 6x14 in
if(save_img)
  pdf( paste0("img/SpeechData_Nt_","all.pdf") )
par(mfrow = c(1,1), mgp=c(2,0.5,0), mar = c(3,3,1,0))
barplot( height = Nt, 
         names.arg = "", las = 2, col = "darkblue", border = NA,
         main = " ", xlab = "Year", ylab = "#words. (x10^-3)", yaxt = "n" )
axis( side = 2, at = ylabs*1e3, labels = ylabs, las = 1)
text( x = bp1, y = par("usr")[3] - 0.02*max(Nt), 
      labels = xlabs, srt = 45, adj = 1, xpd = TRUE, cex = 0.5 )
if(save_img)
  dev.off()



# Words freq.
Ni = rowSums(data)
Ni = sort(Ni, decreasing = TRUE)
# Ni = Ni[1:100]
# Ni = Ni[which(Ni > 10)]
# pdf(NULL)
#   bp1 <- barplot(height = Ni)
# dev.off()
xlabs = sapply(names(Ni), substr,1,4)
ylabs = round(seq(0, max(Ni), by=100),1)

# Saved in custom size: 6x14 in
if(save_img)
  pdf( paste0("../img/SpeechData_Ni_","all.pdf") )
par(mfrow = c(1,1), mgp=c(2,0.5,0), mar = c(3,3,1,0))
barplot( height = Ni, 
         names.arg = "", las = 2, col = "darkblue", border = NA,
         main = " ", xlab = "Word", ylab = "#Occur.", yaxt = "n" )
axis( side = 2, at = ylabs, labels = ylabs, las = 1)
# text( x = bp1, y = par("usr")[3] - 0.02*max(Ni), 
#       labels = xlabs, srt = 45, adj = 1, xpd = TRUE, cex = 0.5 )
if(save_img)
  dev.off()



colnames(data) = Presidents_all[,2]
vocabulary = rownames(data)

# Saved in custom size: 6x14 in
if(save_img)
  pdf( paste0("../img/SpeechData_AllPres_","all.pdf"))
par(mfrow = c(3,3), mgp=c(2,0.5,0), mar = c(3,3,1,0))
for(t in 1:60){
  ft = data[,t]
  ylabs = seq(0,30,by = 10)
  barplot( height = ft, 
           names.arg = "", las = 2, col = "darkblue", border = NA,
           main = Presidents_all[t,2], xlab = "Word", ylab = "#Occur.", yaxt = "n",
           ylim = c(0,30))
  axis( side = 2, at = ylabs, labels = ylabs, las = 1)
  # text( x = bp1, y = par("usr")[3] - 0.02*max(ft), 
  #       labels = vocabulary, srt = 45, adj = 1, xpd = TRUE, cex = 0.5 )
}
if(save_img)
  dev.off()


# Read top r data ---------------------------------------------------------

rgrid = c(3,10,20,50,100,200)

# Matrici
for(r in rgrid){
  data = read.table(paste0("../data/SpeechData_top",r,".txt"))
  colnames(data) = Presidents_all[,2]
  
  # Dimensions
  V = nrow(data)
  Ttot = ncol(data)
  
  mycol = c(rep("white",1), # 0
            rep("#A6D59D",1), # 1
            rep("#004616",9), # 2:10
            rep("darkred",100)  # > 10
  )
  
  if(save_img)
    pdf( paste0(save_name_base,"top",r,".pdf") )
  par(mfrow = c(1,1), mar = c(3.5,3.5,2,8), mgp=c(2,0.5,0))
  image( 1:Ttot, 1:V, 
         t(data),   
         col = mycol,    
         xlab = "Time", 
         ylab = "Words",
         main = paste0("Data - top ",r),
         axes = FALSE )
  axis(1, at = seq(1, Ttot, length.out = min(Ttot, 10)), 
       labels = round(seq(1, Ttot, length.out = min(Ttot, 10))),
       cex.axis = 0.7 )
  axis(2, at = seq(1, V, length.out = min(V, 10)), 
       labels = round(seq(1, V, length.out = min(V, 10))),
       cex.axis = 0.7)
  box()
  fields::image.plot(
    1:Ttot, 1:V, t(data),
    col = mycol,
    legend.only = TRUE,
    horizontal = FALSE,
    legend.width = 1.2,            # controls legend thickness
    legend.shrink = 0.8,           # smaller legend
    legend.mar = 8.5,                # margin from image
    legend.args = list(text = " ", side = 3, line = 1, cex = 0.8)
  )
  if(save_img)
    dev.off()
  
}

# Plot Nt
for(r in rgrid){
  data = read.table(paste0("../data/SpeechData_top",r,".txt"))
  colnames(data) = Presidents_all[,2]
  # Num. words per speech
  Nt = colSums(data)
  pdf(NULL)
    bp1 <- barplot(height = Nt)
  dev.off()
  
  xlabs = sapply(names(Nt), substr,1,4)
  ylabs = round(seq(0, max(Nt), by=500) * 1e-2,1)
  
  # Saved in custom size: 6x14 in
  if(save_img)
    pdf( paste0("img/SpeechData_Nt_","top_",r,".pdf") )
  par(mfrow = c(1,1), mgp=c(2,0.5,0), mar = c(3,3,1,0))
  barplot( height = Nt, 
           names.arg = "", las = 2, col = "darkblue", border = NA,
           main = paste0("r = ",r), xlab = "Year", ylab = "#words. (x10^-2)", yaxt = "n" )
  axis( side = 2, at = ylabs*1e2, labels = ylabs, las = 1)
  text( x = bp1, y = par("usr")[3] - 0.02*max(Nt), 
        labels = xlabs, srt = 45, adj = 1, xpd = TRUE, cex = 0.5 )
  if(save_img)
    dev.off()

}

# Plot Ni
for(r in rgrid){
  data = read.table(paste0("../data/SpeechData_top",r,".txt"))
  colnames(data) = Presidents_all[,2]
  # Words freq.
  Ni = rowSums(data)
  Ni = sort(Ni, decreasing = TRUE)
  # Ni = Ni[1:100]
  pdf(NULL)
    bp1 <- barplot(height = Ni)
  dev.off()
  xlabs = sapply(names(Ni), substr,1,4)
  ylabs = round(seq(0, max(Ni), by=100),1)
  
  # Saved in custom size: 6x14 in
  if(save_img)
    pdf( paste0("img/SpeechData_Ni_","top_",r,".pdf") )
  par(mfrow = c(1,1), mgp=c(2,0.5,0), mar = c(3,3,1,0))
  barplot( height = Ni, 
           names.arg = "", las = 2, col = "darkblue", border = NA,
           main = paste0("r = ",r), xlab = "Word", ylab = "#Occur.", yaxt = "n" )
  axis( side = 2, at = ylabs, labels = ylabs, las = 1)
  text( x = bp1, y = par("usr")[3] - 0.02*max(Ni), 
        labels = xlabs, srt = 45, adj = 1, xpd = TRUE, cex = 0.5 )
  if(save_img)
    dev.off()
}


# Speech specific distributions -------------------------------------------


r = 20
data = read.table(paste0("../data/SpeechData_top",r,".txt"))
colnames(data) = Presidents_all[,2]
vocabulary = rownames(data)

# Saved in custom size: 6x14 in
if(save_img)
  pdf( paste0("../img/SpeechData_AllPres_","top_",r,".pdf") )
par(mfrow = c(3,3), mgp=c(2,0.5,0), mar = c(3,3,1,0))
for(t in 1:ncol(data)){
  ft = data[,t]
  # pdf(NULL)
  #   bp1 <- barplot(height = ft)
  # dev.off()
  # ylabs = round(seq(0, max(ft), by=10),1)
  ylabs = c(0,10,20,30)

  barplot( height = ft, 
           names.arg = "", las = 2, col = "darkblue", border = NA,
           main = Presidents_all[t,2], xlab = "Word", ylab = "#Occur.", yaxt = "n",
           ylim = c(0,30))
  axis( side = 2, at = ylabs, labels = ylabs, las = 1)
  # text( x = bp1, y = par("usr")[3] - 0.02*max(ft), 
  #       labels = vocabulary, srt = 45, adj = 1, xpd = TRUE, cex = 0.5 )
}
if(save_img)
  dev.off()


# Brutta ------------------------------------------------------------------


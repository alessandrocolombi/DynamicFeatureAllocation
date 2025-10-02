wd = "C:/Users/colom/DynamicFeatureAllocation/R/Scripts/TradeData/"
setwd(wd)
list.files()
library(readxl)
library(tidyverse)
library(tibble)
library(igraph)
load("CountryDictionary.rda")
load("TradeData_Sym_list.rda")
p = nrow(CountryDictionary)

## Plot 
save_plot = TRUE
nome_img = "img/TradeData9521_networks.pdf"

# vertex_labs = as.character(CountryDictionary$ISO)
vertex_labs = as.character(CountryDictionary$ID)

if(save_plot)
  pdf(nome_img)
for(ii in 1:length(data_list_sym)){
  g <- graph_from_adjacency_matrix(data_list_sym[[ii]], mode = "undirected", diag = FALSE)
  par(mar = c(0,0,0,0))   # no margins around the plot
  plot(g,
       vertex.size = 15,
       vertex.label = vertex_labs,
       edge.color = "grey85",
       vertex.color = "skyblue",
       vertex.label.cex = 0.8,
       vertex.label.color = "white",
       vertex.frame.color = NA,
       layout = layout_with_kk(g))   
  
}
if(save_plot)
  dev.off()








# Brutta (comm. detection) ------------------------------------------------


cl = cluster_fast_greedy(g)
memb = cl$membership
k <- length(unique(memb))
centers <- layout_in_circle(make_full_graph(k)) * 10

layout_manual <- matrix(NA, nrow = p, ncol = 2)
for (i in 1:k) {
  nodes <- which(memb == i)
  sublay <- layout_in_circle(make_full_graph(length(nodes)))
  layout_manual[nodes, ] <- sweep(sublay, 2, centers[i, ], "+")
}

mycol = c("mediumorchid4", "firebrick","royalblue4","gold1", 
          "brown4", "salmon1","darkslategray3", "tomato4")


# Set nodes sizes, colors and labels
vertex_sizes <- rep(15, p)                 
vertex_colors <- memb
for(i in 1:k){
  vertex_colors[vertex_colors == i] = mycol[i]
}
vertex_labs = as.character(CountryDictionary$ISO)


# Plot the graph
par(mar = c(0,0,0,0)) 
plot(g,
     vertex.size = vertex_sizes,
     vertex.color = vertex_colors,
     vertex.label = vertex_labs,
     vertex.label.cex = 0.8,
     vertex.label.color = "white",
     vertex.frame.color = NA,
     edge.color = "darkgrey",
     edge.width = 2,
     layout = layout_manual,
     main = " ")








g <- graph_from_adjacency_matrix(data_list_sym[[1]], mode = "undirected", diag = FALSE)
cl = cluster_fast_greedy(g)
memb = cl$membership
ord <- order(memb)
g_reordered <- permute(g, ord)
A <- as_adjacency_matrix(g_reordered, sparse = FALSE)







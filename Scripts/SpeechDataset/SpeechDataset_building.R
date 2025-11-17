# Required packages
library(quanteda)
library(tidyverse)
library(tidytext)
wd = "C:/Users/colom/DynamicFeatureAllocation/Scripts/SpeechDataset"
setwd(wd)

# 1) load corpus
data("data_corpus_inaugural")    # built-in in quanteda
corp <- data_corpus_inaugural

# Quick check
ndocs <- ndoc(corp)             # should be 60
cat("Number of documents (T):", ndocs, "\n")

# 2) parameters and custom stopwords
r <- 3   # top-r per document 

# FANBOYS (coordinating conjunctions) + some domain words you may want to drop
fanboys <- c("for", "and", "nor", "but", "or", "yet", "so")  # the classic FANBOYS
extra_stop <- c() 
# combine with quanteda's English stopwords
my_stopwords <- c(stopwords("en"), fanboys, extra_stop) |> unique()

# 3) Tokenize and basic normalization
toks <- tokens(corp,
               remove_punct = TRUE,
               remove_numbers = TRUE,
               remove_symbols = TRUE) |>
  tokens_tolower() |>
  tokens_remove(my_stopwords) 
# Stemming: applies a word stemmer (e.g., "running", "runs", "run" -> "run").
toks <- tokens_wordstem(toks, language = "english")
# Include bigrams: from "united" and "states" to "united_states"
# toks <- tokens_ngrams(toks, n = 1:2)

# 4) Build dfm (document-feature matrix) - counts per doc
dfm_all <- dfm(toks)
dfm_all
# dfm_all <- dfm_sort(dfm_all, margin = "features", decreasing = TRUE)
# dfm_all

# 5) For each document: get top-r tokens (by count inside that doc)
# We'll iterate over documents and store the top tokens
docnames_vec <- docnames(dfm_all)

top_tokens_per_doc <- lapply(docnames_vec, function(dn) {
  # extract feature counts for document dn
  fv <- as.numeric(dfm_all[dn, ])
  names(fv) <- featnames(dfm_all)
  # sort descending and take top r (or fewer if doc shorter)
  top_k <- sort(fv, decreasing = TRUE)
  top_feats <- head(names(top_k), n = r)       # will be <= r if doc has fewer types
  return(top_feats)
})
names(top_tokens_per_doc) <- docnames_vec

# 6) Union across documents -> vocabulary V
vocab_union <- unique(unlist(top_tokens_per_doc))
V <- length(vocab_union)
cat("Vocabulary size after union (V):", V, "\n")

# 7) Build matrix D (rows = words in vocab_union, cols = documents), counts
# dfm_sub <- dfm_select(dfm_all, pattern = vocab_union, selection = "keep", valuetype = "fixed")
dfm_sub <- dfm_all[,vocab_union]

# convert to matrix with rows = words, cols = docs:
mat_counts <- t(as.matrix(dfm_sub)) # now rows = vocab (words), cols = documents

# set rownames and colnames explicitly
rownames(mat_counts) <- vocab_union
colnames(mat_counts) <- docnames_vec

# 8) Quick checks
cat("Dimensions of D: rows (V) =", nrow(mat_counts), ", cols (T) =", ncol(mat_counts), "\n")

# 9) Example: inspect first 10 words and first 6 documents
# mat_counts[which( rownames(mat_counts) == "god" ),]
# x <- as.numeric(dfm_all[, which(colnames(dfm_all) == "god")])
# names(x) <- docnames(dfm_all)
# x   


# 10) Save for later
# write.csv(colnames(mat_counts), "data/Presidents_all.csv")
write.table(mat_counts, paste0("data/SpeechData_top",r,".txt"))
cat("\n DONE \n")


# read.table("SpeechData_top10.txt")
# Presidents_all <- read.csv("C:/Users/colom/DynamicFeatureAllocation/Scripts/Presidents_all.csv")
# Presidents_all[,2]

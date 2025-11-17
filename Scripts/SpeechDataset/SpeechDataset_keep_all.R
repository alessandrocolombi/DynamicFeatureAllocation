# Required packages
library(quanteda)
library(tidyverse)
library(tidytext)

# 1) load corpus
data("data_corpus_inaugural")    # built-in in quanteda
corp <- data_corpus_inaugural

# Quick check
ndocs <- ndoc(corp)             # should be 60
cat("Number of documents (T):", ndocs, "\n")

# 2) parameters and custom stopwords
r <- 100   # top-r per document 

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

# convert to matrix with rows = words, cols = docs:
mat_counts <- t(as.matrix(dfm_all)) # now rows = vocab (words), cols = documents

# set rownames and colnames explicitly
rownames(mat_counts) <- vocab_union
colnames(mat_counts) <- docnames_vec

# 8) Quick checks
cat("Dimensions of D: rows (V) =", nrow(mat_counts), ", cols (T) =", ncol(mat_counts), "\n")

# 10) Save for later
write.table(mat_counts, paste0("../data/SpeechData_all.txt"))
cat("\n DONE \n")

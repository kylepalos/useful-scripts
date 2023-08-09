# I often analyze adjacent genes or overlapping genes and want to know whether their gene expression is correlated
# This script will walk through that process

# load libraries
library(DESeq2)
library(dplyr)
library(magrittr)
library(tibble)
library(ggplot2)
library(corrr)

# make a count matrix

# Set the random seed for reproducibility
set.seed(123)

# make an example DESeq2 dataset so I can extract the counts from them
dds <- makeExampleDESeqDataSet(n = 20000,
                               m = 12)
counts <- counts(dds, normalized=F)

# normalize all gene counts by library depth
# do not need to normalize by gene lengths because these are made up genes
depth.norm.counts <- t(t(counts) * 1e4 / colSums(counts) )

# let's only retain the 25% most variable genes 
x <- apply(depth.norm.counts, 1, mad) #Calculate median absolute deviation for each gene 
var.norm.counts <- depth.norm.counts[x>quantile(x,0.75),] # selecting top 25% most variable genes  

# This next block uses the corrr package because it allows useful data transformations while doing the correlation
counts.cor <- var.norm.counts %>%
  t() %>% #transpose the dataframe
  correlate(method = "pearson") %>%  #run the pearson correlation
  shave() %>% #remove the upper triangle of correlation values
  stretch(na.rm = TRUE) # gather all the columns

# Here is where you would select for your gene pairs of interest for subsequent analyses


# correlating a bunch of random genes with random gene expression values should producea normal distribution of 
# correlation coefficients
counts.cor %>%
  ggplot(aes(x = r)) +
  geom_density() +
  theme_classic(base_size = 16) +
  xlab("Pearson correlation coefficient")

# quite normal
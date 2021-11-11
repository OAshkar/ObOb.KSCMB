#####################################################################################################
# Script to import different TPM's from salmon and return
# 1. Counts matrix in Rdata for DGE
# 2. median table for matlab tINIT.
library(optparse)
option_list <- list(
  make_option(c("-g", "--genome"),
    type = "character", default = NULL,
    help = "Genome GTF file", metavar = "character"
  ),
  make_option(c("-d", "--directory"),
    type = "character", default = "salmon",
    help = "directory of salmon files", metavar = "character"
  ),
  make_option(c("-o", "--out"),
    type = "character", default = "counts/TPMaverage.csv",
    help = "output file name [default= %default]", metavar = "character"
  )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

print(opt)
######################################################################################################
library(tximport)
library(tidyverse)

gtf <- rtracklayer::import(opt$genome)
# gtf <- rtracklayer::import("genome/Mus_musculus.GRCm38.93.gtf")
tx2gene <- as.data.frame(gtf) %>%
  select(transcript_id, transcript_version, gene_name) %>%
  mutate(transcript_id = paste0(transcript_id, ".", transcript_version)) %>%
  select(transcript_id, gene_name)
head(tx2gene)

samples <- list.files(opt$directory)
files <- file.path("salmon", samples, "quant.sf")
names(files) <- samples

txi.salmon <- tximport(files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM"
)
counts <- as.data.frame(txi.salmon$counts)
dim(counts) # 84 samples, 35596 gene

### Raw counts in R objects for DGE
save(counts, file = "./counts/Count.Rdata")

#######################################################################################################
# Manuipulation to add media
counts <- counts %>%
  rownames_to_column("genes") %>%
  pivot_longer(
    values_to = "TPM",
    names_to = "sample", cols = -"genes"
  )
# grouping index
i <- gsub(
  x = counts$sample,
  perl = T,
  pattern = "(WT|ob_ob)(_|-)(HFD|ND|HDF)(_|-)(\\d)(_|-)(\\D\\D)(_|-)(.*)",
  replacement = "\\1_\\3_\\7"
)
ii <- gsub(
  x = i,
  pattern = "HDF",
  replacement = "HFD"
)
counts <- counts %>%
  mutate(sample = ii) %>%
  group_by(sample, genes) %>%
  summarize(average = median(TPM, , na.rm = TRUE)) %>%
  pivot_wider(
    values_from = "average",
    names_from = "sample"
  ) %>%
  column_to_rownames("genes") %>%
  mutate(
    across(everything(), ~ replace_na(.x, 0))
  )

dim(counts) # 28 sample by average, 35596 genes
head(counts)

write.csv(counts, file = opt$out)
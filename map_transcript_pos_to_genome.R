# Load required libraries
library(GenomicFeatures)
library(dplyr)
library(tibble)
library(Repitools)

# Load genome annotation file
txdb <- makeTxDbFromGFF("Athaliana.gff3")

# Create a file that maps transcript IDs (isoforms) to gene IDs
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# Generate a data frame of 10 random transcript IDs and positions within the transcript
set.seed(123) # For reproducibility
ten_transcripts <- slice_sample(tx2gene, n = 10) %>%
  select(TXNAME)

tx_posns <- data.frame(
  transcript_id = ten_transcripts$TXNAME,
  tx_position = sample(1:100, 10)
)

# Add gene IDs for each transcript ID
tx_posns <- tx_posns %>%
  left_join(tx2gene, by = c("transcript_id" = "TXNAME"))

# Get the strand for all genes and attach it to our genes of interest
tx_strand <- transcripts(txdb, columns = c("tx_id", "tx_strand"), use.names = TRUE)
tx_strand_df <- annoGR2DF(tx_strand) %>%
  rownames_to_column(var = "tx") %>%
  select(tx, tx_strand)

tx_posns <- tx_posns %>%
  left_join(tx_strand_df, by = c("transcript_id" = "tx"))

# Convert to GRanges
gr <- GRanges(
  seqnames = tx_posns$GENEID,
  ranges = IRanges(start = tx_posns$tx_position, end = tx_posns$tx_position),
  strand = tx_posns$tx_strand
)

# Isolate transcripts and genes from txdb
transcripts <- transcriptsBy(txdb, by = c("gene"))

# Use mapFromTranscripts in the GenomicFeatures package
map2genome <- mapFromTranscripts(gr, transcripts)

# Convert to a data frame
map2genomedf <- as.data.frame(map2genome)

# Bind it back to transcript positions
combined <- cbind(tx_posns, map2genomedf)

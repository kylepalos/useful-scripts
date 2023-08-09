# useful-scripts

A collection of various BASH, R, and Python commands or scripts that do helpful things.

## Script #1: binary_category_table_converter.R"

If you make venn diagrams, Euler diagrams, or upset plots a lot, it is often useful to determine what values belong to which categories.
For example, if I wanted to investigate the 11 genes that are responsive to all stressors except light in the following plot, this script would make that straightforard.

<img width="515" alt="image" src="https://github.com/kylepalos/useful-scripts/assets/56089443/53b769e9-b59c-4a36-947b-65588b559510">


## Script #2: map_transcript_pos_to_genome.R"

Some algorithms (e.g., Nanopore RNA modification tools) report modified sites as numeric positions on transcripts. Downstream analyses such as comparative genomics or metaplot visualization require a genomic coordinate. This script performs that conversion.


## Script #3: quantify_rnaseq_with_salmon.sh"

Algorithms like [salmon](https://combine-lab.github.io/salmon/) and [kallisto](https://pachterlab.github.io/kallisto/about) quantify RNA-seq data without mapping reads to a reference genome. This seems to be the community's preferred way to quantify RNA-seq data due to it's accuracy, efficiency, and reproducibility. This script goes through making the index for Salmon and quantifying RNA-seq reads.


## Script #4: correlation_gene_expression.R"

I often work with adjacent or overlapping genes and want to know if a subset of them have correlated gene expression patterns. This script walks through making a fake example count matrix and correlates all genes in a pairwise manner which can be easily worked with for downstream analyses.

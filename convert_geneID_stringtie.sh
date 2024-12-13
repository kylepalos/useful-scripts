#!/bin/bash

awk '
{
    # Check if the line contains a gene_name entry
    # Intergenic transcripts will not
    if ($0 ~ /gene_name/) {
        # Extract the string that comes after gene_name
        # This is the gene identifier that we want
        match($0, /gene_name "([^"]+)"/, gene_name_match)
        gene_name = gene_name_match[1]
        
        # Replace the gene_id only if we found a gene_name
        gsub(/gene_id "[^"]+"/, "gene_id \"" gene_name "\"")
    }
    
    # Print the line (modified or original)
    print $0
}' input.gtf > output.gtf

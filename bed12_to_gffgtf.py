import sys

# Function to convert a BED12 file to GFF or GTF format
def bed12_to_gff(input_file, output_file, feature_type="exon"):
    with open(input_file, 'r') as in_file, open(output_file, 'w') as out_file:
        # Dictionary to store information about transcripts
        transcripts = {}
        
        for line in in_file:
            fields = line.strip().split('\t')
            chrom, start, end, name, score, strand, thick_start, thick_end, rgb, block_count, block_sizes, block_starts = fields

            # Calculate exon coordinates
            block_sizes = [int(size) for size in block_sizes.split(',')]
            block_starts = [int(start) + int(pos) for pos in block_starts.split(',')]
            exons = [(start, start + size) for start, size in zip(block_starts, block_sizes)]

            # Create a unique transcript identifier
            transcript_id = f'transcript_{name}'

            if transcript_id not in transcripts:
                transcripts[transcript_id] = {
                    'start': int(start),
                    'end': int(end),
                    'exons': []
                }
                
            transcripts[transcript_id]['exons'].extend(exons)
            transcripts[transcript_id]['start'] = min(transcripts[transcript_id]['start'], int(start))
            transcripts[transcript_id]['end'] = max(transcripts[transcript_id]['end'], int(end))

        # Write the GFF/GTF output
        for transcript_id, data in transcripts.items():
            gff_attributes = f'gene_id "{name}"; transcript_id "{transcript_id}";'
            out_file.write(f'{chrom}\tBED2GFF\ttranscript\t{data["start"]+1}\t{data["end"]}\t{score}\t{strand}\t.\t{gff_attributes}\n')
            for i, (exon_start, exon_end) in enumerate(data['exons']):
                gff_attributes = f'gene_id "{name}"; transcript_id "{transcript_id}"; exon_number {i+1};'
                out_file.write(f'{chrom}\tBED2GFF\texon\t{exon_start+1}\t{exon_end}\t{score}\t{strand}\t.\t{gff_attributes}\n')

if len(sys.argv) != 4:
    print("Usage: python bed12_to_gff.py input.bed12 output.gff/gtf feature_type")
    sys.exit(1)

# Get input file path, output file path, and feature type from command-line arguments
input_file = sys.argv[1]
output_file = sys.argv[2]
feature_type = sys.argv[3]

# Call the function to perform the conversion
bed12_to_gff(input_file, output_file, feature_type)

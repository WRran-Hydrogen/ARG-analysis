#!/bin/bash
# This script takes two arguments: the protein sequence file and the protein database file
# It runs meme and fimo to find motifs and filter the database
# It saves the output in meme_out and fimo_out folders

# Check if the arguments are provided
if [ $# -ne 2 ]; then
  echo "Usage: bash script.sh protein.fasta protein_db.fasta"
  exit 1
fi

# Assign the arguments to variables
protein=$1
protein_db=$2

# Run meme to find motifs
meme $protein -protein -mod zoops -nmotifs 5 -minw 6 -maxw 20 -o meme_out

# Run fimo to filter the database
fimo --oc fimo_out --thresh 1e-5 meme_out/meme.txt $protein_db

# Extract the sequence IDs from fimo output
cut -f 3 fimo_out/fimo.tsv | tail -n +2 > fimo_ids.txt

# Use seqkit grep to retrieve the sequences from protein database
seqkit grep -f fimo_ids.txt $protein_db > filtered_protein.fasta

# Print a message when done
echo "Done. Please check meme_out and fimo_out folders for results!"

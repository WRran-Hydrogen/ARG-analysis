import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Merge csv and tsv files.')
parser.add_argument('csv_file', type=str, help='Path to the csv file')
parser.add_argument('tsv_file', type=str, help='Path to the tsv file')
parser.add_argument('output_file', type=str, help='Path to the output csv file')

args = parser.parse_args()

csv_data = pd.read_csv(args.csv_file, header=None)
tsv_data = pd.read_csv(args.tsv_file, sep='\t')
filter_tsv_data = tsv_data[['genome_id', 'gtdb_lineage']]
rmdup_df = filter_tsv_data.drop_duplicates()
merged_data = pd.merge(csv_data, rmdup_df, left_on=1, right_on='genome_id')
merged_data.to_csv(args.output_file, index=False)


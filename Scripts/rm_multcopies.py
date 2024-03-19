import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Process a csv file. Delete multiple copies!')
parser.add_argument('input_file', type=str, help='The input csv file')
parser.add_argument('output_file', type=str, help='The output csv file')
args = parser.parse_args()

df = pd.read_csv(args.input_file)
df = df[df['pfam'].isin(['PF00144', 'PF00905', 'PF13354'])]
duplicates = df[df.duplicated(['gene_id'], keep=False)]['gene_id'].unique()
df = pd.read_csv(args.input_file)
df = df[~df['gene_id'].isin(duplicates)]
df.to_csv(args.output_file, index=False)


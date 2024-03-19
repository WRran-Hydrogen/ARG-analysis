import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Remove rows with duplicate IDs in the first column of a CSV file and no duplicate items.')
parser.add_argument('input_file', type=str, help='The input CSV file')
parser.add_argument('output_file', type=str, help='The output CSV file')
args = parser.parse_args()

df = pd.read_csv(args.input_file)
df = df.drop_duplicates(subset=df.columns[0], keep=False)
df.to_csv(args.output_file, index=False)


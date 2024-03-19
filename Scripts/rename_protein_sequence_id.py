import argparse
import re

def rename_protein_sequence_id(input_file, output_file):
    with open(input_file, 'r') as f:
        content = f.read()
    pattern = r'^(>.*?)\|.*?\s.*?#\s(\d+)\s#.*$'
    new_content = re.sub(pattern, r'\1 \2', content, flags=re.MULTILINE)
    with open(output_file, 'w') as f:
        f.write(new_content)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Rename protein sequence IDs.')
    parser.add_argument('input_file', type=str, help='Path to the input file.')
    parser.add_argument('output_file', type=str, help='Path to the output file.')
    args = parser.parse_args()
    rename_protein_sequence_id(args.input_file, args.output_file)


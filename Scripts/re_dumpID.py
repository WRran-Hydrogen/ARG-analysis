import argparse

parser = argparse.ArgumentParser(description='Process input and output files. Add extra identifiers in same ID to distinguish them')
parser.add_argument('input_file', type=str, help='Input file name')
parser.add_argument('output_file', type=str, help='Output file name')

args = parser.parse_args()

with open(args.input_file, 'r') as f:
    lines = f.readlines()

ids = {}
output = []

for line in lines:
    if line.startswith('>'):
        id = line.split()[0]
        if id in ids:
            id_parts = id.split('|')
            id_parts[0] += '.' + str(ids[id])
            new_id = '|'.join(id_parts)
            ids[id] += 1
            output.append(new_id + '\n')
        else:
            ids[id] = 2
            output.append(id + '\n')
    else:
        output.append(line)

with open(args.output_file, 'w') as f:
    f.writelines(output)


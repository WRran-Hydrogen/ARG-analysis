import argparse
import pandas as pd

# 创建命令行参数解析器
parser = argparse.ArgumentParser(description='处理csv文件')
parser.add_argument('csv_file', type=str, help='csv文件的名称')
parser.add_argument('output_file', type=str, help='输出文件的名称')

# 解析命令行参数
args = parser.parse_args()

# 读取csv文件
df = pd.read_csv(args.csv_file)

# 筛选出同时存在PF00239和PF00589的seq_id
result1 = df[(df['pfam'] == 'PF00239') | (df['pfam'] == 'PF00589')].groupby('seq_id').filter(lambda x: x['pfam'].nunique() == 2)['seq_id'].unique()

# 输出结果
print(f'在同一个seq_id下同时存在PF00239和PF00589的seq_id有{len(result1)}个')

# 筛选出存在PF00239或PF00589的seq_id
result2 = df[(df['pfam'] == 'PF00239') | (df['pfam'] == 'PF00589')]['seq_id'].unique()

# 输出结果
print(f'存在PF00239或PF00589的seq_id有{len(result2)}个')

# 筛选出存在Resolvase的seq_id
result3 = df[df['feat_id'] == 'Resolvase']['seq_id'].unique()

# 输出结果
print(f'存在Resolvase的seq_id有{len(result3)}个')

# 筛选出存在Phage_integrase的seq_id
result4 = df[df['feat_id'] == 'Phage_integrase']['seq_id'].unique()

# 输出结果
print(f'存在Phage_integrase的seq_id有{len(result4)}个')

# 将result1列表写入文件
with open(args.output_file, 'w') as f:
    for item in result3:
        f.write(f'{item}\n')
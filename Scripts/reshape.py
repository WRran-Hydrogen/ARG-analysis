import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Convert wide format csv to long format.')
parser.add_argument('--input', type=str, required=True, help='Input csv file')
parser.add_argument('--output', type=str, required=True, help='Output csv file')
args = parser.parse_args()

# 读取 csv 文件
df = pd.read_csv(args.input)

# 获取行名列名
row_name = df.columns[0]
col_names = df.columns[1:]

# 将宽表格转换为长表格
long_df = pd.melt(df, id_vars=row_name, value_vars=col_names, var_name='Column', value_name='Value')

# 显示结果
long_df.to_csv(args.output, index=False)


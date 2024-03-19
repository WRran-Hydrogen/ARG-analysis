import pandas as pd
import argparse

def process_csv(input_csv, output_csv):
    # 读取CSV文件
    df = pd.read_csv(input_csv)

    # 对每个源进行分组，并按权重降序排序
    df = df.sort_values(['source', 'weight'], ascending=[True, False])

    # 计算每个源的最大权重，并将其缩减到300
    max_weights = df.groupby('source')['weight'].transform(max)
    reduction_factor = 300 / max_weights

    # 应用缩减因子
    df['weight'] = df['weight'] * reduction_factor

    # 输出CSV文件
    df.to_csv(output_csv, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process a CSV file.')
    parser.add_argument('--input', type=str, help='Input CSV file')
    parser.add_argument('--output', type=str, help='Output CSV file')

    args = parser.parse_args()

    process_csv(args.input, args.output)


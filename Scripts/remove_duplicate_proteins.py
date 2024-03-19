import argparse
from Bio import SeqIO
from collections import Counter

def remove_duplicate_proteins(input_file, output_file):
    # 读取输入文件中的序列
    records = list(SeqIO.parse(input_file, "fasta"))
    # 创建一个字典来存储序列的计数
    record_counts = Counter(record.id.split("|")[0] for record in records)
    # 筛选出唯一ID的序列
    unique_records = [record for record in records if record_counts[record.id.split("|")[0]] == 1]
    # 将唯一的序列写入输出文件中
    SeqIO.write(unique_records, output_file, "fasta")

if __name__ == "__main__":
    # 创建一个解析器对象
    parser = argparse.ArgumentParser(description="Remove duplicate proteins from a faa file.")
    # 添加输入文件和输出文件参数
    parser.add_argument("input_file", help="The input faa file.")
    parser.add_argument("output_file", help="The output faa file.")
    # 解析命令行参数
    args = parser.parse_args()
    # 调用函数来删除冗余序列
    remove_duplicate_proteins(args.input_file, args.output_file)

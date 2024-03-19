import json
import csv

# 假设你的JSON数据存储在一个名为"data.json"的文件中
with open("C:/Users/67477/Desktop/out.json", 'r') as f:
    data = json.load(f)

# 打开一个新的CSV文件并创建一个写入器
with open('output_1.csv', 'w', newline='') as csvfile:
    fieldnames = ['key', 'sub_key', 'Baseline MG94xREV', 'Baseline MG94xREV omega ratio', 'Corrected P-value', 'Full adaptive model', 'Full adaptive model (non-synonymous subs/site)', 'Full adaptive model (synonymous subs/site)', 'LRT', 'Nucleotide GTR', 'Rate Distributions', 'Rate classes', 'Uncorrected P-value', 'original name']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    # 写入表头
    writer.writeheader()

    # 提取"branch attributes"中的"0"键下的所有子键
    sub_keys = data['branch attributes']['0'].keys()

    # 遍历所有子键
    for sub_key in sub_keys:
        # 提取并写入每个子键及其对应的值
        writer.writerow({'key': '0', 'sub_key': sub_key, **data['branch attributes']['0'][sub_key]})
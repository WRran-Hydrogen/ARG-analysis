import pandas as pd

# 读取表格数据
df = pd.read_csv('/mnt/c/Users/67477/OneDrive/ARGs_genes/new/HGT/no_HGT/Distance_noHGT.csv')

# 创建结果表格
result_df = pd.DataFrame(columns=['Gene_id', 'Feat_id', 'Name', 'Other pfam feat_id', 'Other pfam name'])
# 对每一个gene_id进行处理
for gene_id in df['gene_id'].unique():
    gene_df = df[df['gene_id'] == gene_id]
    
    # 筛选出PF00905或PF13354
    pfam_filter = gene_df['pfam'].isin(['PF00905', 'PF13354', 'PF00144'])
    pfam_df = gene_df[pfam_filter]

    # 计算其他pfam到PF00905或PF13354的距离
    def calculate_distance(row):
        distances = []
        for index, pfam_row in pfam_df.iterrows():
            start_diff = abs(row['start'] - pfam_row['start'])
            end_diff = abs(row['end'] - pfam_row['end'])
            distance = (start_diff + end_diff) / 2 / row['length']
            distances.append(distance)
        return distances

    gene_df['distance'] = gene_df.apply(calculate_distance, axis=1)

    # 添加结果到结果表格中
    for index, row in pfam_df.iterrows():
        temp_df = gene_df[~pfam_filter].copy()
        temp_df['Gene_id'] = gene_id
        temp_df['Feat_id'] = row['feat_id']
        temp_df['Name'] = row['name']
        temp_df.rename(columns={'feat_id': 'Other pfam feat_id', 'name': 'Other pfam name'}, inplace=True)
        result_df = pd.concat([result_df, temp_df[['Gene_id', 'Feat_id', 'Name', 'Other pfam feat_id', 'Other pfam name', 'distance']]], ignore_index=True)

# 显示结果
result_df.to_csv('/mnt/c/Users/67477/OneDrive/ARGs_genes/new/HGT/no_HGT/Distance_noHGT_result.csv',index=False)


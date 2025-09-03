import pandas as pd
import sys

# 检查命令行参数
if len(sys.argv) != 4:
    print("Usage: python find_incompatible_snps.py <bim_file1> <bim_file2> <bim_file3>")
    sys.exit(1)

# 读取三个bim文件
try:
    bim1 = pd.read_csv(sys.argv[1], sep='\t', header=None, names=['CHR', 'SNP', 'CM', 'POS', 'A1', 'A2'])
    bim2 = pd.read_csv(sys.argv[2], sep='\t', header=None, names=['CHR', 'SNP', 'CM', 'POS', 'A1', 'A2'])
    bim3 = pd.read_csv(sys.argv[3], sep='\t', header=None, names=['CHR', 'SNP', 'CM', 'POS', 'A1', 'A2'])
except FileNotFoundError as e:
    print(f"Error: {e}")
    sys.exit(1)

# 创建一个唯一的CHR:POS标识符
bim1['pos_id'] = bim1['CHR'].astype(str) + ':' + bim1['POS'].astype(str)
bim2['pos_id'] = bim2['CHR'].astype(str) + ':' + bim2['POS'].astype(str)
bim3['pos_id'] = bim3['CHR'].astype(str) + ':' + bim3['POS'].astype(str)

# 创建一个字典来存储每个位置的等位基因信息
pos_alleles = {}

# 填充字典
for _, row in pd.concat([bim1, bim2, bim3]).iterrows():
    pos_id = row['pos_id']
    alleles = frozenset([row['A1'], row['A2']]) # 使用frozenset来忽略顺序
    if pos_id not in pos_alleles:
        pos_alleles[pos_id] = set()
    pos_alleles[pos_id].add(alleles)

# 找出存在多种不同等位基因组合的位置
incompatible_pos = {pos for pos, allele_sets in pos_alleles.items() if len(allele_sets) > 1}

# 输出不兼容的位置ID，用于后续PLINK排除
for pos in incompatible_pos:
    print(pos)
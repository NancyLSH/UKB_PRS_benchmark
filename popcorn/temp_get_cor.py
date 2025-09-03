import os
import pandas as pd

def parse_popcorn_output(file_path):
    """解析单个Popcorn输出文件并提取所需数据。"""
    results = []
    trait_name = os.path.basename(file_path).replace('.txt', '')
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
    # 跳过表头
    for line in lines[1:]:
        parts = line.split()
        if not parts:
            continue
        
        parameter = parts[0]
        if parameter in ["h1^2", "h2^2", "pge"]:
            value = float(parts[1])
            se = float(parts[2])
            results.append({
                "Trait": trait_name,
                "Parameter": parameter,
                "Value": value,
                "SE": se
            })
            
    return results

def main():
    """主函数，遍历目录，处理文件并保存为TSV。"""
    # *** 请修改为您的Popcorn结果文件所在的目录 ***
    input_directory = '.' 
    
    all_results = []
    
    print(f"开始扫描目录: {os.path.abspath(input_directory)}")
    
    for filename in os.listdir(input_directory):
        if filename.endswith(".txt"):
            file_path = os.path.join(input_directory, filename)
            print(f"正在处理文件: {filename}")
            try:
                file_results = parse_popcorn_output(file_path)
                all_results.extend(file_results)
            except Exception as e:
                print(f"处理文件 {filename} 时出错: {e}")

    if not all_results:
        print("警告: 未找到任何有效的Popcorn结果文件或未能提取任何数据。")
        return

    # 将结果转换为DataFrame并保存为TSV文件
    df = pd.DataFrame(all_results)
    output_path = "popcorn_results.tsv"
    df.to_csv(output_path, sep='\t', index=False)
    
    print(f"\n处理完成！数据已成功保存到: {output_path}")
    print("文件内容预览:")
    print(df)

if __name__ == "__main__":
    main()
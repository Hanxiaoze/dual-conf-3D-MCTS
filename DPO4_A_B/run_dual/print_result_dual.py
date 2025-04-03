import os
import shutil
from openbabel import pybel

file_path_affi_dict = {}
def extract_score_line(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for i in range(len(lines)):
            if '>  <Score>  (1)' in lines[i]:
                if i + 1 < len(lines):
                    score = float(lines[i + 1].strip())
                    # if -1000.0 < score <= -7.0:
                    #     print(f"File: {file_path}, Score: {score}")      
            if '>  <Score_dual>  (1)' in lines[i]:                        
                if i + 1 < len(lines):
                    score_dual = float(lines[i + 1].strip())
                    # if -1000.0 < score_dual <= -7.0:
                    #     print(f"File: {file_path}, Dual_Score: {score_dual}")    
                    file_path_affi_dict[file_path] = str(score) + '     ' + str(score_dual)              
                break

import re
# 定义一个函数用于提取文件名中的数字部分
def extract_number(filename):
    match = re.search(r'record_(\d+)\.sdf', filename)
    if match:
        return int(match.group(1))
    return 0  # 如果没有匹配到数字，返回 0

#  大函数
def process_sdf_files_in_directory(directory):
    files = [f for f in os.listdir(directory) if f.startswith('record_') and f.endswith('.sdf')]
    sorted_files = sorted(files, key=extract_number)
    for filename in sorted_files:        
        file_path = os.path.join(directory, filename)
        extract_score_line(file_path)


def process_file(input_file, v):            
    output_file = os.path.join(output_folder, os.path.basename(input_file).replace('.sdf', f'.smi'))
    # 逐个读取mol2分子并将其转换为xyz格式
    for mol in pybel.readfile("sdf", input_file):
            #mol.addh()  # 添加氢原子
            #mol.make3D()  # 生成3D结构            
            output_stream = pybel.Outputfile('smi', output_file, overwrite=True)  # 创建输出流        
            output_stream.write(mol)  # 将分子写入 smi 文件
            output_stream.close()  # 关闭输出文件
            smiles = mol.write("smi").strip()  # 将分子转换为SMILES字符串并去除多余的空格
            print(f"Generated mol {os.path.splitext(input_file)[0].split('/')[-1]} binding affinity is {v} kcal/mol, SMILES is:  ", smiles)



if __name__ == "__main__":
    work_directory = './record'
    output_folder = './record_smi'
    if os.path.exists(output_folder):        
        shutil.rmtree(output_folder)        
    os.makedirs(output_folder)    
    process_sdf_files_in_directory(work_directory)
    for f, v in file_path_affi_dict.items():
        # process_file(f, v)
        print(f, '   ', v)
    

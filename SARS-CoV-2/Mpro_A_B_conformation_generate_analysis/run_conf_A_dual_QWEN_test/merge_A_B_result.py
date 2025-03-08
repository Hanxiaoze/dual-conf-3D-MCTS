import os
import shutil
from openbabel import pybel

file_path_affi_dict = {}
def extract_gensdf_score_line(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for i in range(len(lines)):
            if '>  <Score>  (1)' in lines[i]:
                if i + 1 < len(lines):
                    score = float(lines[i + 1].strip())
                    # if -1000.0 < score <= -7.0:
                    #     print(f"File: {file_path}, Score: {score}")    
                    file_path_affi_dict[file_path] = str(score)              
                break

def extract_dock_score_line(file_path, txt):
    sss = file_path.split('/')[-1].split('.')[0]
    sss = sss + '.pdbqt'
    with open(txt) as fr:
        for li in fr.readlines():    #['./validate_dock_another_target/record_2.pdbqt', -8.0]
            if sss in li:
                sco = li.split(',')[-1][:-2]
                file_path_affi_dict[file_path] = file_path_affi_dict[file_path] + '   ' + str(sco)
                break
    return

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
        extract_gensdf_score_line(file_path)
        extract_dock_score_line(file_path, '../gen_A_dock_B.txt')




if __name__ == "__main__":
    work_directory = './record'   
        
    process_sdf_files_in_directory(work_directory)
    for f, v in file_path_affi_dict.items():
        
        print(f, '   ', v)
    

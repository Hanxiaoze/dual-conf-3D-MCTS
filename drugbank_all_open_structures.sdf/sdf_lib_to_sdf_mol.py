import os
from openbabel import pybel

# 输入文件夹和输出文件夹
input_file = './open_structures.sdf'

# 获取输入文件夹中的所有  文件

i = 0
for m in pybel.readfile("sdf", input_file):   # 读取此  文件
    i += 1
    output_file = rf'./mol_sdf_3D/drugbank_mol_{i}.sdf'
    m.addh()  # 添加氢原子
    m.make3D()  # 生成3D结构
    
    output_stream = pybel.Outputfile('sdf', output_file, overwrite=True)  # 创建输出流
    output_stream.write(m)  # 将分子写入 pdb 文件
    output_stream.close()  # 关闭输出文件
    print(f"{output_file} sdf file with H generated successfully!")






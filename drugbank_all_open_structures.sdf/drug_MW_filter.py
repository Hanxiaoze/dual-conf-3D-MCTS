import os
import shutil
from rdkit import Chem
from rdkit.Chem import Descriptors

def filter_sdf_files(input_folder, output_folder, mol_weight_threshold=600):
    # 如果输出文件夹不存在，则创建
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # 遍历输入文件夹下的所有 .sdf 文件
    for sdf_file in os.listdir(input_folder):
        if sdf_file.endswith('.sdf'):
            input_file_path = os.path.join(input_folder, sdf_file)
            output_file_path = os.path.join(output_folder, sdf_file)

            # 读取每个 .sdf 文件
            suppl = Chem.SDMolSupplier(input_file_path)            

            for mol in suppl:
                if mol is None:
                    continue

                # 计算分子量并筛选
                mol_weight = Descriptors.MolWt(mol)
                if mol_weight <= mol_weight_threshold:
                    writer = Chem.SDWriter(output_file_path)
                    print(rf'{sdf_file} M.W. =<<<<< 600, save file!')
                    writer.write(mol)
                    writer.close()
                else:
                    print(rf'{sdf_file} M.W. >>>> 600')

            

# 示例用法
input_folder = './mol_sdf_2D'
output_folder = './mol_sdf_2D_600'
filter_sdf_files(input_folder, output_folder)

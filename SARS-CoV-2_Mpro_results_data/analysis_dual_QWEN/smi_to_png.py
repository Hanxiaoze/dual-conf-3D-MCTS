import os
from rdkit import Chem
from rdkit.Chem import Draw

def draw_molecules_from_smiles_files(folder_path, png_path):
    # 检查文件夹是否存在
    if not os.path.exists(folder_path):
        print(f"文件夹 {folder_path} 不存在")
        return

    # 获取文件夹中所有的.smi文件
    smi_files = [f for f in os.listdir(folder_path) if f.endswith('.smi')]

    for smi_file in smi_files:
        file_path = os.path.join(folder_path, smi_file)
        with open(file_path, 'r') as file:
            for line in file:
                smiles = line.strip()
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    # 生成分子图像
                    img = Draw.MolToImage(mol, size=(500, 500))
                    # 保存图像文件
                    img_file = os.path.join(png_path, f"{os.path.splitext(smi_file)[0]}.png")
                    img.save(img_file)
                    print(f"保存图像: {img_file}")
                else:
                    print(f"无效的SMILES: {smiles}")

# 使用示例
folder_path = './filtered_mol'  # 将XXX替换为你的文件夹路径
png_path = './smi_to_png'
draw_molecules_from_smiles_files(folder_path, png_path)

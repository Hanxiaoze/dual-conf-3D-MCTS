import random
import math
import hashlib
import argparse
from rdkit import Chem as ch
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
from rdkit.Chem import rdMolTransforms
import numpy as np
import os
from openbabel.pybel import *
import time
import multiprocessing, traceback
from scipy.spatial.distance import cdist
from func_timeout import func_set_timeout
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import QED
from rdkit.Chem import Descriptors

import shutil

def remove_make_dir(folder_path):
    # 检查文件夹是否存在
    if os.path.exists(folder_path) and os.path.isdir(folder_path):
        # 删除文件夹及其内容
        shutil.rmtree(folder_path)
        print(f"Folder '{folder_path}' has been deleted, making...")
        os.system(rf'mkdir {folder_path}')
    else:
        print(f"Folder '{folder_path}' does not exist, making...")
        os.system(rf'mkdir {folder_path}')

folder_path = './frags_3d_sdf'

remove_make_dir(folder_path)

FRAGMENT_LIB = []
with open(f'./my_frags_2.txt','r') as f:
    lines = f.readlines()
    for line in lines:
        FRAGMENT_LIB.append(line.strip())

i = 0
for smi in FRAGMENT_LIB:
    # 创建分子对象并添加氢原子
    mol = Chem.MolFromSmiles(smi)
    mol_h = Chem.AddHs(mol)

    try:     
        for a in mol_h.GetAtoms():
            if a.GetIsotope() != 0:            
                a.SetIsotope(0)
                a.SetAtomicNum(1)
        
        # 生成 3D 构象
        success = AllChem.EmbedMolecule(mol_h, randomSeed=42)  # randomSeed 设置一个种子，确保结果可复现
        # 优化 3D 构象（UFF 优化方法）
        AllChem.UFFOptimizeMolecule(mol_h)
    except:
        print(rf'frag ID: {i}:  {smi}  failed in 3D conform!')


    # 写入 SDF 文件
    writer = Chem.SDWriter(rf'{folder_path}/frag_{i}.sdf')
    writer.write(mol_h)
    writer.close()

    i += 1
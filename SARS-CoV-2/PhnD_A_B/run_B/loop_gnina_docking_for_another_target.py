import os
import re
from openbabel.pybel import *

GNINA = '/home/zqc-gpu-03/dual-conf-3D-MCTS'
ADFR = '/home/zqc-gpu-03/dual-conf-3D-MCTS/ADFRsuite/ADFRsuite-1.0/bin'


lig_dir = './record_100'
pro_dual_pdb = '../PhnD_A_c9.pdb'
ligand_box_dual = '../ligand_A.sdf'

pro_dual_pdbqt = '../PhnD_A_c9.pdbqt'

os.system(rf'{ADFR}/prepare_receptor -r {pro_dual_pdb} -o {pro_dual_pdbqt}')


if not os.path.exists('./validate_dock_another_target'):
    os.system('mkdir validate_dock_another_target')

if not os.path.exists('./validate_dock_another_target_bindpose'):
    os.system('mkdir validate_dock_another_target_bindpose')


lig_sdf_list = [lig_dir+'/'+f for f in os.listdir(lig_dir) if f.endswith('.sdf')]


def ligand_sdf_to_pdbqt(lig_sdf):

    mol = list(readfile('sdf', rf'{lig_sdf}'))[0]
    lig_pdbqt = './validate_dock_another_target/' + str(lig_sdf).split('.')[1].split('/')[-1] + '.pdbqt'
    mol.write("pdbqt", lig_pdbqt, overwrite=True)
    # print(lig_pdbqt)

    return lig_pdbqt



def dual_dock_score(lig_pdbqt, score=None, smile=None):   # 使用 GNINA 对接    

    lig_name = str(lig_pdbqt).split('.')[1].split('/')[-1]
    try:
        score_dual = float(os.popen(                                                    # 用来定位蛋白口袋位置  --autobox_add 0
            rf'{GNINA}/gnina --receptor {pro_dual_pdbqt} --ligand {lig_pdbqt} --autobox_ligand {ligand_box_dual} --cnn_scoring=none --out ./validate_dock_another_target_bindpose/{lig_name}_dual_dock.sdf').read().split(
            '-----+------------+------------+----------')[1].split()[1])
        if score_dual > -7.5:
            os.system(f'rm ./validate_dock_another_target_bindpose/{lig_name}_dual_dock.sdf')

    except:
        score_dual = 0
 
    return [lig_pdbqt, score_dual]


# 自定义排序键函数，从文件名中提取record_后的数字部分
def extract_record_number(file_path):
    match = re.search(r'record_(\d+)', file_path.split('/')[-1])
    return int(match.group(1)) if match else float('inf')

# 根据提取的数字部分排序文件列表
sorted_lig_sdf_list = sorted(lig_sdf_list, key=extract_record_number)


if __name__ == "__main__":
    result_list = []
    print(sorted_lig_sdf_list)
    for lig_sdf in sorted_lig_sdf_list[: ]:
        # print(lig_sdf)
        lig_pdbqt = ligand_sdf_to_pdbqt(lig_sdf)
        iterm = dual_dock_score(lig_pdbqt)
        print(iterm)
        result_list.append(iterm)
        


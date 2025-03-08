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


import warnings
from rdkit import RDLogger
# 屏蔽 RDKit 的警告信息
RDLogger.DisableLog('*MDLV30*')
warnings.filterwarnings("ignore", category=UserWarning, module="openbabel")
warnings.filterwarnings("ignore")

from contextlib import redirect_stdout


# Please replace the PATH for GNINA and ADFR first.
GNINA = '/home/zhou/3D-MCTS-main'
ADFR = '/home/zhou/3D-MCTS-main/ADFRsuite/ADFRsuite-1.0/bin'

parser = argparse.ArgumentParser(description='3D-MCTS code, for molecular generation')
parser.add_argument('--num_sims', action="store", type=int, default=1000000, help='Number of simulation steps')
parser.add_argument('--ligand', action="store", type=str, help='sdf file to determine the position of pocket',   # 用来定位蛋白口袋位置
                    default='ligand.sdf')
parser.add_argument('--protein', action="store", type=str, help='protein, PDB format',
                    default='protein.pdb')
parser.add_argument('--pocket', action="store", type=str, help='pocket, PDB format',
                    default='pocket.pdb')
parser.add_argument('--score', action="store", type=float, help='threshold for vina score', default=-7)
parser.add_argument('--qed', action="store", type=float, help='threshold for qed', default=0.3)
parser.add_argument('--processor', action="store", type=int, help='number of processor for multiprocessing', default=48)
parser.add_argument('--start', action="store", type=str, help='start fragment', default='1')
parser.add_argument('--frag_lib', action="store", type=str, help='fragment library', default='frags/fragment.txt')
parser.add_argument('--gnina', action="store", type=str, help='the path for GNINA',
                    default= GNINA)
parser.add_argument('--adfr', action="store", type=str, help='the path for adfr',
                    default= ADFR)

args = parser.parse_args()

os.system('mkdir node vinaa record tmp init_pose unique')
os.system(rf'{ADFR}/prepare_receptor -r {args.protein} -o pro_test.pdbqt')

def print_node_info(node):
        
    try:
        print('\n'+rf'node.id: {node.id}')            
    except:
        print('\nnode.id: None') 

    try:
        print(rf'node.state.type: 【{node.state.type}】')            
    except:
        print('node.state.type: 【None】') 

    try:
        print(rf'node.parent.id: {node.parent.id}')            
    except:
        print('node.parent.id: None') 

    try:
        print(rf'node.parent: {node.parent}')            
    except:
        print('node.parent: None') 

    try:
        print(rf'node.visits: {node.visits}')            
    except:
        print('node.visits: None')

    try:
        print(rf'node.reward: {node.reward}')            
    except:
        print('node.reward: None') 

    try:
        print(rf'node.best_score: {node.best_score}')            
    except:
        print('node.best_score: None') 

    try:
        print(rf'node.longest_path: {node.longest_path}')            
    except:
        print('node.longest_path: None') 

    try:
        print(rf'node.qed: {node.qed}')            
    except:
        print('node.qed: None') 

    try:
        print(rf'node.state.Frag_Deg: {node.state.Frag_Deg}')            
    except:
        print('node.state.Frag_Deg: None') 

    try:
        print(rf'node.state.score: {node.state.score}')            
    except:
        print('node.state.score: None') 

    try:
        print(rf'node.state.sdf: {node.state.sdf}')            
    except:
        print('node.state.sdf: None') 

    try:
        print(rf'node.state.h1s_avail: {node.state.h1s_avail}')            
    except:
        print('node.state.h1s_avail: None') 

    try:
        print(rf'node.state.h1: {node.state.h1}')            
    except:
        print('node.state.h1: None') 

    try:
        print(rf'node.state.states: {node.state.states}')            
    except:
        print('node.state.states: None') 

    try:
        print(rf'node.state.choices: {node.state.choices}')            
    except:
        print('node.state.choices: None') 

    try:
        print(rf'node.state.ter_tag: {node.state.ter_tag}')            
    except:
        print('node.state.ter_tag: None') 

    try:
        c = 0
        if len(node.children) == 0:
            print('node.children: None')

        for child in node.children:
            print(rf'node.children_{c}_id: {child.id}, state: {child}')
            c += 1
    except:
        print('node.children: None')

    


# Used to balance exploration and exploitation
SCALAR = 1 / (2 * math.sqrt(2.0))
FRAGMENT_LIB = []
with open(f'{args.frag_lib}','r') as f:
    lines = f.readlines()
    for line in lines:
        FRAGMENT_LIB.append(line.strip())

GNINA = args.gnina
ADFR = args.adfr

# Rules for reassembly of fragments    这是一个对称矩阵，存在 i : j，就存在 j : i
Iso_dic = {
    1: [3, 5, 10],  ###
    2: [],
    3: [1, 6, 4, 8, 13, 14, 15, 16],      # 只有一种 O
    4: [3, 5, 11, ],  ####
    5: [1, 6, 4, 8, 12, 14, 16, 13, 15],
    6: [13, 14, 15, 16, 3, 5, 10],
    7: [],
    8: [3, 5, 11, 9, 10, 13, 14, 15, 16],
    9: [8, 13, 14, 15, 16],
    10: [1, 6, 8, 13, 14, 15, 16],
    11: [4, 8, 13, 14, 15, 16],
    12: [5],
    13: [3, 5, 6, 8, 9, 10, 11, 14, 15, 16],
    14: [3, 5, 6, 8, 9, 10, 11, 13, 14, 15, 16],
    15: [3, 5, 6, 8, 9, 10, 11, 13, 14, 16],
    16: [3, 5, 6, 8, 9, 10, 11, 13, 14, 15, 16]
}
pro = 'pro_test.pdbqt'
# init = ''
NODE_ID = 1
RECORD = 1
SCORES = []
GOOD_SCORES = []
NO_QED = []
NO_QED_GOOD = []
mol_dic = []
mol = Chem.MolFromPDBFile(rf"{args.pocket}")
writer = Chem.SDWriter(rf'./pocket.sdf')
writer.write(mol)
writer.close()
POCKET = Chem.SDMolSupplier('pocket.sdf')[0]
P = POCKET.GetConformer()
POCKET_CORS = [P.GetPositions()[i] for i in range(len(POCKET.GetAtoms()))]



def dock(init):      # 只用来产生第一步的结合 pose

    '''
    Generate the conformation of the starting fragment.    # 只用来产生第一步的结合 pose
    '''
    mymol = list(readfile('sdf', rf'../../../init/{init}.sdf'))[0]      # pybel 下面的函数
    mymol.write("pdbqt", rf'../../../init/{init}.pdbqt', overwrite=True)
    os.system(                                                         # 用来定位蛋白口袋位置
        rf'{GNINA}/gnina --receptor {pro} --ligand ../../../init/{init}.pdbqt --autobox_ligand {args.ligand} --cnn_scoring=none --out ./init_pose/{init}.pdbqt')

    modes = open(rf'./init_pose/{init}.pdbqt', 'r').read().split('ENDMDL\n')   # 打开 gnina 对接的输出 pose 文件
    valids = []
    for i in range(9):

        mode = open(rf'./init_pose/{init}_{i}.pdbqt', 'w')      #  分别将前 9 个 pose 写成 .pdbqt文件
        score = float(modes[i].split('\n')[1].split()[2].replace('REMARK', ''))   # 记录第 i 个 pose 的对接得分
        mode.write(modes[i] + 'ENDMDL\n')
        mode.close()

        os.system(rf'obabel ./init_pose/{init}_{i}.pdbqt -O ./init_pose/{init}_{i}.sdf -D -B -R -C')    #  前 9 个 pose 转成 .sdf

        f = open(rf'./init_pose/{init}_{i}.sdf', 'r').readlines()
        f_ = open(rf'./init_pose/{init}_{i}_repair.sdf', 'w')
        for line in f[:4]:
            f_.write(line)
        for j in range(4, len(f)):
            if f[j].startswith(' ') and len(f[j]) >= 50:
                new_line = f[j][:50] + '0' + f[j][51:]     # 修复 .sdf 文件
                f_.write(new_line)
            else:
                br = j
                break
        for j in range(br, len(f)):
            f_.write(f[j])
        f_.close()

        try:
            mol_h = Chem.AddHs(Chem.SDMolSupplier(rf'./init_pose/{init}_{i}_repair.sdf')[0], addCoords=True)  # 加 H
            writer = Chem.SDWriter(rf'./init_pose/{init}_{i}_repair_H.sdf')
            writer.write(mol_h)
            writer.close()
            valids.append([init, i, score])
        except:
            traceback.print_exc()
            continue
    return valids    # [[init, i, score], ...]   


def assign_iso(frag_id, Del_plus, start=False):
    '''
    Isotope labeling of fragment attachment sites      # 同位素标记的是原子类型
    '''
    if start == False:
        lig_name = rf'{frag_id}_{Del_plus}'
        mol_h = Chem.AddHs(Chem.SDMolSupplier(rf'./vinaa/{lig_name}_repair.sdf')[0], addCoords=True)  # vinna生成的
    else:
        lig_name = rf'{frag_id}_{Del_plus}'
        mol_h = Chem.SDMolSupplier(rf'./init_pose/{lig_name}_repair_H.sdf', removeHs=False)[0]   # 初始
    for at in mol_h.GetAtoms():
        if True:
            label = 0    #  0 代表不可连接新片段
            if at.GetAtomicNum() == 6:               # C
                if at.IsInRing():
                    if at.GetIsAromatic():
                        label = 16                   # C 环芳香连接点
                        neis = at.GetNeighbors()
                        for nei in neis:
                            if nei.GetAtomicNum() == 7 or 8 or 16:   # C 环芳香连接点，邻居有 N，O，S
                                label = 14
                    else:                           # C 环非芳香
                        label = 15
                        neis = at.GetNeighbors()
                        for nei in neis:
                            if nei.GetAtomicNum() == 7 or 8 or 16:
                                label = 13
                else:                               # 非环 C
                    neis = at.GetNeighbors()
                    label = 8                 # 根据后面代码，排除其它C，应该是单键 C
                    for nei in neis:
                        bond = mol_h.GetBondBetweenAtoms(at.GetIdx(), nei.GetIdx())
                        if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and nei.GetAtomicNum() == 8:  # C=O 上的 C
                                label = 6
                                break
                            else:             # 普通双键 C， 0 代表不可连接新片段
                                label = 0     #  0 代表不可连接新片段
                                break
            if at.GetAtomicNum() == 7:       # N
                label = 0
                if at.GetIsAromatic():
                    label = 9                 # 芳香 N
                else:
                    if at.IsInRing():        # 非芳香，是环上 N
                        neis = at.GetNeighbors()
                        for nei in neis:
                            if nei.GetAtomicNum() == 6:
                                neis2 = nei.GetNeighbors()
                                for nei2 in neis2:
                                    if nei2.GetAtomicNum() == 8 and mol_h.GetBondBetweenAtoms(nei.GetIdx(),
                                                                                              nei2.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE:
                                        label = 10      # -N-C=O  上 N
                if label != 10:
                    if at.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                        label = 5   # 普通 sp3  N
            if at.GetAtomicNum() == 8:
                label = 3               # 只有一种 O
            if at.GetAtomicNum() == 16:   # S
                if at.GetHybridization() == Chem.rdchem.HybridizationType.SP3:
                    label = 11  #  单键 S
                else:
                    neis = at.GetNeighbors()
                    count = 0
                    for nei in neis:
                        bond = mol_h.GetBondBetweenAtoms(at.GetIdx(), nei.GetIdx())
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and nei.GetAtomicNum() == 8:  # S=O 键
                            count += 1
                    if count == 2:
                        label = 12    # O=S=O
            if label != 0:
                neis = at.GetNeighbors()
                for at in neis:
                    if at.GetAtomicNum() == 1:
                        at.SetIsotope(label)    # 设置同位素标记
    if start == False:  #不是最开始
        writer = Chem.SDWriter(rf'./tmp/{lig_name}.sdf')
        writer.write(mol_h)    # 写下同位素标记完的分子
        writer.close()
        os.system(rf"grep -Ev 'MDLV30|COLLECTION' ./tmp/{lig_name}.sdf > ./tmp/{lig_name}_tt.sdf")
        os.system(rf'cp ./tmp/{lig_name}_tt.sdf ./tmp/{lig_name}.sdf')
    else:
        writer = Chem.SDWriter(rf'./start.sdf')   # 起始
        writer.write(mol_h)
        writer.close()


def get_id_bysymbol(combo, symbol):   # 输入符号，给出 id
    for at in combo.GetAtoms():
        if at.GetSymbol() == symbol:
            return at.GetIdx()


def get_neiid_bysymbol(combo, symbol):   # 输入 at 的符号，给出 at 的 nei[0] 的 id
    for at in combo.GetAtoms():
        if at.GetSymbol() == symbol:
            at_nei = at.GetNeighbors()[0]
            return at_nei.GetIdx()


def get_neiid_byisotope(combo, symbol):
    for at in combo.GetAtoms():
        if at.GetIsotope() != 0:           # 获得片段中所有 “标记同位素的原子” 的 nei[0] 的 id
            at_nei = at.GetNeighbors()[0]
            return at_nei.GetIdx()


def frag_avail(h1):  

    lib_avail = []
    iso_avail = Iso_dic[h1]    # 查连接规则字典 dic
    for i in iso_avail:
        for smi in FRAGMENT_LIB:
            if rf'[{i}*]' in smi:       # 如果碎片库中有 [i*]
                if smi not in lib_avail:  # bug fixed by zhouzx
                    lib_avail.append(smi)

    return lib_avail    # 返回一个 smi 列表


def mol_connect(par, h1, frag, frag_id, smile):  # par:已生成的部分（在 h1 的 “根” 位置待连，h1 是虚拟原子）， frag:待装上的片段
    '''
    Connect new molecular fragments
    '''
    @func_set_timeout(5)  # 该装饰器用于设置函数的超时时间为 5 秒。这意味着如果 mol_connect2 函数的执行时间超过 5 秒，将引发超时异常
    def mol_connect2(par, h1, frag, frag_id, smile):
        for atom in par.GetAtoms():        # par 应该代表目前 partial 片段
            atom.SetIntProp('FromPar', 1)  # 为原子设置 'FromPar' 属性为 1，未来告诉连接后来自 par
            if atom.GetIsotope() != 0 and atom.GetIdx() == h1:  # 同位素标记非 0，且 待连接的原子 id 值为 h1
                Iso_type = atom.GetIsotope()
                atom.GetNeighbors()[0].SetIntProp('Nei', 1)
                atom.SetAtomicNum(37)      # 37 号元素是 Rb 铷
            if atom.GetIsotope() != 0 and atom.GetIdx() != h1:  # 把其它剩余连接口用 H 封上
                atom.SetAtomicNum(1)       # 设置为 H

        Iso_types = Iso_dic[Iso_type]

        for atom in frag.GetAtoms():
            if atom.GetIsotope() in Iso_types:
                atom.SetAtomicNum(87)     # Fr 钫
                atom.GetNeighbors()[0].SetIntProp('Nei', 1)
                break   # 就只选第一个？？？

        for atom in frag.GetAtoms():
            atom.SetIntProp('FromPar', 0)
            if atom.GetIsotope() != 0 and atom.GetAtomicNum() == 0:   # 把其它剩余连接口用 H 封上
                atom.SetAtomicNum(1)    # 设置为 H

        combo = ch.CombineMols(par, frag)
        Rb_neiid = get_neiid_bysymbol(combo, 'Rb')
        Fr_neiid = get_neiid_bysymbol(combo, 'Fr')
        edcombo = ch.EditableMol(combo)
        edcombo.AddBond(Rb_neiid, Fr_neiid, order=Chem.rdchem.BondType.SINGLE)

        Rb_index = get_id_bysymbol(combo, 'Rb')
        edcombo.RemoveAtom(Rb_index)    # 删掉 Rb
        back = edcombo.GetMol()

        Fr_index = get_id_bysymbol(back, 'Fr')
        edcombo = Chem.EditableMol(back)
        edcombo.RemoveAtom(Fr_index)    # 删掉 Fr
        back = edcombo.GetMol()  # back 是最终片段

        Dihedral = []   # 用于储存二面角的四个原子的 id
        for at in back.GetAtoms():
            if at.GetProp('FromPar') == '1' and at.HasProp('Nei') == 1:  # 是虚拟原子 h1 的邻居， h1 的邻居才是真原子
                for nei in at.GetNeighbors():
                    if nei.GetProp('FromPar') == '1':
                        Dihedral.append(str(nei.GetIdx()))  # 用于储存二面角的四个原子的 id
                        Dihedral.append(str(at.GetIdx()))
                        break
                break

        for at in back.GetAtoms():
            if at.GetProp('FromPar') == '0' and at.HasProp('Nei') == 1:
                for nei in at.GetNeighbors():
                    if nei.GetProp('FromPar') == '0':
                        Dihedral.append(str(at.GetIdx()))
                        Dihedral.append(str(nei.GetIdx()))
                        break
                break

        par_back = par   # par 一直在通过 atom 被修改
        par_back.GetAtomWithIdx(get_id_bysymbol(par_back, 'Rb')).SetAtomicNum(1)  # 把 Rb 换成 H
        par_back_rm_H = Chem.RemoveHs(par_back)  # par 骨架
        for at in par_back_rm_H.GetAtoms():
            at.SetIsotope(0)   # 都设为不可连接
        par_back_rm_H = Chem.RemoveHs(par_back_rm_H)
        par_index = par_back.GetSubstructMatch(par_back_rm_H)  # 匹配出来原子 id  ???
        combo_index = back.GetSubstructMatch(par_back_rm_H)    # back 是最终片段
        try:
            cmap = {combo_index[j]: par_back.GetConformer().GetAtomPosition(par_index[j]) for j in
                    range(len(par_index))}    #  cmap 是 coordMap 坐标映射，用于叠合
        except:
            traceback.print_exc()
            print(
                [par_index, combo_index, ch.MolToSmiles(par_back), ch.MolToSmiles(back), ch.MolToSmiles(par_back_rm_H)]
            )                           # back 是最终片段
        cids = AllChem.EmbedMultipleConfs(back, numConfs=30, coordMap=cmap, maxAttempts=1000, numThreads=4,
                                          randomSeed=1)    # 嵌入 30 个构象
        tag = 1
        if len(cids) > 0:

            rms_min = 10
            for i in range(len(cids)):   # 对每个构象 i:

                rms = rdMolAlign.AlignMol(back, par_back, prbCid=i, atomMap=list(zip(combo_index, par_index)))   # 叠合
                if rms < 2:
                    if rms < rms_min:
                        rms_min = rms                        
                        writer = Chem.SDWriter(rf'tmp/{frag_id}.sdf')
                        writer.SetProps(['DihAtoms', 'DihDeg'])
                        back.SetProp('DihAtoms', ','.join(Dihedral))  # 记录二面角原子 id
                        Deg = rdMolTransforms.GetDihedralDeg(back.GetConformer(id=i), int(Dihedral[0]),
                                                             int(Dihedral[1]),
                                                             int(Dihedral[2]), int(Dihedral[3]))   # 获得二面角
                        back.SetProp('DihDeg', str(Deg))
                        writer.write(back, confId = i )
                        writer.close()
                        os.system(rf"grep -Ev 'MDLV30|COLLECTION' ./tmp/{frag_id}.sdf > ./tmp/{frag_id}_tt.sdf")
                        os.system(rf'cp ./tmp/{frag_id}_tt.sdf ./tmp/{frag_id}.sdf')
                        tag = 0       # 成功 0
        return [frag_id, tag, smile]

    try:
        a = mol_connect2(par, h1, frag, frag_id, smile)
        return a
    except:
        return [frag_id, 1, smile]   # 出现异常 1


def mol_rotate(frag_id, Del_plus, smile):    # Del_plus 是增加的角度
    '''
    Rotate newly introduced dihedral angles to get multiple molecular conformations
    '''
    mol = ch.SDMolSupplier(rf'tmp/{frag_id}.sdf', removeHs=False)[0]
    Dihedral_atoms = mol.GetProp("DihAtoms").split(',')
    Deg = float(mol.GetProp("DihDeg")) + Del_plus       # DihDeg 是原二面角大小
    Comf = mol.GetConformer()
    rdMolTransforms.SetDihedralDeg(Comf, int(Dihedral_atoms[0]), int(Dihedral_atoms[1]), int(Dihedral_atoms[2]),
                                int(Dihedral_atoms[3]), Deg)
    writer = Chem.SDWriter(rf'tmp/{frag_id}_{Del_plus}.sdf')
    writer.write(mol)
    writer.close()


def score0(frag_id, Del_plus, smile):
    '''
    Determine whether newly connected fragments can cause collisions between atoms
    '''
    time1 = time.time()
    lig = ch.SDMolSupplier(rf'tmp/{frag_id}_{Del_plus}.sdf', removeHs=False)[0]
    for atom in lig.GetAtoms():
        if atom.GetAtomicNum() == 37:
            atom.SetAtomicNum(1)         # 如果还有未连接的，则先换成 H
    lig_name = rf'{frag_id}_{Del_plus}'
    L = lig.GetConformer()
    lig_cors = [L.GetPositions()[i] for i in range(len(lig.GetAtoms()))]

    dis_matrix = cdist(POCKET_CORS, lig_cors, metric='euclidean')            #  距离矩阵
    pt = Chem.GetPeriodicTable()
    radis_matrix = np.zeros((len(POCKET.GetAtoms()), len(lig.GetAtoms())))   #  储存蛋白中原子 i 和配体原子 j 之间的共价原子半径之和
    for i in range(len(POCKET.GetAtoms())):
        at1 = pt.GetRcovalent(POCKET.GetAtomWithIdx(i).GetAtomicNum())
        for j in range(len(lig.GetAtoms())):
            at2 = pt.GetRcovalent(lig.GetAtomWithIdx(j).GetAtomicNum())
            radis_matrix[i][j] = at1 + at2     # 共价半径之和
    judge = (radis_matrix > dis_matrix).sum()

    lig_intra_matrix = cdist(lig_cors, lig_cors, metric='euclidean')
    covalent = np.ones((len(lig.GetAtoms()), len(lig.GetAtoms())))   # 都是 1
    bonds = lig.GetBonds()
    for bd in bonds:
        idx1 = bd.GetBeginAtomIdx()
        idx2 = bd.GetEndAtomIdx()
        covalent[idx1][idx2] = 0      # covalent 储存成键原子对， 成键是 0，不影响原子之间距离小
        covalent[idx2][idx1] = 0
        covalent[idx1][idx1] = 0
        covalent[idx2][idx2] = 0
    lig_intra_radis = np.zeros((len(lig.GetAtoms()), len(lig.GetAtoms())))
    for i in range(len(lig.GetAtoms())):
        at1 = pt.GetRcovalent(lig.GetAtomWithIdx(i).GetAtomicNum())
        for j in range(len(lig.GetAtoms())):
            at2 = pt.GetRcovalent(lig.GetAtomWithIdx(j).GetAtomicNum())
            lig_intra_radis[i][j] = at1 + at2       #  原子对共价半径之和
    judge2 = ((lig_intra_matrix < lig_intra_radis) * covalent).sum()    # 距离 < 半径
    # print('score0:',time.time()-time1)

    return [frag_id, Del_plus, judge + judge2, smile]


def score1(frag_id, Del_plus, smile):   # 将配体的 .sdf 中 Rb 换成 H，存成 .sdf， 再转换成 .pdbqt
    '''
    Prepare files for scoring.
    '''
    time1 = time.time()
    pro = 'pro_test.pdbqt'
    mol = ch.SDMolSupplier(rf'tmp/{frag_id}_{Del_plus}.sdf', removeHs=False)[0]
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 37:
            atom.SetAtomicNum(1)
    lig_name = rf'{frag_id}_{Del_plus}'
    writer = Chem.SDWriter(rf'vinaa/{lig_name}_score.sdf')
    writer.write(mol)
    writer.close()
    mymol = list(readfile('sdf', rf'vinaa/{lig_name}_score.sdf'))[0]
    mymol.write("pdbqt", 'vinaa/' + lig_name + '.pdbqt', overwrite=True)   #  {frag_id}_{Del_plus}.pdbqt


def score3(frag_id, Del_plus, smile):   # 使用 GNINA 打分 亲和能
    '''
    Obtain the binding affinity of the molecule
    '''
    time1 = time.time()
    pro = 'pro_test.pdbqt'
    lig_name = rf'{frag_id}_{Del_plus}'

    try:
        score = float(os.popen(
            rf'{GNINA}/gnina --minimize --receptor {pro} --ligand vinaa/{lig_name}.pdbqt --cnn_scoring=none --out ./vinaa/{lig_name}_minimize.sdf').read().split(
            'Affinity: ')[1].split(' (kcal/mol)')[0].split()[0])
    except:
        score = 0

    return [frag_id, Del_plus, score, smile]


def repair(frag_id, Del_plus, score, smile):
    '''
    Repair the chemical valence of atoms after minimization
    '''
    lig_name = rf'{frag_id}_{Del_plus}'
    f = open(rf'./vinaa/{lig_name}_minimize.sdf', 'r').readlines()
    f_ = open(rf'./vinaa/{lig_name}_repair.sdf', 'w')
    for line in f[:4]:                                    # 前四行不变
        f_.write(line)
    for i in range(4, len(f)):
        if f[i].startswith(' ') and len(f[i]) >= 50:
            new_line = f[i][:50] + '0' + f[i][51:]        # 在第 50 和第 51 个字符之间插入一个 '0'
            f_.write(new_line)
        else:
            br = i
            break
    for i in range(br, len(f)):                           # br 到最后行，不变
        f_.write(f[i])
    f_.close()

    try:
        assign_iso(frag_id, Del_plus)        # 赋 同位素 标记
    except:          
        mol_parent = Chem.SDMolSupplier(rf'./tmp/{lig_name}.sdf', removeHs=False)[0]
        for at in mol_parent.GetAtoms():
            if at.GetIsotope() != 0:
                at.SetIsotope(0)    # 同位素标记全部删除
        writer = Chem.SDWriter(rf'./tmp/{lig_name}.sdf')
        writer.write(mol_parent)
        writer.close()
        os.system(rf"grep -Ev 'MDLV30|COLLECTION' ./tmp/{lig_name}.sdf > ./tmp/{lig_name}_tt.sdf")
        os.system(rf'cp ./tmp/{lig_name}_tt.sdf ./tmp/{lig_name}.sdf')
        traceback.print_exc()


def qed_score(frag_id, Del_plus, score, smile):
    '''
    Obtain the drug-like properties of a molecule
    '''
    lig_name = rf'{frag_id}_{Del_plus}'
    mol = Chem.SDMolSupplier(rf'./tmp/{lig_name}.sdf', removeHs=False)[0]
    qed = QED.qed(mol)    # 计算 QED
    return [frag_id, Del_plus, score, smile, qed]


def roulette(select_list):    # 轮盘赌算法，我觉得原代码写的不对
    '''
    roulette algorithm
    '''
    sum_val = sum(select_list)
    random_val = random.random()
    probability = 0
    if sum_val != 0:
        for i in range(len(select_list)):
            probability += select_list[i] / sum_val
            if probability >= random_val:
                return i      # 返回第一个 i 后，就会跳出函数，结束
            else:
                continue
    else:
        return random.choice(range(len(select_list)))




class State():  # 态里面存有： (self.type, self.sdf, self.Frag_Deg, self.h1s_avail/self.h1, self.score=sco/0, self.states = sta+[self.score],
                #              self.ter_tag = 0 0代表未终结, self.frags_avail(里面是smi),) 在 State().next_states() 里面： 连接，旋转，打分，操作全局变量mol_dic
    '''        # 节点的状态。在3D-MCTS中有两种类型的节点, 一种用于确定片段连接位置(类型 0), 另一种用于确定片段的类型和片段的构象(类型 1).
    The status of the node. There are two types of nodes in 3D-MCTS, one is used to determine the connection position of fragments (type 1),
    and the other is used to determine the type and conformation of fragments (type 0).
    '''                                                                                # sco 是打分 score   # ter_tag 是 MCTS node 的 terminate tag = 1 终结
    def __init__(self, state_type=0, sdf='start.sdf', h1=None, Frag_Deg=None, frag=None, sco=-3.55479, ter_tag=None,
                 sta=[], choices=[]):
        self.type = state_type
        self.sdf = sdf

        if self.type == 0:               # 确定片段连接位置 h1s

            self.Frag_Deg = Frag_Deg
            self.h1s_avail = []
            mol = ch.SDMolSupplier(self.sdf, removeHs=False)[0]
            if self.sdf == 'start.sdf':
                os.system(rf'cp start.sdf state0.sdf')
            self.score = sco
            for atom in mol.GetAtoms():
                if atom.GetIsotope() != 0:
                    self.h1s_avail.append(atom.GetIdx())    # 存储可连接的位置
            self.states = sta + [self.score]           # 储存得分
            self.choices = choices + [self.Frag_Deg]   # 储存片段 id 和角度

        elif self.type == 1:           # 确定 h1 可连的 片段的类型和片段的构象   

            self.h1 = h1

            self.score = 0
            self.states = sta + [self.score]        # sta=[]
            self.choices = choices + [self.h1]      # choices=[]
            self.ter_tag = 0

    def next_states(self):
        '''
        Get all possibilities for the next state of the current state.

        '''
        if self.type == 0:      #  确定头片段连接位置 h1s
            pass

        elif self.type == 1:    # 确定 h1 可连的 片段的类型和片段的构象
            mol = ch.SDMolSupplier(self.sdf, removeHs=False)[0]
            h1_isotop = mol.GetAtomWithIdx(self.h1).GetIsotope()    # self.h1: 头片段连接点 id 
            try:
                self.frags_avail = frag_avail(h1_isotop)     # 查连接规则字典 Iso_dic，self.frags_avail：是一个 可连接的片段的smi 列表
            except:
                print(self.sdf, self.h1, h1_isotop)     # ( 头片段sdf, 头片段连接点id, 头片段连接点同位素类型)

            frags = []
            mols = []
            h1s = []
            ids = []
            j = 0
            time1 = time.time()
            for i in self.frags_avail:      # 对一个 h1 的多个可连接的片段的smi
                frags.append(ch.AddHs(ch.MolFromSmiles(i)))    # 尾片段
                mols.append(mol)            # 头片段
                h1s.append(self.h1)         # 头片段连接点 id 
                j += 1
                ids.append(j)
            # fragment connection
            pool = multiprocessing.Pool(args.processor)     # 创建了一个多进程池：进程数量 = 处理器的个数 args.processor
            ids_tags_smiles = pool.starmap(mol_connect, zip(mols, h1s, frags, ids, self.frags_avail)) # 将 mol_connect 函数应用到由 zip 生成的迭代器中的每个元组上
            # ids_tags_smiles：[[frag_id, tag, smile],...]   tag=0 成功， tag=1 失败
            pool.close()    # 关闭进程池，表示不需要继续提交新任务，防止额外的资源浪费，避免出现潜在的内存泄漏或资源泄漏
            pool.join()     # 等待所有的子进程结束，确保进程池中的所有任务都已完成
            # print(time.time()-time1)
            time1 = time.time()
            ids_degs_smiles = []
            for id_tag_smile in ids_tags_smiles:
                if id_tag_smile[1] == 0:
                    for deg in np.arange(0, 360, 15):
                        ids_degs_smiles.append([id_tag_smile[0], deg, id_tag_smile[2]])

            # Rotate dihedral angles to obtain multiple conformations
            pool = multiprocessing.Pool(args.processor)
            pool.starmap(mol_rotate, ids_degs_smiles)     # 应用 mol_rotate 函数进行旋转，结果直接生成在：tmp/{frag_id}_{Del_plus}.sdf
            pool.close()
            pool.join()
            # print(time.time() - time1)
            # nextmove = random.choice([x for x in self.frags_avail])

            time1 = time.time()
            pool = multiprocessing.Pool(args.processor)
            ids_degs_judge_smiles = pool.starmap(score0, ids_degs_smiles)  # 应用 score0 进行原子碰撞评估
            pool.close()
            pool.join()
            # print(time.time() - time1)

            ids_degs_smiles = [[i[0], i[1], i[3]] for i in ids_degs_judge_smiles if i[2] == 0]  # 如果没有原子冲突，i[2] = 0
            if len(ids_degs_smiles) > 0:   # 存在不冲突的构象
                pool = multiprocessing.Pool(args.processor)
                pool.starmap(score1, ids_degs_smiles)    # 将配体的 .sdf 中 Rb 换成 H，存成 .sdf， 再转换成 {frag_id}_{Del_plus}.pdbqt
                pool.close()
                pool.join()

                pool = multiprocessing.Pool(args.processor)
                ids_degs_scores_smiles = pool.starmap(score3, ids_degs_smiles)   # 使用 GNINA 最小化(但不是对接！)， 打分 亲和能
                pool.close()
                pool.join()
                ids_degs_scores = [i for i in ids_degs_scores_smiles if i[2] <= (self.states[-2] + 0.5)]   # 选择能量上升不高的

                # Sort multiple conformations according to binding affinity
                ids_degs_scores_smiles = ids_degs_scores
                ids_degs_scores_smiles_selected = []
                ids_degs_scores_sorted = sorted(ids_degs_scores, key=lambda x: x[2])   # 对 ids_degs_scores 中 score 排序
                ids_degs_scores_smiles_dic = {}
                for i in ids_degs_scores_sorted:
                    try:
                        if len(ids_degs_scores_smiles_dic[i[0]]) < 1:      # 如果存在键 id: i[0]，但是值不存在
                            ids_degs_scores_smiles_dic[i[0]].append(i)     # 添加到字典的键为 i[0] 的值中，值本身是一个列表
                        else:
                            continue                                       # 已经存在键和值，那么已经记录了最小的能量
                    except:
                        ids_degs_scores_smiles_dic[i[0]] = [i]             # 不存在键，即未记录这个分子的 [id, deg, score]
                for value in ids_degs_scores_smiles_dic.values():
                    ids_degs_scores_smiles_selected.extend(value)    # extend() 方法用于在列表末尾一次性追加另一个可迭代对象（通常是另一个列表）中的所有元素


                # Repair the structure after minimization
                pool = multiprocessing.Pool(args.processor)
                pool.starmap(repair, ids_degs_scores_smiles_selected)   # score3 使用 GNINA 最小化(但不是对接！)，也输出了结构，所以需要修复一下
                pool.close()
                pool.join()

                # Save the molecules that match the criteria 判据：args.score，args.qed
                for i in ids_degs_scores_smiles_selected:
                    if i[2] < (args.score):  # 小于用户设定的 gnina 能量 -7.0 kcal/mol
                        mol = ch.SDMolSupplier(rf'tmp/{i[0]}_{i[1]}.sdf', removeHs=False)[0]  # .sdf 有同位素标记的
                        global RECORD
                        for at in mol.GetAtoms():
                            at.SetIsotope(0)   # 把同位素标记全删了
                        writer = Chem.SDWriter(rf'tmp/{i[0]}_{i[1]}_.sdf')                    # _.sdf 没有同位素标记
                        writer.write(mol)
                        writer.close()
                        mol = ch.SDMolSupplier(rf'tmp/{i[0]}_{i[1]}_.sdf', removeHs=False)[0]
                        if Chem.MolToSmiles(mol) not in mol_dic:
                            mol_dic.append(Chem.MolToSmiles(mol))    # 记录下 smi
                            qed = QED.qed(mol)
                            NO_QED.append(i[2])   # 非 QED 的 gnina 分
                            if qed > args.qed:     # qed 越大越好
                                writer = Chem.SDWriter(rf'record/record_{RECORD}.sdf')   # 记录下 gnina 分，qed 都好的分子结构
                                writer.SetProps(['Score'])
                                mol.SetProp('Score', str(i[2]))  # 分子 .sdf 文件中也记录了 score3（gnina分）
                                writer.write(mol)
                                writer.close()
                                RECORD += 1     # 递增
                                SCORES.append(i[2])      # 里面是二面角旋转排名第一的，且满足 args.score，args.qed 的分子的得分
                                if i[2] < (args.score + 0):
                                    GOOD_SCORES.append(i[2])  # GOOD_SCORES主要是数据分析时用，代码中没有引用
                            if i[2] < (args.score + 0):
                                NO_QED_GOOD.append(i[2])    # NO_QED_GOOD主要是数据分析时用，代码中没有引用

                if False:
                    self.ids_degs_scores = ids_degs_scores_smiles[:50]
                else:     # 必定进入一下这段代码：
                    pool = multiprocessing.Pool(args.processor)
                    ids_degs_scores = pool.starmap(qed_score, ids_degs_scores_smiles_selected)   # return [frag_id, Del_plus, score, smile, qed]
                    pool.close()
                    pool.join()


                    self.ids_degs_scores = [i for i in ids_degs_scores if i[-1] > args.qed]   #  [frag_id, Del_plus, score, smile, qed]

                if len(self.ids_degs_scores) == 0:
                    self.ter_tag = 1  # 态终结
                    print('node_end')    # 这个 MCTS node 的分子都过大，qed 全部不能满足 > args.qed
                    return []
                else:
                    self.ter_tag = 0
                    print('expand')
                    return self.ids_degs_scores
            else:
                self.ter_tag = 1     # 态终结，全是冲突的构象
                return []

    def terminal(self, tag):

        # Determine whether the current state meets the termination conditions
        if self.type == 1:
            if self.ter_tag == 1:
                return True
            else:
                return False


        else:
            mol = ch.SDMolSupplier(self.sdf, removeHs=False)[0]

            if len(self.h1s_avail) == 0:
                return True

            if self.score > 0:
                self.stop_normal = 0   # 不正常停止
                return True
            else:
                try:                 # 不一定能有 states[-3]                     #  i[2] <= (self.states[-2] + 0.5)]
                    if (self.score - 1.5) > self.states[-3]:    # 能量过高，相当于连续 3 次能量不下降: 1.5 = 3 * 0.5
                        return True
                    else:                                            # CalcNumLipinskiHBA用于计算分子中符合Lipinski规则的氢键受体数目
                        if rdMolDescriptors.CalcNumLipinskiHBA(      # CalcNumLipinskiHBD用于计算分子中符合Lipinski规则的氢键供体数目
                                mol) > 9 or rdMolDescriptors.CalcNumLipinskiHBD(mol) > 4 or Descriptors.MolWt(
                            mol) > 500:
                            self.stop_normal = 1
                            return True
                        else:
                            return False   # 不终结
                except:
                    if rdMolDescriptors.CalcNumLipinskiHBA(
                            mol) > 9 or rdMolDescriptors.CalcNumLipinskiHBD(mol) > 4 or Descriptors.MolWt(mol) > 500:
                        self.stop_normal = 1
                        return True
                    else:
                        return False

    def reward(self):   # 没被调用过 ！
        r1 = 0            # 里面是二面角旋转排名第一的，且满足 args.score，args.qed 的分子的得分
        r2 = max((np.mean(SCORES) - self.score2) / abs(np.mean(SCORES + [self.score2])),    0)   # 主要是看 (np.mean(SCORES) - self.score2) 的 +/-
        if self.score2 < np.mean(SCORES):
            SCORES.append(self.score2)
        print(self.score, self.score2)
        print(r1, r2)
        return r1, r2

    def __hash__(self):
        if self.type == 0:
            return int(hashlib.md5(str(self.Frag_Deg).encode('utf-8')).hexdigest(), 16)
        elif self.type == 1:
            return int(hashlib.md5(str(self.h1).encode('utf-8')).hexdigest(), 16)
    def __eq__(self, other):
        if hash(self) == hash(other):
            return True
        return False

class Node():                                       # best_score 在每次调用时都会给，而且会不断在模拟中更新，最初始值给的是起始片段 gnina 对接打分
    def __init__(self, state, parent=None, node_id=1, best_score=-3.55479, reward=0, qed=1):
        self.visits = 0
        self.reward = reward
        self.state = state
        self.children = []
        self.parent = parent
        self.best_score = best_score
        self.longest_path = 0
        self.id = node_id
        self.qed = qed
                         # child_state 是一个 State() 类的实例化对象
    def add_child(self, child_state, node_id, bestscore, qed=1):  # bestscore 是以这个节点往下延生出的节点的最小的 score
        sdf_name = child_state.sdf     # State() 类的属性 sdf
        os.system(rf'cp {sdf_name} node/{node_id}.sdf')
        child_state.sdf = rf'node/{node_id}.sdf'
        child = Node(child_state, node_id=node_id, parent=self, best_score=bestscore, qed=qed)
        self.children.append(child)

    def update(self, reward):   # 没调用过 ！
        self.reward += reward
        self.visits += 1
                           # num_moves_lambda=None
    def fully_expanded(self, num_moves_lambda):
        if self.state.type == 0:
            num_moves = len(self.state.h1s_avail)    # 选哪个 h1
        elif self.state.type == 1:
            num_moves = len(self.state.ids_degs_scores)    # 连哪个片段，什么角度
        if num_moves_lambda != None:
            num_moves = num_moves_lambda(self)
        if len(self.children) == num_moves:
            return True
        return False

    def __repr__(self):
        s = "Node; children: %d; visits: %d; reward: %f" % (len(self.children), self.visits, self.reward)
        return s

            # budget = args.num_sims
def UCTSEARCH(budget, root, start_score=0, num_moves_lambda=None):  # UCTSEARCH(args.num_sims, current_node, start_score=results[i][2]) 
    # Begin the MCTS
    start = time.time()
    global GOOD_SCORES
    global NO_QED_GOOD

    for iter in range(int(budget)):   # budget=args.num_sims
        print('\n\n\n\n')
        print('='*30)
        print('='*30)
        print('='*30)
        print('\n'+rf'########## UCTSEARCH iter/args.num_sims : {iter} / {budget}, current node: {root.id}')
        front = TREEPOLICY(root, start_score, num_moves_lambda)
        BACKUP2(front)

                    # start_score = results[i][2]
def TREEPOLICY(node, start_score, num_moves_lambda):
    # Choose whether to expand the node based on the status of the current node

    print('\n'+rf'''
          #######################################################
          ############# TREEPOLICY while loop start #############
          #######################################################
          ''')
    cc = 0
    while node.state.terminal('state1.sdf') == False:     # while node 状态一直不终结， 我怀疑这里不要 'state1.sdf' 也行 ？？？
        
        print('\n'+'############# TREEPOLICY node.state.terminal() == False')
        print(rf'''############# TREEPOLICY while 【loop {cc}】:''')
        cc += 1
        print_node_info(node)        
            
        if len(node.children) == 0:  
            print('\n'+rf'------------TREEPOLICY (Run EXPAND) :')          
            node = EXPAND(node, start_score)    # 更新 node
            
        else:       
            print('\n'+rf'----------TREEPOLICY (CHOOSE BESTCHILD Using Roulette) :')     
            node = BESTCHILD(node, SCALAR, start_score)  # 更新 node
            print_node_info(node)
    print('\n'+'############# TREEPOLICY node.state.terminal() = False 终结')
    print('\n'+rf'''
          #######################################################
          ############# TREEPOLICY while loop end ###############
          #######################################################
          ''')
    return node


def EXPAND(node, start_score):

    # Get the children of a node and add them to the tree    感觉像在 模拟
    global NODE_ID
    if node.state.type == 0:    # 把每一个可用的 h1 都选一遍，做好准备，不需要大量计算！！！
        print(rf'node.state.h1s_avail: {node.state.h1s_avail}')
        for nextmove in node.state.h1s_avail:    # 从 h1s 中选取下一个可用的 h1         State().states = sta + [self.score]  # 储存得分
            next = State(state_type=1, sdf=rf'{node.state.sdf}', h1=nextmove, sta=node.state.states) # 从 type=0 跳到 =1           
            NODE_ID += 1
            node.add_child(next, node_id=NODE_ID, bestscore=start_score)     # start_score = results[i][2]
            print('\n----expanding out the node with state_type=【1】by choosing h1  【fast!】:')
            print_node_info(node.children[-1])   
        return node.children[-1]
    elif node.state.type == 1:
        new_states = node.state.next_states()     # .next_states()需要大量计算！！！（确定h1可连的 所有片段 和 旋转角度）
        if len(new_states) == 0:
            return node  # 节点本身状态已更新，在下次 TREEPOLICY 中  while node.state.terminal()终结
        else:
            scores = []
            print(rf'for nextmove in new_states # ids_degs_scores : {new_states}')
            for nextmove in new_states:    # ids_degs_scores                
                os.system(rf'cp tmp/{nextmove[0]}_{nextmove[1]}.sdf state0.sdf')
                next = State(state_type=0, sdf=rf'state0.sdf', Frag_Deg=nextmove[:2], sco=nextmove[2],   # 从 type=1 跳到 =0
                             sta=node.state.states)                
                NODE_ID += 1
                best_score = min(start_score, nextmove[2])   # best_score 是以这个节点往下延生出的节点的最小的 score
                scores.append(abs(nextmove[4]))
                node.add_child(next, node_id=NODE_ID, bestscore=best_score, qed=abs(nextmove[4]))
                print('\n----expanding out the node with state_type=【0】by node.state.next_states()  【slow!】:')
                print_node_info(node.children[-1])
                
            rt = node.children[roulette(scores)]

            print('\n'+rf'&&&& Rouletting out the node with state_type=【0】:')
            print_node_info(rt)

            return rt

def BESTCHILD(node, scalar, start_score):

    # Select child nodes based on the node's UCB
    scores = []
    scores2 = []
    for c in node.children:

        exploit = start_score - c.best_score   # 正值
        explore = math.sqrt(2.0 * math.log(node.visits + 0.000001) / float(c.visits + 0.000001))   # UCB 探索项
        score = exploit + scalar * explore    # UCB 得分
        scores.append(score)
        scores2.append(c.qed)
    if True:
        idx1 = roulette(scores)
        idx2 = roulette(scores2)
        idx = random.choice([idx1, idx2])    #  在 UCB 得分 和 QED 两个因素影响下，随机选取 子节点
    else:
        idx = random.choice(range(len(scores)))
    return node.children[idx]


def DEFAULTPOLICY(node):    # 并没有被调用过 ！！！，是模拟步，直接被 TREEPOLICY 中的 EXPAND 模拟取代了
    state = node.state
    num_states = 0

    while state.terminal('state1.sdf') == False:   # 这里的 'state1.sdf' 或许可以不要 ？？？
        state = state.next_state()     # 不停的延生，模拟，跟新态，直到达到终止条件，成为一个叶节点
        num_states += 1
    if state.type == 1:       # 确定 构象 和 旋转角度
        if num_states != 0:   # 这个节点不是叶节点
            num_states -= 1
    num_nodes = len(state.states) - num_states    # state.states = sta + [self.score]    # 储存得分
    print(state.type)
    return state.states, num_nodes, num_states


def BACKUP2(node):   # 反向跟新节点的访问次数
    '''
    After a search is completed, backpropagation is performed to update the number of visits to the node.
    '''
    parent_node = node
    son_node = None  # 初始化子节点
    while parent_node != None:        # 不停的回溯父节点
        parent_node.visits += 1       # 给访问记录加 1
        if len(parent_node.children) == 0:  # 如果是叶节点，应该每次调用 BACKUP2() 里面放的都是叶节点
            x = parent_node
            parent_node = node.parent  # 向上寻找 父节点
            son_node = x
        else:
            if son_node is not None:          # 如果存在子节点
                if parent_node.best_score > son_node.best_score:
                    parent_node.best_score = son_node.best_score   # 更新父节点的最佳分数
                x = parent_node
                parent_node = parent_node.parent  # 向上寻找 父节点
                son_node = x

        # node 应该是一个叶节点
def BACKUP(node, states, num_nodes, num_states):    # 没被调用过 ！！！
    print(states)
    i = 1
    if node.longest_path == 0:
        node.longest_path = len(states)
    while node != None:
        node.visits += 1                           # 回溯 node, 分别计算访问次数, reward
        best_score = min(states[num_nodes - i:])   # states[num_nodes: ], states[num_nodes-1: ], states[num_nodes-2: ]
        i += 1
        if best_score not in SCORES:
            if best_score < node.best_score:
                node.best_score = best_score
                reward = max(-3.55479 - best_score, 0)   # 只有比 -3.55479 更低才有奖励
            else:
                reward = 0
            if best_score < np.mean(SCORES):   # 小于 SCORES的平均
                SCORES.append(best_score)      # 添加到 SCORES
        else:
            if best_score < node.best_score:
                node.best_score = best_score
                reward = max(-3.55479 - best_score, 0)
            else:
                reward = 0
        node.reward += reward
        node = node.parent
    return


if __name__ == "__main__":

    results = []
    i = args.start
    # Obtain binding modes of starting fragments using molecular docking  只在第一步通过 gnina 对接获得构象
    results.extend(dock(i))      #  results = [[args.start, j, score], ...]    j: 0 - 9  一共产生 9 个对接构象
    # Choose the best docking conformation
    results.sort(key=lambda results: results[2])  # 根据 gnina 对接分排序
    for i in range(len(results)):   # 对 9 个起始对接 pose 都延生
        # Label starting fragments with isotopes for new fragment ligation
        assign_iso(results[i][0], results[i][1], start=True)     # (args.start, j, start=True)
        # 经过 assign_iso，对 args.start 结构中各原子打上同位素标记，并存成 .sdf 文件
        current_node = Node(State(sta=[], sco=results[1][2]), best_score=results[i][2])   # 能比 sco=(第 2 个对接 pose 的分高)，就能有 reward，只有 BACKUP 里面真正用到了 reward，别的地方的 reward 都没被调用
        result = UCTSEARCH(args.num_sims, current_node, start_score=results[i][2])   # 返回新的一轮的

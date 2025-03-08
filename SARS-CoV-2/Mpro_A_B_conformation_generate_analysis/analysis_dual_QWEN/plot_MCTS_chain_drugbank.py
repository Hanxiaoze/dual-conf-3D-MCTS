import os
import shutil
from openbabel import pybel
import ast
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
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
from openbabel import openbabel 
import time
import multiprocessing, traceback
from scipy.spatial.distance import cdist
from func_timeout import func_set_timeout
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import QED
from rdkit.Chem import Descriptors


parser = argparse.ArgumentParser(description='3D-MCTS code, for molecular generation')
parser.add_argument('--frag_lib', action="store", type=str, help='fragment library', default='../../../frags/my_frags_2.txt')
args = parser.parse_args()

FRAGMENT_LIB = []
FRAGMENT_LIB_no_iso = []
with open(f'{args.frag_lib}','r') as f:
    lines = f.readlines()
    for line in lines:
        FRAGMENT_LIB.append(line.strip())

mcts_fragment_to_idx = {fragment: idx for idx, fragment in enumerate(FRAGMENT_LIB)}
mcts_idx_to_fragment = {idx: fragment for idx, fragment in enumerate(FRAGMENT_LIB)}

for frag_smi in FRAGMENT_LIB:
    frag_mol = Chem.MolFromSmiles(frag_smi)
    frag_mol_H = Chem.AddHs(frag_mol, addCoords=True)
    try:
        for a in frag_mol_H.GetAtoms():
            if a.GetIsotope() != 0:            
                a.SetIsotope(0)
                a.SetAtomicNum(1)
    except:
        print(rf'frag : {frag_smi}  failed in 3D conform!')
    smi_new = Chem.MolToSmiles(frag_mol_H, isomericSmiles=False)
    smi_new = Chem.MolToSmiles(Chem.MolFromSmiles(smi_new))
    FRAGMENT_LIB_no_iso.append(smi_new)

print(FRAGMENT_LIB_no_iso)
print(len(FRAGMENT_LIB_no_iso))

mcts_noiso_fragment_to_idx = {fragment: idx for idx, fragment in enumerate(FRAGMENT_LIB_no_iso)}
mcts_noiso_idx_to_fragment = {idx: fragment for idx, fragment in enumerate(FRAGMENT_LIB_no_iso)}

print(mcts_noiso_idx_to_fragment)


# print(f'\nmcts_fragment_to_idx: \n{mcts_fragment_to_idx}')
# print(f'len of mcts_fragment_to_idx: {len(mcts_fragment_to_idx)}')

datasmi_mol_elem = []
datasmi_mol_2nd = []
datasmi_mol_total = []
data_mol_elem = []
data_mol_2nd = []
data_mol_total = []
data_1_path = './frags/mol_in_chem_space.txt'
data_2_path = './frags/mol_not_in_chem_space.txt'


elem_frag_smi_list = []
frag_smi_to_2nd_break_list = []
frag_smi_2nd_break_dict = {}


with open("../../../frags/my_frags_3.txt", "r") as f:
    frag_smi_to_2nd_break_list = [line.split(':')[0].strip() for line in f.readlines()[:109] if ':' in line]

with open("../../../frags/my_frags_3.txt", "r") as f:
    elem_frag_smi_list = [line.strip() for line in f.readlines()[:109] if ':' not in line]
    elem_frag_smi_list.append('xxxx')   # 65 个片段库  + EOS
    elem_frag_smi_list.append('EOS')

# print(f'\nelem_frag_smi_list: \n{elem_frag_smi_list}')
# print(f'len of elem_frag_smi_list: {len(elem_frag_smi_list)}')
# print(f'\nfrag_smi_to_2nd_break_list: \n{frag_smi_to_2nd_break_list}')
# print(f'len of frag_smi_to_2nd_break_list: {len(frag_smi_to_2nd_break_list)}')

def read_txt_to_2nd_break_dict(txt_file):
    with open(txt_file, "r") as f:
        for line in f.readlines()[:109]:
            if ':' in line:
                frag_smi_2nd_break_dict[line.split(':')[0].strip()] = [ s.strip() for s in line.split(':')[-1].strip().split(',')]
    return frag_smi_2nd_break_dict

frag_smi_2nd_break_dict = read_txt_to_2nd_break_dict("../../../frags/my_frags_3.txt")
# print(f'\nfrag_smi_2nd_break_dict: \n{frag_smi_2nd_break_dict}')
# print(f'len of frag_smi_2nd_break_dict: {len(frag_smi_2nd_break_dict)}')

# 创建分子片段到索引的映射
train_fragment_to_idx = {fragment: idx for idx, fragment in enumerate(elem_frag_smi_list)}
train_idx_to_fragment = {idx: fragment for fragment, idx in train_fragment_to_idx.items()}
train_idx_to_fragment.update({66: 'EOS'})

with open('../../../frags/mol_in_chem_space.txt', 'r') as f:
    for line in f.readlines():
        datasmi_mol_elem.append(ast.literal_eval(line.split(':')[2].strip()))
        

with open('../../../frags/mol_not_in_chem_space.txt', 'r') as f:
    for line in f.readlines():        
        # 将字符串按 ':' 分割为列表
        parts = line.split(':')
        # 去掉前两个字段，保留后面的部分，并重新组合为一个字符串
        result = ':'.join(parts[2:])        
        # word_freq_vector = word_freq_func(ast.literal_eval(result))
        # data_mol_2nd.append(word_freq_vector)
        li = []
        for i in result.strip()[1:-1].split(','):
            li.append(i.strip()[1:-1])
        datasmi_mol_2nd.append(li)

datasmi_mol_total = datasmi_mol_elem + datasmi_mol_2nd  

# print(datasmi_mol_total)
# print(len(datasmi_mol_total))
plot_list_smi = []

for lis in datasmi_mol_total:
    skip = False
    li_new = []
    for fg_smi in lis:
        if 'xxxx' in fg_smi:
            skip = True
            break
        else:
            li_new.append(fg_smi)
    if skip == True:
        continue
    else:
        plot_list_smi.append(li_new)

# print(plot_list)
# print(len(plot_list))
        
plot_list_id = []
for li_smi in plot_list_smi:
    flg = False
    li_id = []
    for fg_ in li_smi:
        if fg_ in mcts_noiso_fragment_to_idx:
            li_id.append(mcts_noiso_fragment_to_idx[fg_])
        else:
            flg = True
            break
    if flg == True:
        continue
    else:
        plot_list_id.append(li_id)

print(len(plot_list_id))


MCTS_chain = []
def extract_score_line(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for i in range(len(lines)):
            if '>  <Choices>  (1)' in lines[i]:
                if i + 1 < len(lines):
                    lx = lines[i+1]
                    MCTS_chain.append(lx)
                break

def process_sdf_files_in_directory(directory):
    for filename in os.listdir(directory):
        if filename.endswith(".sdf"):
            file_path = os.path.join(directory, filename)
            extract_score_line(file_path)


def txt_2_list(txt_path):
    out_list = []
    with open(txt_path, 'r') as file:
        for line in file:
            line = line.strip()
            if not line:
                continue
            try:
                # 将字符串转换为Python列表
                data = ast.literal_eval(line)
                # 提取子列表的第一个元素
                extracted = [item[0].split('-')[0] for item in data if isinstance(item, list)]
                out_list.append(extracted)
            except Exception as e:
                print(f"解析错误（跳过此行）: {line}")
    # print(out_list)
    return out_list


def plot_MCTS_chain_list(chain_list):
    # 示例数据
    chains = chain_list

    # 统计边的连接次数
    edge_counts = defaultdict(int)
    for chain in chains:
        for i in range(len(chain) - 1):
            edge = tuple(sorted([chain[i], chain[i + 1]]))  # 边按顺序排序，避免重复统计
            edge_counts[edge] += 1

    # 创建图
    G = nx.Graph()

    # 添加节点和边
    for edge, count in edge_counts.items():
        G.add_edge(edge[0], edge[1], weight=count)

    # 设置边的粗细
    edge_widths = [G[u][v]['weight'] for u, v in G.edges()]

    # 绘制图
    pos = nx.spring_layout(G)  # 布局算法
    nx.draw_networkx_nodes(G, pos, node_size=500, node_color='cyan')  # 绘制节点
    nx.draw_networkx_edges(G, pos, width=edge_widths, edge_color='purple', alpha=0.5)  # 绘制边
    nx.draw_networkx_labels(G, pos, font_size=12, font_color='black')  # 绘制节点标签

    # 显示图
    plt.title('MCTS Chain Visualization\n\n')
    plt.axis('off')  # 关闭坐标轴
    plt.show()


import matplotlib.pyplot as plt
from collections import defaultdict
import numpy as np

def plot_MCTS_chain_list_inline(chain_list):
    # 统计边的连接次数
    edge_counts = defaultdict(int)
    node_counts = defaultdict(int)
    for chain in chain_list:
        for i in range(len(chain) - 1):
            edge = tuple(sorted([chain[i], chain[i + 1]]))  # 边按顺序排序，避免重复统计
            edge_counts[edge] += 1
            node_counts[chain[i]] += 1

    # 获取所有节点
    nodes = set()
    for edge in edge_counts:
        nodes.add(edge[0])
        nodes.add(edge[1])
    nodes = sorted(nodes, key=lambda x: int(x))  # 按 id 排序

    # 计算节点位置
    pos = {node: (i, 0) for i, node in enumerate(nodes)}  # 所有节点在 y=0 的水平线上，x 坐标按顺序排列

    # 创建绘图
    fig, ax = plt.subplots()

    # 绘制节点
    for node, (x, y) in pos.items():
        ax.scatter(x, y, s=int(0.333*2*node_counts[node]), c='cyan', edgecolors='blue')  # 绘制节点
        ax.text(x, y-0.2, node, fontsize=12, ha='center', va='bottom')  # 绘制节点标签

    # 绘制边
    for (u, v), count in edge_counts.items():
        u_pos = pos[u]
        v_pos = pos[v]
        u_id = int(u)
        v_id = int(v)
        if abs(u_id - v_id) == 1:  # 如果节点 id 相差 1，画直线
            ax.plot([u_pos[0], v_pos[0]], [u_pos[1], v_pos[1]], color='blue', alpha=0.5, linewidth=0.2*count/30)
        else:  # 否则画弯曲边，跳出 y=0 的水平线
            # 计算控制点
            mid_y = 0.5*count/30  # 控制点的高度
            control_point = ((u_pos[0] + v_pos[0]) / 2, mid_y)  # 控制点在中间位置
            # 使用贝塞尔曲线绘制弯曲边
            bezier_path = [
                (u_pos[0], u_pos[1]),  # 起点
                control_point,  # 控制点
                (v_pos[0], v_pos[1])  # 终点
            ]
            t = np.linspace(0, 1, 100)
            curve = np.outer((1 - t) ** 2, bezier_path[0]) + \
                    np.outer(2 * (1 - t) * t, bezier_path[1]) + \
                    np.outer(t ** 2, bezier_path[2])
            ax.plot(curve[:, 0], curve[:, 1], color='blue', alpha=0.5, linewidth=0.2*count/30)

    # 显示图
    plt.subplots_adjust(left=0.0001, right=0.9999)  # 调整左右边距
    ax.set_title('MCTS Chain in DrugBank Visualization\n\n')
    ax.set_axis_off()  # 关闭坐标轴
    plt.show()


if __name__ == "__main__":
    work_directory = './dual_AB_QWEN_500'
    process_sdf_files_in_directory(work_directory)
    with open('./dual_AB_QWEN_500_MCTS_chain.txt', "w") as fw:
        for i in MCTS_chain:
            fw.write(i)

    chain_list = txt_2_list('./dual_AB_QWEN_500_MCTS_chain.txt')
    # plot_MCTS_chain_list(chain_list)
    plot_MCTS_chain_list_inline(plot_list_id)
    
    

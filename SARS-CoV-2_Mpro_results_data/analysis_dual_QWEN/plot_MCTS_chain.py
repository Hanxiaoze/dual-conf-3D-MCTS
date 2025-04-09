import os
import shutil
from openbabel import pybel
import ast
import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict

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
    print(out_list)
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
        ax.scatter(x, y, s=int(2*node_counts[node]), c='cyan', edgecolors='blue')  # 绘制节点
        ax.text(x, y-0.2, node, fontsize=12, ha='center', va='bottom')  # 绘制节点标签

    # 绘制边
    for (u, v), count in edge_counts.items():
        u_pos = pos[u]
        v_pos = pos[v]
        u_id = int(u)
        v_id = int(v)
        if abs(u_id - v_id) == 1:  # 如果节点 id 相差 1，画直线
            ax.plot([u_pos[0], v_pos[0]], [u_pos[1], v_pos[1]], color='blue', alpha=0.5, linewidth=0.2*count)
        else:  # 否则画弯曲边，跳出 y=0 的水平线
            # 计算控制点
            mid_y = 0.5*count  # 控制点的高度
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
            ax.plot(curve[:, 0], curve[:, 1], color='blue', alpha=0.5, linewidth=0.2*count)

    # 显示图
    plt.subplots_adjust(left=0.0001, right=0.9999)  # 调整左右边距
    ax.set_title('MCTS Chain Visualization\n\n')
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
    plot_MCTS_chain_list_inline(chain_list)
    
    

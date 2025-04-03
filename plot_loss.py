import matplotlib.pyplot as plt
import pandas as pd
import os
import glob
import numpy as np

# 设置中文字体（如果标签含中文）
# plt.rcParams['font.sans-serif'] = ['SimHei']  # Windows系统        # 使用衬线字体族
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['axes.unicode_minus'] = False    # 解决负号显示问题
plt.rcParams['font.size'] = 16

# 指定CSV文件路径（例如：读取所有以 "data_" 开头的CSV文件）
csv_folder = "./"  # 替换为你的文件夹路径
csv_files = glob.glob(os.path.join(csv_folder, "Average_Loss_*.csv"))

# 存储所有数据的DataFrame列表
dataframes = []

# 读取每个CSV文件
for file in csv_files:
    df = pd.read_csv(file)
    df = df.sort_values(by="Step")  # 按x列排序（避免绘图时线条混乱）
    df = df.iloc[:150]
    dataframes.append(df)

plt.figure(figsize=(12, 6))  # 设置画布大小

# 定义颜色和线型（自动循环）
colors = plt.cm.tab10(np.linspace(0, 1, len(csv_files)))  # 使用tab10色系
line_styles = ['-', '--', '-.', ':']  # 线型循环

lable_list = ['Training', 'Validation']

for idx, (df, file) in enumerate(zip(dataframes, csv_files)):
    # 提取文件名作为标签（去除路径和后缀）
    # label = os.path.splitext(os.path.basename(file))[0]
    label = lable_list[idx]
    
    # 绘制线图
    plt.plot(
        df["Step"], 
        df["Value"], 
        color=colors[idx],          # 自动分配颜色
        linestyle=line_styles[idx % len(line_styles)],  # 循环线型
        linewidth=2,
        marker='o',                # 数据点标记
        markersize=5,
        label=label                 # 图例标签
    )

# 图表装饰
# plt.title("多文件数据对比图", fontsize=14)
plt.ylabel("Average Loss")
plt.xlabel("Epoch")
# plt.grid(True, linestyle='--', alpha=0.6)  # 网格线
plt.legend(loc='upper right') # 图例位置

# 自动调整刻度密度
# plt.xticks(rotation=45)  # 如果x轴标签过长，旋转45度

# 紧凑布局
plt.tight_layout()

# 显示图表
plt.show()

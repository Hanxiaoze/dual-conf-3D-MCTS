import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rcParams
# 设置全局字体为Times
rcParams['font.family'] = 'Times New Roman'
# 设置全局字体大小
plt.rc('font', size=16)

with open('./dual_conform_affinity.txt', 'r') as file:
    lines = file.readlines()


x_vals = []
y_vals = []

for line in lines:
    parts = line.split()
    x_vals.append(float(parts[-2]))
    y_vals.append(float(parts[-1]))

# 创建2D直方图
plt.figure(figsize=(10, 8))
plt.hist2d(x_vals, y_vals, bins=(40, 40), cmap=plt.cm.jet)
plt.colorbar(label='\nFrequency')
plt.xlabel('\nConformation 1 affinity')
plt.ylabel('Conformation 2 affinity\n')
plt.title('2D Frequency Histogram\n')
plt.show()

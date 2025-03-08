import numpy as np
import matplotlib.pyplot as plt

# 将数据解析为数组
lines = data.strip().split('\n')
x_vals = []
y_vals = []

for line in lines:
    parts = line.split()
    x_vals.append(float(parts[-2]))
    y_vals.append(float(parts[-1]))

# 创建2D直方图
plt.figure(figsize=(10, 8))
plt.hist2d(x_vals, y_vals, bins=(20, 20), cmap=plt.cm.jet)
plt.colorbar(label='Frequency')
plt.xlabel('Value 1')
plt.ylabel('Value 2')
plt.title('2D Frequency Histogram')
plt.show()

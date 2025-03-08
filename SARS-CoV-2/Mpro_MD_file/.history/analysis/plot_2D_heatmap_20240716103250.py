import numpy as np
import matplotlib.pyplot as plt



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
plt.colorbar(label='Frequency')
plt.xlabel('Conformation 1 affinity')
plt.ylabel('Conformation 1 affinity')
plt.title('2D Frequency Histogram')
plt.show()

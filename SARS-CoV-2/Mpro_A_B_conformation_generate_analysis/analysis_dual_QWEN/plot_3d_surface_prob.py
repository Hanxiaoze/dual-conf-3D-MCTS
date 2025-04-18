import sys

# arg1 = str(sys.argv[1])
# arg2 = str(sys.argv[2])

import matplotlib.pyplot as plt
import numpy as np

from matplotlib import rcParams
# 设置全局字体为Times
rcParams['font.family'] = 'Times New Roman'
# 设置全局字体大小
plt.rc('font', size=16)


# 从文件中读取数据


file_path_2 = './dual_conform_affinity_5000.txt'

# data = np.loadtxt(file_path, skiprows=0, usecols=(1, 2))
# x = data[:, 0]
# y = data[:, 1]

# # 手动设置 x 和 y 的范围
x_range = (-10, -2)
y_range = (-10, -2)

# # 网格离散化
bins = 50
# hist, x_edges, y_edges = np.histogram2d(x, y, bins=bins, density=True, range=[x_range, y_range])

# # 创建3D图形
fig = plt.figure()
# ax = fig.add_subplot(221, projection='3d')

# # 获取网格的中心点坐标
# x_centers = (x_edges[:-1] + x_edges[1:]) / 2
# y_centers = (y_edges[:-1] + y_edges[1:]) / 2

# # 生成网格坐标
# x_grid, y_grid = np.meshgrid(x_centers, y_centers)


# # 使用scatter绘制3D散点图，高度信息为频数
# # ax.scatter(x_grid, y_grid, hist.T, c=hist.T, cmap='Reds', marker='o', s=50, alpha=0.8)    # cmap='coolwarm'

# # 使用plot_surface绘制平滑曲面图
# surf = ax.plot_surface(x_grid, y_grid, hist.T, cmap='Reds', rstride=1, cstride=1, antialiased=True, alpha=0.8, vmin=0, vmax=0.5)



# # 添加标题和标签
# ax.set_title('Joint probability density of the traning molecules', fontsize=12)
# ax.set_xlabel('Epro binding affinity  (kcal/mol)', fontsize=10)
# ax.set_ylabel('Mpro binding affinity  (kcal/mol)', fontsize=10)
# ax.set_zlabel('Probability density', fontsize=10)

# # 在3D图形底面投影二维轮廓
# # contour = ax.contour(x_grid, y_grid, hist.T, zdir='z', offset=-600, cmap='coolwarm', linewidths=1)
# # 在3D图形底面投影二维轮廓并填充轮廓
# contour = ax.contourf(x_grid, y_grid, hist.T, zdir='z', offset=-0.5, cmap='Reds', levels=20, alpha=0.8, vmin=0, vmax=0.5)

# # 添加颜色条，调整 Colorbar 的大小
# cbar1 = fig.colorbar(surf, ax=ax, label=' ', shrink=0.5)


# # 设置 x 轴的绘制范围为 5 到 10
# ax.set_xlim(-14, -2)
# ax.set_ylim(-14, -2)
# ax.set_zlim(-0.5, 0.5)
# ax.view_init(elev=12, azim=255)


# # 在3D图形右方绘制二维投影热图
# ax2 = fig.add_subplot(222)
# # im = ax2.imshow(hist.T, cmap='coolwarm', extent=[x_centers.min(), x_centers.max(), y_centers.min(), y_centers.max()])
# ax2.hist2d(x, y, bins=bins, density=True, cmap='Reds', range=[x_range, y_range], vmin=0, vmax=0.5)
# ax2.set_title('2D projection probability density of the traning molecules', fontsize=12)
# ax2.set_xlabel('Epro binding affinity  (kcal/mol)', fontsize=10)
# ax2.set_ylabel('Mpro binding affinity  (kcal/mol)', fontsize=10)
# # 添加颜色条，调整 Colorbar 的大小
# cbar2 = fig.colorbar(surf, ax=ax2, label=' ', shrink=0.5)

# # 设置 x 轴的绘制范围为 5 到 10
# ax2.set_xlim(-14, -2)
# ax2.set_ylim(-14, -2)

##########
##########
x_vals = []
y_vals = []
with open(file_path_2, 'r') as file:
    for line in file:
        parts = line.split()
        x_vals.append(float(parts[-2]))
        y_vals.append(float(parts[-1]))

# 转换为numpy数组
x_vals = np.array(x_vals)
y_vals = np.array(y_vals)



x_2 = x_vals
y_2 = y_vals

# 网格离散化
hist_2, x_edges_2, y_edges_2 = np.histogram2d(x_2, y_2, bins=bins, density=True, range=[x_range, y_range])

# 创建3D图形
ax_3 = fig.add_subplot(111, projection='3d')

# 获取网格的中心点坐标
x_centers_2 = (x_edges_2[:-1] + x_edges_2[1:]) / 2
y_centers_2 = (y_edges_2[:-1] + y_edges_2[1:]) / 2

# 生成网格坐标
x_grid_2, y_grid_2 = np.meshgrid(x_centers_2, y_centers_2)


# 使用scatter绘制3D散点图，高度信息为频数
# ax.scatter(x_grid, y_grid, hist.T, c=hist.T, cmap='Reds', marker='o', s=50, alpha=0.8)    # cmap='coolwarm'

# 使用plot_surface绘制平滑曲面图
surf_2 = ax_3.plot_surface(x_grid_2, y_grid_2, hist_2.T, cmap='coolwarm', rstride=1, cstride=1, antialiased=True, alpha=0.8, vmin=0, vmax=2.0)



# 添加标题和标签
# ax_3.set_title('Generated molecules with affinity conditions Epro= -{}.0 / Mpro= -{}.0 kcal/mol'.format(arg1, arg2))
ax_3.set_title('Dual conformation generated molecules binding affinities distribution')
ax_3.set_xlabel('\nConf_1 binding affinity  (kcal/mol)')
ax_3.set_ylabel('\nConf_2 binding affinity  (kcal/mol)')
ax_3.set_zlabel('\nProbability density')

# 在3D图形底面投影二维轮廓
# contour = ax.contour(x_grid, y_grid, hist.T, zdir='z', offset=-600, cmap='coolwarm', linewidths=1)
# 在3D图形底面投影二维轮廓并填充轮廓
contour_2 = ax_3.contourf(x_grid_2, y_grid_2, hist_2.T, zdir='z', offset=-1.1, cmap='coolwarm', levels=20, alpha=0.8, vmin=0, vmax=2.0)

# 添加颜色条，调整 Colorbar 的大小
cbar3 = fig.colorbar(surf_2, ax=ax_3, label=' ', shrink=0.5)

# 设置 x 轴的绘制范围为 5 到 10
ax_3.set_xlim(-10, -2)
ax_3.set_ylim(-10, -2)
ax_3.set_zlim(-1.0, 2.0)
ax_3.view_init(elev=12, azim=235)

# # 在3D图形右方绘制二维投影热图
# ax_4 = fig.add_subplot(224)
# # im = ax_4.imshow(hist_2.T, cmap='coolwarm', extent=[x_centers_2.min(), x_centers_2.max(), y_centers_2.min(), y_centers_2.max()])
# ax_4.hist2d(x_2, y_2, bins=bins, density=True, cmap='Reds', range=[x_range, y_range], vmin=0, vmax=0.5)
# ax_4.set_title('2D projection probability density of the generated molecules', fontsize=12)
# ax_4.set_xlabel('Epro binding affinity  (kcal/mol)', fontsize=10)
# ax_4.set_ylabel('Mpro binding affinity  (kcal/mol)', fontsize=10)
# # 添加颜色条，调整 Colorbar 的大小
# cbar4 = fig.colorbar(surf_2, ax=ax_4, label=' ', shrink=0.5)

# # 设置 x 轴的绘制范围为 5 到 10
# ax_4.set_xlim(-14, -2)
# ax_4.set_ylim(-14, -2)


# 调整布局，紧凑排列子图
plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.03, hspace=0.03)
# plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.3, hspace=0.3)

# 调整布局
plt.tight_layout()


# 显示图形
plt.show()

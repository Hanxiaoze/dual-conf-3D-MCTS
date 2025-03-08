import matplotlib.pyplot as plt
import numpy as np

from matplotlib import rcParams
# 设置全局字体为Times
rcParams['font.family'] = 'Times New Roman'

# 设置全局字体大小
plt.rc('font', size=20)

# 从文件中读取数据

file_path = './E_M_6W_affinity.txt'
file_path_2 = './E_M_gen_affinity_warm_up.txt'

data = np.loadtxt(file_path, skiprows=0, usecols=(1, 2))
x = data[:, 0]
y = data[:, 1]

data_2 = np.loadtxt(file_path_2, skiprows=0, usecols=(1, 2))
x_2 = data_2[:, 0]
y_2 = data_2[:, 1]


file_path = './dual_conform_affinity_5000.txt'

file_path_2 = '/home/zhou/3D-MCTS-main/nsp5_conf_1/run_1/validate_results.txt'

file_path_3 = '/home/zhou/3D-MCTS-main/nsp5_conf_2/run_1/validate_results.txt'






# 手动设置 x 和 y 的范围
x_range = (-10, -2)
y_range = (-10, -2)

# 网格离散化
bins = 50
hist, x_edges, y_edges = np.histogram2d(x, y, bins=bins, density=True, range=[x_range, y_range])
hist_2, x_edges_2, y_edges_2 = np.histogram2d(x_2, y_2, bins=bins, density=True, range=[x_range, y_range])

# 创建3D图形
fig = plt.figure()
ax = fig.add_subplot(121, projection='3d')

# 获取网格的中心点坐标
x_centers = (x_edges[:-1] + x_edges[1:]) / 2
y_centers = (y_edges[:-1] + y_edges[1:]) / 2
# 生成网格坐标
x_grid, y_grid = np.meshgrid(x_centers, y_centers)
# 使用plot_surface绘制平滑曲面图
surf = ax.plot_surface(x_grid, y_grid, hist.T, color='orange', rstride=1, cstride=1, antialiased=True, alpha=0.6)


# 获取网格的中心点坐标
x_centers_2 = (x_edges_2[:-1] + x_edges_2[1:]) / 2
y_centers_2 = (y_edges_2[:-1] + y_edges_2[1:]) / 2
# 生成网格坐标
x_grid_2, y_grid_2 = np.meshgrid(x_centers_2, y_centers_2)
# 使用plot_surface绘制平滑曲面图
surf_2 = ax.plot_surface(x_grid_2, y_grid_2, hist_2.T, color='cyan', rstride=1, cstride=1, antialiased=True, alpha=0.2)


# 添加标题和标签
ax.set_title('Joint probability density of \nthe traning molecules (orange) and generated molecules (cyan)')
ax.set_xlabel('\n\nEpro binding affinity  (kcal/mol)')
ax.set_ylabel('\n\nMpro binding affinity  (kcal/mol)')
ax.set_zlabel('\n\nProbability density')


# 设置 x 轴的绘制范围为 5 到 10
ax.set_xlim(-14, -2)
ax.set_ylim(-14, -2)
ax.set_zlim(0, 0.7)
ax.view_init(elev=3, azim=255)


##########
##########

ax_2 = fig.add_subplot(122, projection='3d')


# 创建掩码，将不满足条件的值设为 np.ma.masked
masked_data = np.ma.masked_where(hist_2.T - hist.T < 0, hist_2.T - hist.T)
masked_data_1 = np.ma.masked_where(hist.T - hist_2.T < 0, hist_2.T - hist.T)
masked_data_2 = np.ma.masked_where(hist.T - hist_2.T == 0, hist_2.T - hist.T)

surf_3 = ax_2.plot_surface(x_grid, y_grid, masked_data, color='red', rstride=1, cstride=1, antialiased=True, alpha=0.6)
surf_4 = ax_2.plot_surface(x_grid, y_grid, masked_data_1, color='green', rstride=1, cstride=1, antialiased=True, alpha=0.5)
surf_5 = ax_2.plot_surface(x_grid, y_grid, masked_data_2, color='red', rstride=1, cstride=1, antialiased=True, alpha=0.6)
# surf_3.set_facecolors(np.where((hist_2 - hist).T > 0, 'red', 'green'))

# 添加标题和标签
ax_2.set_title('Joint probability density difference of \nthe generated molecules - training molecules (Red: +, Green: -)')
ax_2.set_xlabel('\n\nEpro binding affinity  (kcal/mol)')
ax_2.set_ylabel('\n\nMpro binding affinity  (kcal/mol)')
ax_2.set_zlabel('\n\nProbability density difference')

ax_2.set_xlim(-14, -2)
ax_2.set_ylim(-14, -2)
ax_2.set_zlim(-0.5, 0.5)
ax_2.view_init(elev=3, azim=255)






# 调整布局
plt.tight_layout()


# 显示图形
plt.show()

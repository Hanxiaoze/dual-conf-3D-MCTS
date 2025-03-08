import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.ndimage import gaussian_filter
from matplotlib import rcParams

# 设置全局字体为Times
rcParams['font.family'] = 'Times New Roman'
f_size = 25
# 设置全局字体大小
plt.rc('font', size=f_size, weight='bold')

# 从文件中读取数据
filename = 'hists_1-2.gnu'  # 假设数据文件名为 hists_1-2.gnu
tmp_filename = 'tmp.txt'

# 使用 loadtxt 读取数据，并跳过第2-5行和第15-20行
def load_filtered_data(filename):
    filtered_lines = []
    with open(filename, 'r') as f:
        lines = f.readlines()

    max_lines = len(lines)

    for i, line in enumerate(lines):
        # 过滤掉第2-5行（索引1-4）和第15-20行（索引14-19）
        if 1 <= i+1 <= 7 or max_lines-3 <= i+1 <= max_lines:
            continue
        filtered_lines.append(line)

    # 将过滤后的行写入临时文件
    with open(tmp_filename, 'w') as tmp_file:
        tmp_file.writelines(filtered_lines)

    # 使用 np.loadtxt 从临时文件读取数据
    return np.loadtxt(tmp_filename)

# 读取并过滤数据
data_filtered = load_filtered_data(filename)

data = data_filtered


# 获取 x, y, z 列数据
x = data[0:, 0]
y = data[0:, 1]
z = data[0:, 2]

# 创建网格数据（用于绘制曲面图）
x_grid = np.unique(x)
y_grid = np.unique(y)
X, Y = np.meshgrid(x_grid, y_grid)

# 插值 Z 值
Z = np.zeros_like(X)

scale = 1

# 填充 Z 值（对应于 X, Y 坐标的 Z 值）
for i in range(len(x)):
    ix = np.where(x_grid == x[i])[0][0]
    iy = np.where(y_grid == y[i])[0][0]
    Z[iy, ix] = z[i] * scale - np.max(z) * scale

# 高斯平滑
sigma = 1.5  # 设置高斯平滑的标准差
Z_smooth = gaussian_filter(Z, sigma=sigma)

# 绘制 3D 曲面图
fig = plt.figure(figsize=(15, 15))
ax = fig.add_subplot(111, projection='3d')

# 绘制平滑后的曲面
surf = ax.plot_surface(X, Y, Z_smooth, cmap='viridis', alpha=0.5)

# 添加颜色条，并调整位置
cbar = fig.colorbar(surf, ax=ax, fraction=0.02, pad=0.1)  # fraction 控制颜色条的宽度，pad 控制颜色条和图形的间距

# 设置标签
ax.set_xlabel('\n\nPC1', fontweight='bold')
ax.set_ylabel('\n\nPC2', fontweight='bold')
ax.set_zlabel('\n\n\nFree Energy Surface (Kcal/mol)', fontweight='bold')

# 2D颜色投影图：在(x, y, 0)平面上绘制
# 创建一个与 X, Y 匹配的平面，并将其 Z 值设置为 0
Z_zero = np.zeros_like(X) - np.max(z) * scale

# 绘制平面投影，将其显示在 3D 图的底部
ax.plot_surface(X, Y, Z_zero, facecolors=surf.to_rgba(Z_smooth), rstride=1, cstride=1, alpha=0.5)

ax.view_init(elev=33, azim=-56)

# 保存图形
plt.savefig('./FEL_with_projection_on_plane_smoothed.png', dpi=600, bbox_inches='tight')

# 显示图形
plt.show()

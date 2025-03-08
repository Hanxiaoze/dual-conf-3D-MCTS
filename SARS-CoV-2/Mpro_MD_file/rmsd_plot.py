import pandas as pd
import matplotlib.pyplot as plt

from matplotlib import rcParams
# 设置全局字体为Times
rcParams['font.family'] = 'Times New Roman'


# 读取文本文件
df = pd.read_csv('./rmsd_ref_C0_pocket.dat', delim_whitespace=True, skiprows=1, header=None, names=['Frame', 'RMSD'])

# 绘制折线图
plt.figure(figsize=(10, 6))
plt.plot(df['Frame']/1, df['RMSD'], marker='', color='r', linestyle='-', label='Mpro')

# 添加标题和标签
# plt.title('RMSD vs Frame')
plt.xlabel('Time (ns)')
plt.ylabel('RMSD (\u212B)')

#plt.ylim(0, 10) 

#plt.xticks([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100])

plt.legend()

# 显示图形
# plt.grid(True)
plt.show()

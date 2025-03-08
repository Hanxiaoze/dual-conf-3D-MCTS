import numpy as np

# 从文件中读取数据
filename = 'projection.dat'  # 假设数据文件名为 data.txt

# 使用 numpy 读取文件，跳过标题行
data = np.loadtxt(filename, skiprows=1)

# 输入 Mode1 和 Mode2 的范围
mode1_min, mode1_max = map(float, input("请输入 Mode1 的最小值和最大值，用空格分隔: ").split())
mode2_min, mode2_max = map(float, input("请输入 Mode2 的最小值和最大值，用空格分隔: ").split())

# 过滤满足条件的行
filtered_data = data[(data[:, 1] >= mode1_min) & (data[:, 1] <= mode1_max) &
                     (data[:, 2] >= mode2_min) & (data[:, 2] <= mode2_max)]

# 打印满足条件的行
if filtered_data.size > 0:
    print(f"满足条件的行信息: {mode1_min} < model1 < {mode1_max}, {mode2_min} < model2 < {mode2_max} ")
    print("Frame    Mode1    Mode2    Mode3")
    for row in filtered_data:
        print(f"{int(row[0]):<8} {row[1]:<8.3f} {row[2]:<8.3f} {row[3]:<8.3f}")
else:
    print("没有满足条件的数据")


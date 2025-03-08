# 打开原始文件和新文件
with open('/home/zhou/3D-MCTS-main/SARS-CoV-2/8A4Y.pdb', 'r') as input_file, open('/home/zhou/3D-MCTS-main/SARS-CoV-2/8A4Y_noABres.pdb', 'w') as output_file:
    # 逐行读取原始文件
    lines = input_file.readlines()

    # 遍历每一行
    for i, line in enumerate(lines):
        if len(line) >= 17:        
            if line[16] is not "A" and line[16] is not "B":
                output_file.write(line)
            if line[16] is "A":
                line = line[0 : 16] + " " + line[17 : ]
                output_file.write(line)
        else:
            output_file.write(line)
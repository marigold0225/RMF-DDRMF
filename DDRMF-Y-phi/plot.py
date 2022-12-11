# *************************************************************************
#      > File Name: plot.py
#      > Author: marigold
#      > Mail: mflovelky418@gmail.com
#      > Created Time: 2022年01月07日 星期五 21时49分56秒
# ************************************************************************
#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

file = './output/M_R.txt'
a = np.loadtxt(file)
x = a[:, 3]  # 取第一列数据
y = a[:, 2]  # 取第二列数据
# 画图
plt.plot(x, y)
plt.legend(['DDME2', 'Cubic-Spline'])
plt.xlabel('R(km)')
plt.ylabel('M(M_sun)')
plt.title('M-R')
# 保存图片
plt.savefig('./output/plot/M-R.jpg')
plt.xlim(xmin=10, xmax=18)
plt.ylim(ymin=0, ymax=2.5)
plt.show()

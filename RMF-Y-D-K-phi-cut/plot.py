#*************************************************************************
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
x = a[:,3]  # 取第一列数据
y = a[:,2]  # 取第二列数据
# 画图
plt.plot(x,y)
plt.xlim(11,15)
plt.ylim(0,2.5)
plt.legend(['FSUGOLD','Cubic-Spline'])
plt.xlabel('R(km)') 
plt.ylabel('M(M_sun)')
plt.title('M-R')  
# 保存图片  
plt.savefig('./output/M-R.jpg')
plt.xlim(xmin = 11 , xmax = 15)
plt.ylim(ymin = 0 , ymax = 3)
plt.show()

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
R = a[:,3]
M = a[:,2]
# 画图
#plt.figure()
#ax=plt.gca()
#ax.set_title(r"NS Plot: Love Number vs Mass")
#ax.set_xlabel(r"Mass of the Star [$M_\odot$]")
#ax.set_ylabel(r"Love Number [$k_2$]")
#plt.plot(M,k2,'r--',label="DD2 npeuYD")
##plt.plot(M_f,k2_f,'b--',label="DD2 -120")
##plt.plot(M_a,k2_a,'g--',label="DD2 -140")
#plt.legend(loc="best")
#plt.savefig("k2vM_TD.png")
#
plt.figure()
ax=plt.gca()
ax.set_title(r"M-R relationship")
ax.set_xlabel(r"Radius[$km$]")
ax.set_ylabel(r"M/M$sun$")
ax.set_xlim(8,15)
ax.set_ylim(0,3)
plt.plot(R,M,'r--',label="M-R")
#plt.plot(C_f,k2_f,'b--',label="DD2 -120")
#plt.plot(C_a,k2_a,'g--',label="DD2 -140")
plt.legend(loc="best")
plt.savefig("M_R.png")


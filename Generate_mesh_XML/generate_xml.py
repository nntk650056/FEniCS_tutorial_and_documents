"""
@author: 章肖阳
"""
'''****************************第一部分*******************************'''
import numpy as np
import pandas as pd
from sklearn.model_selection import ShuffleSplit
import csv
import matplotlib.pyplot as plt
from pandas import DataFrame
from scipy.linalg import solve
from scipy.interpolate import griddata
import matplotlib as mpl

n_nodes = 3759
n_elements = 14860

# read data
f=open('3d_gel_0_5.txt','r')
n=0
#Data = f.readlines()  #直接将文件中按行读到list里，效果与方法2一样
Data={} #记录文件所有行数据信息的字典，检索方式为行数
Num={}  #记录关键字所在行数的字典，检索方式为关键字字符串
for line in f:
    Data[n]=l = list(map(str, line.split()))   #将整个line 分隔 
    n=n+1
    
#定义数组储存应变应力信息
linenum=0

nodes = np.zeros(shape=(n_nodes,4))
elements = np.zeros(shape=(n_elements,5),dtype = int)


#从文件中获取最大主应变、最小主应变以及von-Mises应力            
for i in range(n_nodes):
    for j in range(4):
        #print(j)
        nodes[i,j] = float(Data[i+1][j])
print(Data[n_nodes])        
for i in range(n_elements):
    for j in range(5):
        #print(j)
        elements[i,j] = int(Data[n_nodes+2+i][j]) - 1

f.close()

with open("E:\\3D_gel_size0_5.xml","w") as f:
    
    f.write('<?xml version="1.0" encoding="UTF-8"?>')   #设置文件对象
    f.write('\n')
    f.write('\n')
    f.write('<dolfin xmlns:dolfin="http://fenicsproject.org">')
    f.write('\n')
    
    f.write('  <mesh celltype="tetrahedron" dim="3">')  #设置文件对象
    f.write('\n')
    f.write('    <vertices size="')
    f.write(str(n_nodes))
    f.write('">')
    f.write('\n')

    for i in range(n_nodes):   #对于双层列表中的数据
        nnn = str('      ') + str('<vertex index="') + str(i) + str('" x="') + str(nodes[i,1]) + str('" y="') + str(nodes[i,2]) + str('" z="') + str(nodes[i,3]) +str('" />')
        f.write(nnn)
        f.write('\n')

    f.write('    </vertices>')
    f.write('\n')
    f.write('    <cells size="')
    f.write(str(n_elements))
    f.write('">')
    f.write('\n')
    for i in range(n_elements):
        nnn = str('      ') + str('<tetrahedron index="') + str(i) + str('" v0="') + str(elements[i,1]) + str('" v1="') + str(elements[i,2]) + str('" v2="') + str(elements[i,3]) + str('" v3="') + str(elements[i,4]) + str('" />')                                                                #对于双层列表中的数据
        f.write(nnn)
        f.write('\n')
    f.write('    </cells>')
    f.write('\n')
    f.write('  </mesh>')
    f.write('\n')
    f.write('</dolfin>')
    f.write('\n')
        
f.close

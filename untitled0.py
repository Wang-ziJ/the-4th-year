# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 08:12:28 2020

@author: 王子健
"""

import math 
import matplotlib.pyplot as plt

#质量
m_1 = 1
m_2 = 2
m_3 = 3
#距离
x_2 = 5
x_3 = 10
length = 5    #弹簧长度
#劲度系数
k = 10
#初始速度
v_0 = input('输入一个整数：')
print (v_0)
v_0 = int(v_0)
#初始位移
x_1 = input('输入一个整数：')
print (x_1)
x_1 = int(x_1)
#微分
#时间
t = 5
n = 5000
dt = t / n

X_1 = []
X_2 = []
X_3 = []
T = []

v_1 = v_0
v_2 = 0
v_3 = 0

#对物块m_1
'''for i in range(0,501):
    T.append(T_0 + i * dt)
    dx = 0
    if dx == 0:
        a = -k * (x_0 + dx) / m_1
        v = v_0 + a * dt
        x = x_0 + v * dt
        X.append(x)
        dx = x - x_0
    else:
        a = -k * (x_0 + dx) / m_1
        v = v_0 + a * dt
        x = x_0 + v * dt
        X.append(x)
        dx = x - x_0'''
        
for i in range(0,5001):
    T.append(i * dt)
    a_1 = -k * (length - (x_2 - x_1)) / m_1    #a_1 物块一的加速度
    v_1 = v_1 + a_1 * dt
    x_1 = x_1 + v_1 * dt
    X_1.append(x_1)
    a_2 = k * ((length - (x_2 - x_1)) - (length - (x_3 - x_2))) / m_2
    v_2 = v_2 + a_2 * dt
    x_2 = x_2 + v_2 * dt
    X_2.append(x_2)
    a_3 = k * (length - (x_3 - x_2)) / m_3    #a_3 物块二的加速度
    v_3 = v_3 + a_3 * dt
    x_3 = x_3 + v_3 * dt
    X_3.append(x_3)

plt.scatter(T,X_1,s = 10)
plt.show
plt.scatter(T,X_2,s = 10)
plt.show
plt.scatter(T,X_3,s = 10)
plt.show      
    
    
    
    
# -*- coding: utf-8 -*-
import math 
import matplotlib.pyplot as plt

x1 = -math.pi 
x2 = math.pi
num = 1000
dx = (x2 - x1) / num
h = dx * 10
x = []
y = []
for m in range(0,1001):
    x.append(x1 + m * dx)
for j in range(0,1001):
    delta = 0
    a = 0
    b = 0
    c = 0
    for i in range(j-10,j+11):
        if i - j == 0:
            c = 0
        else:
            R = abs((i - j) * dx) / h
            a = (1 - R) ** 2
            b = -12 * R
            c = ((i - j) * dx) / (h * abs((i - j) * dx))
            delta = delta - 1.25 * a * b * c * math.sin(x1 + i * dx) * dx / h 
    y.append(delta)
plt.scatter(x,y,s = 10)
plt.show


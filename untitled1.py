# -*- coding: utf-8 -*-
import math 
#from mpl_toolkits.mplot3d import axes3d
#from matplotlib.widgets import Slider, Cursor
import matplotlib.pyplot as plt
import numpy as np
#from glasstone.overpressure import brode_overpressure, soviet_overpressure
"""
Created on Sun Sep 13 03:18:42 2020

@author: 王子健
"""
r=int(input("爆炸距离："))
Q=int(input("当量："))
omega=0
b=0
beta=0
t1=1                                   #冲击波到达时间
t0=1                                   #正相持续时间
dp=0                                   #超压
dpt=0

if r/pow(Q,1/3)<=76:
    dp=1.0812*10**4*pow(Q,2/3)/pow(r,2)+3.654*pow(10,6)*Q/pow(r,3)
if 76<r/pow(Q,1/3)<=860:
    dp=0.8678*10**3*pow(Q,1/3)/r+0.8351*10**4*pow(Q,2/3)/pow(r,2)+3.3*10**6*Q/pow(r,3)-0.016
if r/pow(Q,1/3)>860:
    dp=0.5102*10**2*pow(Q,1/3)/r+3.023*10**4*pow(Q,2/3)/pow(r,2)
    
c=8.71+0.184*dp-104/(10+dp)


if dp>=16:                      
    t0=0.1645*Q**(1/3)                 #正相持续时间
elif dp<16 and dp>=2:
    t0=0.1411*Q**(1/3)*dp**(0.07813*math.log(dp)-0.1604)
elif dp<2 and dp>0.1:
    t0=0.172*dp**(-0.3431)*Q**(1/3)
elif dp<=0.1 and dp>0.01:
    t0=0.2579*dp**(-0.1775)*Q**(1/3)
    

    

if dp<=4:
    beta=0.735+0.974*dp

if dp>4:
    beta=0.0759*dp**2+0.646*dp+1.286



#到达时间
if r/pow(Q,1/3)<=42:
    t1=6.45*10**(-7)*pow(r**5/Q,1/2)   
if 42<r/pow(Q,1/3)<=126:
    t1=-1.0377*10**(-4)*r+6.4026*10**(-6)*r**2/pow(Q,1/3)+3.8984*10**(-9)*r**3/pow(Q,2/3)
if 126<r/pow(Q,1/3)<=422:
    t1=8.2548*(pow(0.121141*10**(-6)*r**2+0.37757*10**(-2)*pow(Q,2/3),1/2)-0.64031/10*pow(Q,1/3))
if r/pow(Q,1/3)>422:
    t1=-0.272819*pow(Q,1/3)+0.294014/100*r-0.222095*pow(Q,1/3)*pow(math.log(r/pow(Q,1/3))-5.3685,1/2)


print(beta)
print(t1)
print(t0)
print(dp)

fig,ax1 = plt.subplots()

#动压随时间变化
t=np.linspace(0,10,1000) 
y=[]
qo=2.5*dp**2/(dp+7)  #动压
print(qo)
for i in t:
   if i<=t1:
      y.append(0)
   if i>t1:
      qt=qo*(1-(i-t1)/t0)**2*math.exp(-beta*(i-t1)/t0)
      y.append(qt)
ax1.set_title('Nuclear effect manual')
ax1.set_xlabel('arrival time($s$)')
ax1.set_ylabel('dynamic pressure($psi$)')
plt.plot(t,y)
plt.show()

fig,ax2 = plt.subplots()


#超压随时间变化
#负压峰值
dp_=0
dp_t=0
if r/pow(Q,1/3)<=110:
    dp_=0.23   
if r/pow(Q,1/3)>110:
    dp_=25.3*pow(Q,1/3)/r
print(dp_)
#核爆冲击波负压作用时间
t_0=1.069*pow(Q,1/3)



t=np.linspace(0,10,1000) 
y=[]
for i in t:
   if i<=t1:
      y.append(0)
   if i>t1:
      if dp<50: 
          if dp<=4:
              omega=0.5+0.8*dp
          elif 4<dp<50:
              if 4<dp<10:
                  b=dp*(0.88+0.072*dp)
              if 10<dp<50:
                  b=dp*(1.67-0.011*dp)
              omega=b/(1+c*(i-t1)/t0)
          dp_t=dp_*14*(i-t1)/t_0*(1-(i-t1)/t_0)*math.exp(-4*(i-t1)/t_0)
          dpt=dp*(1-(i-t1)/t0)**2*math.exp(-omega*(i-t1)/t0)-dp_t
      elif dp>=50:
          dp_t=dp_*14*(i-t1)/t_0*(1-(i-t1)/t_0)*math.exp(-4*(i-t1)/t_0)
          dpt=dp*(0.37*((1+(i-t1)/t0)**(-5/6)+0.63*(1-(i-t1)/t0)**(-7)))-dp_t
      y.append(dpt)
ax2.set_title('Nuclear effect manual')
ax2.set_xlabel('arrival time($s$)')
ax2.set_ylabel('overpressure($psi$)')
plt.plot(t,y)
plt.show()

# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 16:34:32 2020

@author: songyunjia
"""

import numpy as np
#import pylab as plt
import pickle as pk 
import matplotlib.pyplot as plt 
'''
定义初值，
pace是步长
points代表计算的总格点数
lam_0是储存lambda的list
tao_0是储存tau的list
'''
epsilon=1.0
pace=0.0001
points=50000
lam_0=[]
tao_0=[]
# dxk_0=[]
# xk_0=[]


#for tao>=1

def Lam(tau):        #the big Lambda
    return tau**((2/3)+2/(9*epsilon))
##define nu
'''
定义nu函数
因为计算量过大，电脑会算崩，这里调小了nu的精度
'''
def nu(x,lambd,y):
    integrate=0
    for i in range(int(points/1000)):
        if abs(x) >= abs(lambd[i*1000]/Lam(y[i*1000])):
           integrate+=1000*pace/(y[i*1000]**(1+(2/epsilon/3)))
    return integrate*2/(3*epsilon)

'''
lambda的一阶导、二阶导
delta1为lambda的一阶微分
delta2为lambda的二阶微分
'''
def delta1(dlam,tao):  
    return dlam
def delta2(lamb,tao,n): #second derivative of lambda  
    return -(np.pi**2)/8*(tao**(2/3/epsilon))/(lamb**2)*n
'''
猜了nu的初值这里设为x^3
'''
nu_i=[]
def nu_0(x):
    return x


'''
这里计算了nu的initial guess
以及tao_0(在之后的计算中不变)
'''
lam0=1
dlam0=0
for i in range(points):
    t=1+i*pace
    tao_0.append(t)
    lam_0.append(lam0)
    temp_nu=nu_0(lam0/Lam(t))
    nu_i.append(temp_nu)
    dlam0+=delta2(lam0,t,nu_i[i])*pace
    lam0+=delta1(dlam0,t)*pace
'''
判断nu是否收敛
收敛return True
'''
def judge(nu1,nu2):
    standard=0.1
    total=0
    for i in range(int(points)):
        total+=(nu1[i]-nu2[i])**2
    if np.sqrt(total)>standard:
        return False
    return True

fig,ax=plt.subplots() #用于之后作图

#开始迭代
nu1=[]
nu2=[]
nu1=nu_i
j=False   #这个j没有实际意义，只是为了让它循环迭代
k=0    #k只是作为隔几次迭代保存数据的依据
while j==False:
    k+=1
    for i in range(points):
        '''
        算出新的nu
        '''
        temp_nu=nu(lam_0[i]/Lam(tao_0[i]),lam_0,tao_0)
        nu2.append(temp_nu) 
    # nu2.append(temp_nu)
    if judge(nu1,nu2)==True:
        '''
        如果迭代完成就画个图
        '''
        
        plt.xlabel(r'$ \tau $')
        plt.ylabel(r'$ \lambda $')
        plt.plot(tao_0,np.abs(lam_0))
        plt.show()
        break
    else:
        '''
        迭代没有完成就把新nu带进去在迭代一次
        '''
        nu1=nu2
        lam_0=[0]*points
        lam0=1
        dlam0=0
        sign=1
        dlam0_ex=dlam0
        lam0_ex=lam0
        i=0
        while i<=points-1:
            lam_0[i]=lam0
            xk1=delta1(dlam0,tao_0[i])            
            dxk1=delta2(lam0,tao_0[i],nu2[i])
            xk2=delta1(dlam0+dxk1*pace/2.,tao_0[i]+pace/2.)
            dxk2=delta2(lam0+xk1*pace/2.,tao_0[i]+pace/2.,nu2[i])
            xk3=delta1(dlam0+dxk2*pace/2.,tao_0[i]+pace/2.)
            dxk3=delta2(lam0+xk2*pace/2.,tao_0[i]+pace/2.,nu2[i])
            xk4=delta1(dlam0+dxk3*pace,tao_0[i]+pace)
            dxk4=delta2(lam0+xk3*pace,tao_0[i]+pace,nu2[i])
            # dxk_0.append((dxk1+2*dxk2+2*dxk3+dxk4)/6.)
            # xk_0.append((xk1+2*xk2+2*xk3+xk4)/6.)
            dxk=(dxk1+2*dxk2+2*dxk3+dxk4)/6.
            xk=(xk1+2*xk2+2*xk3+xk4)/6.
            dlam0+=sign*pace*dxk
            lam0+=pace*xk
            if abs(dlam0-dlam0_ex)>0.5:
                '''
                判断是否速度跳动过于厉害，
                dlam0_ex是上一个时间的速度。
                逻辑上是，如果跳动过于厉害就退回上一个时间
                使速度等于上一时间速度
                lambda等于上一时间lambda的相反数
                并翻转加速度的方向
                '''
                lam0=-lam0_ex
                dlam0=dlam0_ex
                i-=1
                i+=int(2*abs(lam0)/abs(dlam0))
                sign=-sign
            i+=1
            dlam0_ex=dlam0
            lam0_ex=lam0
            # if lam0<0.1 and lam0>-0.1:
            #     sign=-sign
            #     lam0=-lam0
            # else:
            #     continue
        print(lam_0[1200])   #为了证明这个程序仍然存活
        nu2=[]  #初始化nu
        '''
        下面这段做了个动画
        可以算一轮看一次图像
        如果显示未响应需要等一会儿
        '''
        ax.cla()
        ax.set_xlim(0,11)
        ax.set_ylim(0,3.)
        ax.plot(tao_0,np.abs(lam_0))
        plt.pause(0.05)

        '''
        这段是用来保存数据的
        '''
        if k%10==0:
            with open('A'+str(k),'wb') as f:
               pk.dump(lam_0,f)
            with open('nA'+str(k),'wb') as g:
                pk.dump(nu1,g) 



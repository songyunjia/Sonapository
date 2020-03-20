# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 14:40:49 2020

@author: songyunjia
"""

import numpy as np
import pylab as plt

epsilon=0.8
pace=0.001
points=8000
lam_0=[]
tao_0=[]


#for tao>=1
def gama(tau):        #the big Lambda
    return tau**((2/3)+2/(9*epsilon))
##define nu
def nu(x,lambd,y):
    integrate=0
    for i in range(points):
        if x >= (lambd[i]/gama(y[i])):
           integrate+=pace/(y[i]**(1+(2/epsilon/3)))
    return integrate*2/(3*epsilon)


##
def delta1(dlam,tao):  #d\lambda/d\tau
    return dlam

##initial guess of nu
nu_i=[]
def nu_0(x):
    return x**2
def delta2_0(lamb,tao): #second derivative of lambda under initial guess of nu 
    return -(np.pi**2)/8*(tao**(2/3/epsilon))/(lamb**2)*nu_0(lamb/gama(tao))
lam0=1
dlam0=0
lam=lam0
dlam=dlam0
for i in range(points):
    t=1+i*pace
    tao_0.append(t)
    lam_0.append(lam0)
    nu_i.append(nu_0(lam0/t))
    dlam+=delta2_0(lam0,t)*pace
    lam+=delta1(dlam0,t)*pace
    lam0=lam
    dlam0=dlam

##    
def delta2(lam,tao,n):           #second derivative of lambda
    return -(np.pi**2)/8*(tao**(2/3/epsilon))/(lam**2)*n

##judgement
def judge(nu1,nu2):
    standard=0.01
    for i in range(int(points/5)):
        if abs(nu1[i*5]-nu2[i*5])>standard:
            return False
    return True

##iteration
nu1=[]
nu2=[]
nu1=nu_i
j=False
while j==False:
    for i in range(points):
        nu2.append(nu(lam_0[i]/tao_0[i],lam_0,tao_0)) #to get new nu
    if judge(nu1,nu2)==True:
        plt.plot(tao_0,lam_0)
        plt.xlabel(r'$ \tau $')
        plt.ylabel(r'$ \lambda $')
        break
    else:
        nu1=nu2
        lam_0=[]
        lam0=1
        dlam0=0
        lam=lam0
        dlam=dlam0
        for i in range(points):
            lam_0.append(lam0)
            dxk1=delta2(lam0,tao_0[i],nu2[i])
            dxk2=delta2(lam0+dxk1*pace/2.,tao_0[i]+pace/2.,nu2[i])
            dxk3=delta2(lam0+dxk2*pace/2.,tao_0[i]+pace/2.,nu2[i])
            dxk4=delta2(lam0+dxk3*pace,tao_0[i]+pace,nu2[i])
            dlam+=pace*(dxk1+2*dxk2+2*dxk3+dxk4)/6.
#            dlam+=delta2(lam0,tao_0[i],nu2[i])*pace  
            xk1=delta1(dlam0,tao_0[i])
            xk2=delta1(dlam0+xk1*pace/2.,tao_0[i]+pace/2.)
            xk3=delta1(dlam0+xk2*pace/2.,tao_0[i]+pace/2.)
            xk4=delta1(dlam0+xk3*pace,tao_0[i]+pace)
            lam+=pace*(xk1+2*xk2+2*xk3+xk4)/6.
#            lam+=delta1(dlam0,tao_0[i])*pace
#            if lam<0:
#                lam=0.01
#                dlam=-dlam
            lam0=lam
            dlam0=dlam
        nu2=[]


'''

    
    
#for tao<1
epsilon=0.2
EPS=0.0001
maxi=10000
tao=1
lam0=1
dlam0=0
lam_d=[lam0]
tao_d=[tao]
def nu(t):
    return t**(-2/3/epsilon)
def delta2(lam,tao):
    return -np.pi**2/8*tao**(2/3/epsilon)/lam**2*nu(tao)

def delta1(dlam,tao):
    return dlam
lam=lam0
dlam=dlam0
while tao>0:
    
    dlam-=delta2(lam0,tao)*EPS
    lam-=delta1(dlam0,tao)*EPS
    lam0=lam
    dlam0=dlam
    tao-=EPS
    lam_d.append(lam)
    tao_d.append(tao)

plt.plot(tao_d,lam_d)    
plt.xlim((0,12))
plt.xlabel(r'$ \tau $')
plt.ylabel(r'$ \lambda $')
plt.ylim((0,1.1))
 '''   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
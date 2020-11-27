# -*- coding: utf-8 -*-
"""Untitled7.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/13bcbGmkwflZYw2_XottDuGMJVHLhoZmx
"""

#Author(s) : Hrithik Agarwal, Hrishikesh Kotwal
#This code solve Q3 of Group Assignment 1 of ME608 course

#import libraries 
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy import interpolate
 
hf = 4*pow(10,7)          #The enthalpy of formation of fuel which is ethane in our case
cp = 1200                 #specific heat of fuel, air and products (constant and equal)
Ru = 8315                 #Universal gas constant in J/Kmol.K
Ru1 = 287                 #Universal gas constant in J/Kg.K

P0 = 25.12*101325         #initial pressure
T0 = 753                  #initial temperature
F0 = (1/17)*P0/(Ru*T0)    #molar concentration of fuel at t=0 ms
Ox0 = (16/17)*P0/(Ru*T0)  #molar concentration of Oxidizer at t=0 ms
Pr0 = 0                   #molar concentration of Product at t=0 ms

t1 = np.linspace(0,0.003,10)   # Creating vector t, time points

#function that returns dF/dt,dOx/dt,dPr/dt,dT/dt,dP/d
#Function that returns derivative values at requested F,Ox,Pr,T,P and t values
def module1(z,t,hf,cp,Ru,Ru1):
    F = z[0]
    Ox= z[1]
    Pr = z[2]
    T = z[3]
    P = z[4]

    dF_dt = -1*(6.19*pow(10,9))*np.exp(-1*15098/T)*pow(F,0.1)*pow(0.21*Ox,1.65) #expression for dF/dt
    dOx_dt = 16*dF_dt                                                           #expression for dOx/dt
    dPr_dt = -17*dF_dt                                                          #expression for dPr/dt
    dT_dt = (-1*dF_dt*hf*Ru*T)/((cp-Ru1)*P)                                     #expression for dT/dt
    dP_dt = (P0/T0)*dT_dt                                                       #expression for dP/dt
    

    return [dF_dt,dOx_dt,dPr_dt,dT_dt,dP_dt]


# initial condition, condition at t=0 ms
z0 = [F0,Ox0,Pr0,T0,P0]

# solve ODE
z = odeint(module1, z0, t1, args=(hf,cp,Ru,Ru1))

F1 = z[:,0]                #vector that stores concentration of fuel at different time period between 0-3 ms     
Ox1 = z[:,1]               #vector that stores concentration of Oxidizer at different time period between 0-3 ms
Pr1 = z[:,2]               #vector that stores concentration of Product at different time period between 0-3 ms
T1 = z[:,3]                #vector that stores concentration of Temperature at different time period between 0-3 ms
P1 = z[:,4]                #vector that stores concentration of Pressure at different time period between 0-3 ms


F0 = z[:,0][-1]           #value of concentration of fuel at 3 ms
Ox0 = z[:,1][-1]          #value of concentration of Oxidizer at 3 ms
Pr0 = z[:,2][-1]          #value of concentration of product at 3 ms
T0= z[:,3][-1]            #value of concentration of temperature at 3 ms
P0 = z[:,4][-1]           #value of concentration of pressure at 3 ms

t2 = np.linspace(0.003,0.00309,5)    # Creating vector t, time points

#function that returns dF/dt,dOx/dt,dPr/dt,dT/dt,dP/d
#Function that returns derivative values at requested F,Ox,Pr,T,P and t values
def module2(z,t,hf,cp,Ru,Ru1):
    F = z[0]
    Ox= z[1]
    Pr = z[2]
    T = z[3]
    P = z[4]

    dF_dt = -1*(6.19*pow(10,9))*np.exp(-1*15098/T)*pow(F,0.1)*pow(0.21*Ox,1.65) #expression for dF/dt
    dOx_dt = 16*dF_dt                                                           #expression for dOx/dt
    dPr_dt = -17*dF_dt                                                          #expression for dPr/dt
    dT_dt = (-1*dF_dt*hf*Ru*T)/((cp-Ru1)*P)                                     #expression for dT/dt
    dP_dt = (P0/T0)*dT_dt                                                       #expression for dP/dt
    
    return [dF_dt,dOx_dt,dPr_dt,dT_dt,dP_dt]

# initial condition, condition at t=3 ms
z0 = [F0,Ox0,Pr0,T0,P0]

# solve ODE
z = odeint(module2, z0, t2, args=(hf,cp,Ru,Ru1))


F2 = z[:,0]                 #vector that stores concentration of fuel at different time period between 3-3.09 ms
Ox2 = z[:,1]                #vector that stores concentration of Oxidizer at different time period between 3-3.09 ms
Pr2 = z[:,2]                #vector that stores concentration of Product at different time period between 3-3.09 ms
T2 = z[:,3]                 #vector that stores concentration of Temp at different time period between 3-3.09 ms
P2 = z[:,4]                 #vector that stores concentration of pressure at different time period between 3-3.09 ms

F_pre = np.concatenate((F1,F2))                   #vector that stores concentration of fuel at different time period between 0-3.09 ms
Ox_pre = np.concatenate((Ox1,Ox2))                #vector that stores concentration of Oxidizer at different time period between 0-3.09 ms
Pr_pre = np.concatenate((Pr1,Pr2))                #vector that stores concentration of Product at different time period between 0-3.09 ms
T_pre = np.concatenate((T1,T2))                   #vector that stores concentration of Temp at different time period between 0-3.09 ms
P_pre = np.concatenate((P1,P2))                   #vector that stores concentration of pressure at different time period between 0-3.09 ms
t_pre = np.concatenate((t1,t2))                   #time points between 0-3.09ms

F0 = z[:,0][-1]           #value of concentration of fuel at 3.09 ms
Ox0 = z[:,1][-1]          #value of concentration of Oxidizer at 3.09 ms
Pr0 = z[:,2][-1]          #value of concentration of product at 3.09 ms
T0= z[:,3][-1]            #value of concentration of temperature at 3.09 ms
P0 = z[:,4][-1]           #value of concentration of pressure at 3.09 ms

t3 = np.linspace(0.00309,0.0031,5)   # Creating vector t, time points

#function that returns dF/dt,dOx/dt,dPr/dt,dT/dt,dP/d
#Function that returns derivative values at requested F,Ox,Pr,T,P and t values
def module3(z,t,hf,cp,Ru,Ru1):
    F = z[0]
    Ox= z[1]
    Pr = z[2]
    T = z[3]
    P = z[4]

    dF_dt = -1*(6.19*pow(10,9))*np.exp(-1*15098/T)*pow(F,0.1)*pow(0.21*Ox,1.65) #expression for dF/dt
    dOx_dt = 16*dF_dt                                                           #expression for dOx/dt
    dPr_dt = -17*dF_dt                                                          #expression for dPr/dt
    dT_dt = (-1*dF_dt*hf*Ru*T)/((cp-Ru1)*P)                                     #expression for dT/dt
    dP_dt = (P0/T0)*dT_dt                                                       #expression for dP/dt

    return [dF_dt,dOx_dt,dPr_dt,dT_dt,dP_dt]

# initial condition, condition at t=3.09 ms
z0 = [F0,Ox0,Pr0,T0,P0]

# solve ODE
z = odeint(module3, z0, t3, args=(hf,cp,Ru,Ru1))

F3 = z[:,0]                 #vector that stores concentration of fuel at different time period between 3.09-3.1 ms
Ox3 = z[:,1]                #vector that stores concentration of Oxidizer at different time period between 3.09-3.1 ms
Pr3 = z[:,2]                #vector that stores concentration of Product at different time period between 3.09-3.1 ms
T3 = z[:,3]                 #vector that stores concentration of Temp at different time period between 3.09-3.1 ms
P3 = z[:,4]                 #vector that stores concentration of pressure at different time period between 3.09-3.1 ms

F = np.concatenate((F_pre,F3))                   #vector that stores concentration of fuel at different time period between 0-3.1 ms
Ox = np.concatenate((Ox_pre,Ox3))                #vector that stores concentration of Oxidizer at different time period between 0-3.1 ms
Pr = np.concatenate((Pr_pre,Pr3))                #vector that stores concentration of Product at different time period between 0-3.1 ms
T = np.concatenate((T_pre,T3))                   #vector that stores concentration of Temp at different time period between 0-3.1 ms
P = np.concatenate((P_pre,P3))                   #vector that stores concentration of pressure at different time period between 0-3.1 ms
t = np.concatenate((t_pre,t3))                   #time points between 0-3.1 ms

#calculating dP/dt

dpdt1=[]             #vector that stores dP/dt at different time period between 0-3 ms
dpdt2 = []           #vector that stores dP/dt at different time period between 3-3.09 ms
dpdt3 = []           #vector that stores dP/dt at different time period between 3.09-3.1 ms
dpdt = []            #vector that stores dP/dt at different time period between 0-3.1 ms
for i in range(0,10):
    if round ((t[i+1]-t[i]),10)==0:
      dpdt1.append(1)
    else:
      dpdt1.append((P[i+1]-P[i])/(t[i+1]-t[i]))

for i in range(10,15):
    if round ((t[i+1]-t[i]),10)==0:
      dpdt2.append(1)
    else:
      dpdt2.append((P[i+1]-P[i])/(t[i+1]-t[i]))

for i in range(15,19):
    if round ((t[i+1]-t[i]),10)==0:
      dpdt3.append(1)
    else:
      dpdt3.append((P[i+1]-P[i])/(t[i+1]-t[i]))
dpdt3.append ((P[-1]-P[-2])/(t[-1]-t[-2]))

for i in range(0,19):
    if round ((t[i+1]-t[i]),10)==0:
      dpdt.append(1)
    else:
      dpdt.append((P[i+1]-P[i])/(t[i+1]-t[i]))
dpdt.append ((P[-1]-P[-2])/(t[-1]-t[-2]))

#if you want to see values of F,Ox,Pr,T,P,t vector, uncomment that line accordingly
#print(F)    #printing vector F 
#print(Ox)    #printing vector Ox
#print(Pr)    #printing vector Pr
#print(T)    #printing vector T
#print(P)    #printing vector P
#print(t)    #printing vector t


# plot results
plt.plot(t1,F1,color = '#FF0800',label = 'Fuel concentration')
plt.xlabel('Time in second')
plt.legend()
plt.title('Fuel Concentration vs Time')
plt.show()

plt.plot(t2,F2,color = '#FF0800',label = 'Fuel concentration')
plt.xlabel('Time in second')
plt.legend()
plt.title('Fuel Concentration vs Time (scale change)')
plt.show()

plt.plot(t3,F3,color = '#FF0800',label = 'Fuel concentration')
plt.xlabel('Time in second')
plt.legend()
plt.title('Fuel Concentration vs Time (scale change)')
plt.show()

plt.plot(t,F,color = '#FF0800',label = 'Fuel concentration')
plt.xlabel('Time in second')
plt.legend()
plt.title('Complete graph of Fuel Concentration vs Time')
plt.show()

plt.plot(t1,Pr1,color = '#FF0800',label = 'Product Concentration')
plt.xlabel('Time in second')
plt.title('Product Concentration vs Time')
plt.legend()
plt.show()

plt.plot(t2,Pr2,color = '#FF0800',label = 'Product Concentration')
plt.xlabel('Time in second')
plt.title('Product Concentration vs Time(scale change)')
plt.legend()
plt.show()

plt.plot(t3,Pr3,color = '#FF0800',label = 'Product Concentration')
plt.xlabel('Time in second')
plt.title('Product Concentration vs Time(scale change)')
plt.legend()
plt.show()

plt.plot(t,Pr,color = '#FF0800',label = 'Product Concentration')
plt.xlabel('Time in second')
plt.title('Complete graph of Product Concentration vs Time')
plt.legend()
plt.show()

plt.plot(t1,dpdt1,color = '#FF0800',label = 'Product Concentration')
plt.xlabel('Time in second')
plt.legend()
plt.title('dP/dt vs Time')
plt.show()

plt.plot(t2,dpdt2,color = '#FF0800',label = 'Product Concentration')
plt.xlabel('Time in second')
plt.legend()
plt.title('dP/dt Time (scale change)')
plt.show()

plt.plot(t3,dpdt3,color = '#FF0800',label = 'Product Concentration')
plt.xlabel('Time in second')
plt.legend()
plt.title('dP/dt vs Time(scale change)')
plt.show()

plt.plot(t,dpdt,color = '#FF0800',label = 'Product Concentration')
plt.xlabel('Time in second')
plt.legend()
plt.title('Complete graph of dP/dt vs Time')
plt.show()

#plt.plot(t,P,color ='#FF0800',label = 'Pressure')
#plt.xlabel('Time in second')
#plt.legend()
#plt.title('Pressure vs Time')
#plt.show()

plt.plot(t1,T1,color = '#FF0800',label = 'Temprature')
plt.xlabel('Time in second')
plt.legend()
plt.title('Temprature vs Time')
plt.show()

plt.plot(t2,T2,color = '#FF0800',label = 'Temprature')
plt.xlabel('Time in second')
plt.legend()
plt.title('Temprature vs Time(scale change)')
plt.show()

plt.plot(t3,T3,color = '#FF0800',label = 'Temprature')
plt.xlabel('Time in second')
plt.legend()
plt.title('Temprature vs Time(scale change)')
plt.show()

plt.plot(t,T,color = '#FF0800',label = 'Temprature')
plt.xlabel('Time in second')
plt.legend()
plt.title('Complete graph of Temprature vs Time')
plt.show()
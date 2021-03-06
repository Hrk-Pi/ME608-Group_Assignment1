# -*- coding: utf-8 -*-
"""Q2.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1aARF2kHHMxxGa2ii15saiuF0m7BYSdXr
"""

#Author(s) : Hrithik Agarwal, Hrishikesh Kotwal
#This code solve Q2 of Group Assignment 1 of ME608 course

#import libraries
import numpy as np
import sympy  as sp
import scipy as sc
from sympy import *
import math as mp
import matplotlib.pyplot as plt


m = 3                           #It is in form of CmHn, where C is carbon and H is hydrogen
n = 8                          
phi = 1                         #equivalence ratio
del_H = 4.65*(10**7)            # heat of combustion
cp = 1230                       # specific Heat at constant pressure
astoich = m+n/4                 # stichometric A/F
mmix = 29                       # mass of mixture
Earv=15098                      # activation energy/universal gas constant
vol = (3.14*((0.08)**3))/6      # volume of well     
Af =4.836*(10**9)               # constant from mass flow rate of fuel   
X=0.1;                          # constant from mass flow rate of fuel
Y=1.65;                         # constant from mass flow rate of fuel
t1 = 298                        # inlet temprature
p1 = 101325                     # inlet pressure
mf= 44                          # mass o fuel
mo = 32                         # mass of oxygen
mn= 28                          # mass of nitrogen
m_u = mo*astoich/mf             
ru = 8315                       # universal gas constant
m_r = 0.1                       # mass flow rate


# initial conditions 

yfi = mf/(mf + (astoich/phi)*(mo + 3.76*mn))                         # mass fraction fuel                            
yoi=(astoich/phi)*(mo/(mf + (astoich/phi)*(mo + 3.76*mn)))           # mass fraction oxygen
yni= 1-yfi-yoi                                                       # mass fraction nitrogen

# graph ploting

T1 = np.linspace(0,3000,500)
xaxis = [0]*(500)
f= -T1 + t1 + (del_H*vol/(cp*m_r)) * (Af * np.exp(-Earv/T1) * ((p1*mmix/(ru*T1))**(X+Y))* (mf**(1-X)) * (mo**(-Y)) * ((yfi -(cp/del_H)*(T1-t1))**(X)) * ((yoi -(m_u*cp/del_H)*(T1-t1))**(Y)) )
plt.title('Plot of F(T) (mass flow rate=0.1)' )
plt.plot(T1,f)
plt.plot(T1,xaxis)
plt.show()

print("We can see from graph that F(T) becomes 0 for 3 Temperatures(T1) and using Newton Raphson Method and different initial guess, we will try to find out those temperature", "\n")
print("Tf we get a complex anwer then F(T) dosent cross x axis givin complex solution and reason might be blow off condtioin is reached or reaction dosent take place ","\n")

# sumbolic initialisation
T = symbols('T')
f1= -T + t1 + (del_H*vol/(cp*m_r)) * (Af * exp(-Earv/T) * ((p1*mmix/(ru*T))**(X+Y))* (mf**(1-X)) * (mo**(-Y)) * ((yfi -(cp/del_H)*(T-t1))**(X)) * ((yoi -(m_u*cp/del_H)*(T-t1))**(Y)) )
c1=Derivative(f1, T).doit()         # finda derivative

initial_guess = 2400                  # initial guess, if you want to change initial guess for Temperature then change here

error=100                            # initail place taker
tol = 10**(-4)                        #tolerance

ig=initial_guess


while error > tol: 
    h=float(f1.subs({T:initial_guess}))                   # value of F(T)
    k=float(c1.subs({T:initial_guess}))                   # value of d(F(T))/d(T) 
    titr =initial_guess - h/k                             # titer = T at iteration   .... newton raphson formula 
    error = abs(titr - initial_guess)                     # error calculation
    initial_guess = titr                                  # changing initial guess to preious iteration answer   

print("Temperature where solution will converge with given initial guess (" + str(ig) +") is :" + str(titr))

############################## mass flow rate = 0.2

m_r = 0.2

T1 = np.linspace(0,3000,500)
xaxis = [0]*(500)
f= -T1 + t1 + (del_H*vol/(cp*m_r)) * (Af * np.exp(-Earv/T1) * ((p1*mmix/(ru*T1))**(X+Y))* (mf**(1-X)) * (mo**(-Y)) * ((yfi -(cp/del_H)*(T1-t1))**(X)) * ((yoi -(m_u*cp/del_H)*(T1-t1))**(Y)) )
plt.title('Plot of F(T) (mass flow rate=0.2)' )
plt.plot(T1,f)
plt.plot(T1,xaxis)
plt.show()

f1= -T + t1 + (del_H*vol/(cp*m_r)) * (Af * exp(-Earv/T) * ((p1*mmix/(ru*T))**(X+Y))* (mf**(1-X)) * (mo**(-Y)) * ((yfi -(cp/del_H)*(T-t1))**(X)) * ((yoi -(m_u*cp/del_H)*(T-t1))**(Y)) )
c1=Derivative(f1, T).doit() 

initial_guess = 2400                  # initial guess, if you want to change initial guess for Temperature then change here

error=100                            # initail place taker
tol = 10**(-4)                        #tolerance

ig=initial_guess


while error > tol: 
    h=float(f1.subs({T:initial_guess}))                   # value of F(T)
    k=float(c1.subs({T:initial_guess}))                   # value of d(F(T))/d(T) 
    titr =initial_guess - h/k                             # titer = T at iteration   .... newton raphson formula 
    error = abs(titr - initial_guess)                     # error calculation
    initial_guess = titr                                  # changing initial guess to preious iteration answer   

print("Temperature where solution will converge with given initial guess (" + str(ig) +") is :" + str(titr))


############################## mass flow rate = 0.3

m_r = 0.3

T1 = np.linspace(0,3000,500)
xaxis = [0]*(500)
f= -T1 + t1 + (del_H*vol/(cp*m_r)) * (Af * np.exp(-Earv/T1) * ((p1*mmix/(ru*T1))**(X+Y))* (mf**(1-X)) * (mo**(-Y)) * ((yfi -(cp/del_H)*(T1-t1))**(X)) * ((yoi -(m_u*cp/del_H)*(T1-t1))**(Y)) )
plt.title('Plot of F(T) (mass flow rate=0.3)' )
plt.plot(T1,f)
plt.plot(T1,xaxis)
plt.show()

f1= -T + t1 + (del_H*vol/(cp*m_r)) * (Af * exp(-Earv/T) * ((p1*mmix/(ru*T))**(X+Y))* (mf**(1-X)) * (mo**(-Y)) * ((yfi -(cp/del_H)*(T-t1))**(X)) * ((yoi -(m_u*cp/del_H)*(T-t1))**(Y)) )
c1=Derivative(f1, T).doit() 

initial_guess = 2400                  # initial guess, if you want to change initial guess for Temperature then change here

error=100                            # initail place taker
tol = 10**(-4)                        #tolerance

ig=initial_guess


while error > tol: 
    h=float(f1.subs({T:initial_guess}))                   # value of F(T)
    k=float(c1.subs({T:initial_guess}))                   # value of d(F(T))/d(T) 
    titr =initial_guess - h/k                             # titer = T at iteration   .... newton raphson formula 
    error = abs(titr - initial_guess)                     # error calculation
    initial_guess = titr                                  # changing initial guess to preious iteration answer   

print("Temperature where solution will converge with given initial guess (" + str(ig) +") is :" + str(titr))


############################## mass flow rate = 0.4

m_r = 0.4

T1 = np.linspace(0,3000,500)
xaxis = [0]*(500)
f= -T1 + t1 + (del_H*vol/(cp*m_r)) * (Af * np.exp(-Earv/T1) * ((p1*mmix/(ru*T1))**(X+Y))* (mf**(1-X)) * (mo**(-Y)) * ((yfi -(cp/del_H)*(T1-t1))**(X)) * ((yoi -(m_u*cp/del_H)*(T1-t1))**(Y)) )
plt.title('Plot of F(T) (mass flow rate=0.4)' )
plt.plot(T1,f)
plt.plot(T1,xaxis)
plt.show()

f1= -T + t1 + (del_H*vol/(cp*m_r)) * (Af * exp(-Earv/T) * ((p1*mmix/(ru*T))**(X+Y))* (mf**(1-X)) * (mo**(-Y)) * ((yfi -(cp/del_H)*(T-t1))**(X)) * ((yoi -(m_u*cp/del_H)*(T-t1))**(Y)) )
c1=Derivative(f1, T).doit() 

initial_guess = 2400                  # initial guess, if you want to change initial guess for Temperature then change here

error=100                            # initail place taker
tol = 10**(-4)                        #tolerance

ig=initial_guess


while error > tol: 
    h=float(f1.subs({T:initial_guess}))                   # value of F(T)
    k=float(c1.subs({T:initial_guess}))                   # value of d(F(T))/d(T) 
    titr =initial_guess - h/k                             # titer = T at iteration   .... newton raphson formula 
    error = abs(titr - initial_guess)                     # error calculation
    initial_guess = titr                                  # changing initial guess to preious iteration answer   

print("Temperature where solution will converge with given initial guess (" + str(ig) +") is :" + str(titr))

################################# mass flow rate = 0.41

m_r = 0.41

T1 = np.linspace(0,3000,500)
xaxis = [0]*(500)
f= -T1 + t1 + (del_H*vol/(cp*m_r)) * (Af * np.exp(-Earv/T1) * ((p1*mmix/(ru*T1))**(X+Y))* (mf**(1-X)) * (mo**(-Y)) * ((yfi -(cp/del_H)*(T1-t1))**(X)) * ((yoi -(m_u*cp/del_H)*(T1-t1))**(Y)) )
plt.title('Plot of F(T) (mass flow rate=0.41)' )
plt.plot(T1,f)
plt.plot(T1,xaxis)
plt.show()

f1= -T + t1 + (del_H*vol/(cp*m_r)) * (Af * exp(-Earv/T) * ((p1*mmix/(ru*T))**(X+Y))* (mf**(1-X)) * (mo**(-Y)) * ((yfi -(cp/del_H)*(T-t1))**(X)) * ((yoi -(m_u*cp/del_H)*(T-t1))**(Y)) )
c1=Derivative(f1, T).doit() 

initial_guess = 2400                  # initial guess, if you want to change initial guess for Temperature then change here

error=100                            # initail place taker
tol = 10**(-4)                        #tolerance

ig=initial_guess


while error > tol: 
    h=float(f1.subs({T:initial_guess}))                   # value of F(T)
    k=float(c1.subs({T:initial_guess}))                   # value of d(F(T))/d(T) 
    titr =initial_guess - h/k                             # titer = T at iteration   .... newton raphson formula 
    error = abs(titr - initial_guess)                     # error calculation
    initial_guess = titr                                  # changing initial guess to preious iteration answer   

print("Temperature where solution will converge with given initial guess (" + str(ig) +") is :" + str(titr))


############################## mass flow rate = 0.42

m_r = 0.42

T1 = np.linspace(0,3000,500)
xaxis = [0]*(500)
f= -T1 + t1 + (del_H*vol/(cp*m_r)) * (Af * np.exp(-Earv/T1) * ((p1*mmix/(ru*T1))**(X+Y))* (mf**(1-X)) * (mo**(-Y)) * ((yfi -(cp/del_H)*(T1-t1))**(X)) * ((yoi -(m_u*cp/del_H)*(T1-t1))**(Y)) )
plt.title('Plot of F(T) (mass flow rate=0.42)' )
plt.plot(T1,f)
plt.plot(T1,xaxis)
plt.show()

#f1= -T + t1 + (del_H*vol/(cp*m_r)) * (Af * exp(-Earv/T) * ((p1*mmix/(ru*T))**(X+Y))* (mf**(1-X)) * (mo**(-Y)) * ((yfi -(cp/del_H)*(T-t1))**(X)) * ((yoi -(m_u*cp/del_H)*(T-t1))**(Y)) )
#c1=Derivative(f1, T).doit() 

#initial_guess = 2400                  # initial guess, if you want to change initial guess for Temperature then change here

#error=100                            # initail place taker
#tol = 10**(-4)                        #tolerance

#ig=initial_guess


#while error > tol: 
    #h=float(f1.subs({T:initial_guess}))                   # value of F(T)
    #k=float(c1.subs({T:initial_guess}))                   # value of d(F(T))/d(T) 
    #itr =initial_guess - h/k                             # titer = T at iteration   .... newton raphson formula 
    #error = abs(titr - initial_guess)                     # error calculation
    #initial_guess = titr                                  # changing initial guess to preious iteration answer   

#print("Temperature where solution will converge with given initial guess (" + str(ig) +") is :" + str(titr))


print("As Solution is complex we got our blow off limit")


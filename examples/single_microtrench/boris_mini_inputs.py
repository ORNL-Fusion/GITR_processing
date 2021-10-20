#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 15:34:09 2021

@author: Alyssa
"""

from boris_mini import *
import numpy as np
from numpy.linalg import norm, inv
from matplotlib import pyplot as plt
#%% example inputs

# # ----------example1 inputs----------
# tmax0 = 10 #milliseconds
# tsteps0 = 1000
# m0 = 184
# q0 = 4

# # initial velocity
# v0 = np.zeros(3)
# v0[1] = 100 #m/s

# # E-field unit vector
# ee = np.zeros(3)
# ee[1] = 1 + ee[1]
# Emag = 1e3 #V/m
# E0 = Emag*ee

# # B-field unit vector
# b = np.ones(3)
# b[0] = (sqrt(3)/4) * b[0]
# b[1] = -0.5 * b[1]
# b[2] = -0.75 * b[2]
# Bmag = 2 #Tesla
# B0 = Bmag*b

# -------example2 inputs-------
mp = 1.6605E-27
me = 9.109383701528e-31
qe = 1.602176634e-19

m = 184*mp
q = 4*qe

# define B field
B0 = 1 #T
alpha = 5 * np.pi/180
#B = B0*np.array([0,0,-1])
B = B0*np.array([np.cos(alpha),0,-np.sin(alpha)])

# rotate B=[Bx,By,Bz] into B'=[0,0,Bz']
theta = np.pi/2-alpha
M1 = rotate_matrix(np.array([0,1,0]), theta)
Bprime = np.dot(M1,B).astype('float16')
print('B before =', B)
print('Bprime = ', Bprime)

# define timing
omega = np.linalg.norm(B)*q/m #plasma frequency
t=np.arange(0,0.00001,1/(1e2*omega))
vmag=1e4
phi=0

# Larmor radius [m] = sqrt(T_D[J]/m_D[kg]) / cyclotron frequency [1/s]
Larmor = np.sqrt(1.6022e-18/3.3211e-27)/omega

# define E field
# z_axis = np.arange(0,0.001,1e-6)
# Ex = np.zeros(len(z_axis))
# Ey = np.zeros(len(z_axis))
# Ez = -3 * (0.8/Larmor) * np.exp(-(0.8/Larmor)*z_axis)
# E = np.dstack([Ex,Ey,Ez])

# E_prime = np.zeros(E[0].shape)
# for i in range(0,len(E[0])): E_prime[i] = np.dot(M1,E[0][i])

E = np.array([0,1e3,0]) #V/m
Eprime = np.dot(M1,E).astype('float16')
print('E before =', E)
print('Eprime =', Eprime)

vE_prime = np.cross(Eprime,Bprime)/(np.linalg.norm(Bprime)**2)

# initial velocity
v0 = np.array([vmag,0,0]) #m/s
v0prime = np.dot(M1,v0)

# calculate primed velocity
vsign = Bprime[-1]/norm(Bprime)
vinit = v0prime - vE_prime
vx_prime = norm(vinit) * np.cos(omega*t+phi) + vE_prime[0]
vy_prime = -vsign*norm(vinit) * np.sin(omega*t+phi) + vE_prime[1]
v_prime = np.dstack([vx_prime,vy_prime,np.zeros(len(vx_prime))])

# calculate primed x,y,z
xxx_prime = norm(vinit)/omega * np.sin(omega*t+phi) + vE_prime[0]*t
yyy_prime = vsign*norm(vinit)/omega * (np.cos(omega*t+phi)-1) + vE_prime[1]*t
xyz_prime = np.dstack([xxx_prime,yyy_prime,np.zeros(len(xxx_prime))])

# rotate velocity into lab space
v_soln = np.zeros(v_prime[0].shape)
for i in range(0,len(v_prime[0])): v_soln[i] = np.dot(inv(M1),v_prime[0][i])
vx = v_soln[:,0]
vy = v_soln[:,1]
vz = v_soln[:,2]

# rotate x,y,z into lab space
xyz_soln = np.zeros(xyz_prime[0].shape)
for i in range(0,len(xyz_prime[0])): xyz_soln[i] = np.dot(inv(M1),xyz_prime[0][i])
xxx = xyz_soln[:,0]
yyy = xyz_soln[:,1]
zzz = xyz_soln[:,2]

# Add electrostatic force on particle by projecting E-field onto B-field and adding the parallel Fe
Eparallel = np.dot(E,B)
v_electrostatic = q*Eparallel*t/m



v0 = np.array([vx[0],vy[0],vz[0]])
vout1, xyz1 = ParticlePusher(t, m,q, v0,E,B)
vout2, xyz2 = SimplePusher(t, m,q, v0,E,B)

#%% plot vx(t) and vy(t)

# vx(t)
fig,ax = plt.subplots(1)
ax.plot(t,vx, t,vout1[:,0], t,vout2[:,0])
ax.set_xlim(([0,8*np.pi/omega]))
ax.set_title("vx(t)")
plt.legend(['analytic','boris','simple push'])

# vy(t)
fig,ax = plt.subplots(1)
ax.plot(t,vy, t,vout1[:,1], t,vout2[:,1])
ax.set_xlim(([0,8*np.pi/omega]))
ax.set_title("vy(t)")
plt.legend(['analytic','boris','simple push'])

# vz(t)
fig,ax = plt.subplots(1)
ax.plot(t,vz, t,vout1[:,2], t,vout2[:,2])
ax.set_xlim(([0,8*np.pi/omega]))
ax.set_title("vz(t)")
plt.legend(['analytic','boris','simple push'])

#%% plot x(t) and y(t)

# x(t)
fig,ax = plt.subplots(1)
ax.plot(t,xxx, t,xyz1[:,0], t,xyz2[:,0])
ax.set_xlim(([0,8*np.pi/omega]))
ax.set_title("x(t)")
plt.legend(['analytic','boris','simple push'])

# y(t)
fig,ax = plt.subplots(1)
ax.plot(t,yyy, t,xyz1[:,1], t,xyz2[:,1])
ax.set_xlim(([0,8*np.pi/omega]))
ax.set_title("y(t)")
plt.legend(['analytic','boris','simple push'])

#%% plot z(t)

# z(t)
fig,ax = plt.subplots(1)
ax.plot(t,zzz, t,xyz1[:,2], t,xyz2[:,2])
ax.set_xlim(([0,8*np.pi/omega]))
ax.set_title("z(t)")
plt.legend(['analytic','boris','simple push'])

#%% find error

boris_error = np.array([rmse(vout1[:,0],vx), rmse(vout1[:,1],vy)])
simple_error = np.array([rmse(vout2[:,0],vx), rmse(vout2[:,1],vy)])

print("boris rmse = ", boris_error)
print("simple rmse = ", simple_error)

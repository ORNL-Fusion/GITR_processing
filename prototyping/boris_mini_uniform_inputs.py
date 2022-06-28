#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  8 15:34:09 2021

@author: Alyssa
"""

from boris_mini_uniform import *
import numpy as np
from numpy.linalg import norm, inv
from matplotlib import pyplot as plt

#%% example inputs

mp = 1.6605E-27
me = 9.109383701528e-31
qe = 1.602176634e-19

# assume tungsten +4
m = 184*mp
q = 4*qe

# define B field
Bmag = 1 #T
#B = Bmag*np.array([0,0,-1])
alpha = 20 * np.pi/180 #angle of incidence
B = Bmag*np.array([np.cos(alpha),0,-np.sin(alpha)])

# rotate B=[Bx,By,Bz] into B'=[0,0,Bz']
theta = np.pi/2-alpha
M1 = rotate_matrix(np.array([0,1,0]), theta)
Bprime = np.dot(M1,B).astype('float16')
print('B before =', B)
print('Bprime = ', Bprime)

# define timing
omega = Bmag*q/m #plasma frequency
t=np.arange(0,0.0001,1/(1e2*omega))

# define E field
E = np.array([0, 0, 0])
#E = np.array([5e2,1e3,1e3]) #V/m
Eprime = np.dot(M1,E).astype('float16')
print('E before =', E)
print('Eprime =', Eprime)

# define initial velocity
vmag=1e4
v0 = np.array([vmag,0,0]) #m/s
v0prime = np.dot(M1,v0)
vE_prime = np.cross(Eprime,Bprime)/(Bmag**2)

# calculate primed velocity
vsign = Bprime[-1]/Bmag
vperp = v0prime - vE_prime
vx_prime = -vsign*norm(vperp) * np.cos(omega*t) + vE_prime[0]
vy_prime = -vsign*norm(vperp) * np.sin(omega*t) + vE_prime[1]
v_prime = np.dstack([vx_prime,vy_prime,np.zeros(len(vx_prime))])

# calculate primed x,y,z
xxx_prime = norm(vperp)/omega * np.sin(omega*t) + vE_prime[0]*t
yyy_prime = vsign*norm(vperp)/omega * (np.cos(omega*t)-1) + vE_prime[1]*t
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

# add electrostatic force on particle by projecting E-field onto B-field and adding the parallel F_E
Epara = np.dot(E,B)
vpara = q*Epara/Bmag * t/m
vx += vpara * B[0]/Bmag
vy += vpara * B[1]/Bmag
vz += vpara * B[2]/Bmag

xxx += vpara * B[0]/Bmag * t/2
yyy += vpara * B[1]/Bmag * t/2
zzz += vpara * B[2]/Bmag * t/2

print("\n")
print("velocity and position initial conditions")
print("v0 = [", vx[0], vy[0], vz[0], "]")
print("pos0 = [", xxx[0], yyy[0], zzz[0], "]")

# run boris and simple pushers using same initial velocity as analytic example
v0 = np.array([vx[0],vy[0],vz[0]])
vout1, xyz1 = ParticlePusher(t, m,q, v0,E,B)
vout2, xyz2 = SimplePusher(t, m,q, v0,E,B)

#%% plot vx(t),  vy(t), and vz(t)

# vx(t)
fig,ax = plt.subplots(1)
ax.plot(t,vx, t,vout1[:,0], t,vout2[:,0])
ax.set_title("vx(t)")
plt.legend(['analytic','boris','simple push'])

# vy(t)
fig,ax = plt.subplots(1)
ax.plot(t,vy, t,vout1[:,1], t,vout2[:,1])
ax.set_title("vy(t)")
plt.legend(['analytic','boris','simple push'])

# vz(t)
fig,ax = plt.subplots(1)
ax.plot(t,vz, t,vout1[:,2], t,vout2[:,2])
ax.set_title("vz(t)")
plt.legend(['analytic','boris','simple push'])

#%% plot x(t), y(t), and z(t)

# x(t)
fig,ax = plt.subplots(1)
ax.plot(t,xxx, t,xyz1[:,0], t,xyz2[:,0])
ax.set_title("x(t)")
plt.legend(['analytic','boris','simple push'])

# y(t)
fig,ax = plt.subplots(1)
ax.plot(t,yyy, t,xyz1[:,1], t,xyz2[:,1])
ax.set_title("y(t)")
plt.legend(['analytic','boris','simple push'])

# z(t)
fig,ax = plt.subplots(1)
ax.plot(t,zzz, t,xyz1[:,2], t,xyz2[:,2])
ax.set_title("z(t)")
plt.legend(['analytic','boris','simple push'])

#%% 3D plot
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot(xxx, yyy, zzz)
ax.plot(xyz1[:,0], xyz1[:,1], xyz1[:,2])
ax.set_title("position")
plt.legend(['analytic','boris'])

#%% find error

boris_v_error = np.array([rmse(vout1[:,0],vx), rmse(vout1[:,1],vy), rmse(vout1[:,2],vz)])
boris_xyz_error = np.array([rmse(xyz1[:,0],xxx), rmse(xyz1[:,1],yyy), rmse(xyz1[:,2],zzz)])

simple_v_error = np.array([rmse(vout2[:,0],vx), rmse(vout2[:,1],vy), rmse(vout2[:,2],vz)])
simple_xyz_error = np.array([rmse(xyz2[:,0],xxx), rmse(xyz2[:,1],yyy), rmse(xyz2[:,2],zzz)])

print("\n")
print("rmse errors")
print("boris v rmse =", boris_v_error)
print("boris xyz rmse =", boris_xyz_error)
print("simple v rmse =", simple_v_error) 
print("simple xyz rmse =", simple_xyz_error)

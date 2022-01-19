#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 12 09:04:47 2021

@author: guterlj
"""
from pyGITR.PostProcess import *
Post = PostProcess('/home/guterlj/simulations/GITR/data/single_microtrench_D_1_25__20000/last.npy')
Post.CumulateParticleData()
xe = Post.CumulativeData['ParticleEndData']['Data']['x']
ye = Post.CumulativeData['ParticleEndData']['Data']['y']
ze = Post.CumulativeData['ParticleEndData']['Data']['z']
xs = Post.CumulativeData['ParticleStartData']['Data']['x']
ys = Post.CumulativeData['ParticleStartData']['Data']['y']
zs = Post.CumulativeData['ParticleStartData']['Data']['z']
fig,ax = plt.subplots()
#ax.scatter(x,y,z)
plt.figure()
plt.hist2d(x, y, bins=20)
plt.figure()
L = 30e-6/2
plt.hist2d(xs, ys, bins=30,range=[[-L, L], [-L, L]],cmax=3000)
plt.colorbar()

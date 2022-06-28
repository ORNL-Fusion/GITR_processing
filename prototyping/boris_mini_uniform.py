#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 11:02:37 2021

@author: Alyssa, Jerome
"""

import numpy as np
from numpy import cross, eye, dot, sqrt, cos, sin
from numpy.linalg import norm
from scipy.linalg import expm

class particle:
    def __init__(self):
        self.x = 0
        self.y = 0
        self.z = 0


def ParticlePusher(t, m,q, v0,E,B):
    #void pusher_charged_particle(ParticleStruct *particle, Numeric *num_info)
    # /*
    #  update velocity for charged particles in E and B by solving Newton equation with Boris scheme
    #  w1 [e/amu] to be used in vmin = v- = v0 + (q/m) Edt/2
    #  w2 [e/amu] to be used in t = (q/m) Bdt/2
    #  in SI: q[C], m[kg], B[T], t[s], v[m/s], E[V/m], ...
    #  in SI: e/amu(proton) = 1.602e-19[Coulomb]/1.6605e-27[kg] = 9.6477e+07[C/kg]
    #  in ERO: v[cm/s = 1e-02 m/s] and E[V/mm = 1e+03 V/m] -> w1=(1e+05)*w2
    #  */
    
    dt = t[1]-t[0]
    
    #define unit normalization weightings
    Bmag = np.linalg.norm(B)
    w1 = dt * q / (2. * m);
    w2 = dt * q / (2. * m);
    w3 = 2 / (1 + w2**2 * Bmag**2)

    #define output history arrays
    vhist = np.zeros([len(t), len(v0)])
    v = v0
    
    ion = particle()
    spacehist = np.zeros([len(t), len(v)])
    
    for i in range(0,len(t)):
        vhist[i] = v
        
        vm_Efield = v + w1*E
        vm_Bfield = vm_Efield + w2*cross(vm_Efield, B) 
        vp_Bfield = vm_Efield + w2*w3*cross(vm_Bfield, B)
        vp_Efield = vp_Bfield + w1*E
        v = vp_Efield
        
        # track movement in space
        ion.x, ion.y, ion.z = SpacePusher(ion.x,ion.y,ion.z, vhist[i], dt)
        spacehist[i] = np.array([ion.x,ion.y,ion.z])
        
    return vhist, spacehist
        
        
def SimplePusher(t, m,q, v0,E,B):
   	# particle->vx=particle->vx+(particle->efield.Ex*100000.+particle->vy*particle->bfield.Bz-particle->vz*particle->bfield.By)*particle->dt*((double)(particle->charge)*EELEKTRON)/((double)(particle->mass)*MPROTON);
   	# particle->vy=particle->vy+(particle->efield.Ey*100000.+particle->vz*particle->bfield.Bx-particle->vx*particle->bfield.Bz)*particle->dt*((double)(particle->charge)*EELEKTRON)/((double)(particle->mass)*MPROTON);
   	# particle->vz=particle->vz+(particle->efield.Ez*100000.+particle->vx*particle->bfield.By-particle->vy*particle->bfield.Bx)*particle->dt*((double)(particle->charge)*EELEKTRON)/((double)(particle->mass)*MPROTON);
    dt = t[1]-t[0]
    vhist = np.zeros([len(t), len(v0)])
    v = v0
    
    ion = particle()
    spacehist = np.zeros([len(t), len(v)])
    
    for i in range(0,len(t)):
        vhist[i] = v
        v = v + q/m*(E+cross(v,B))*dt
        
        # track movement in space
        ion.x, ion.y, ion.z = SpacePusher(ion.x,ion.y,ion.z, vhist[i], dt)
        spacehist[i] = np.array([ion.x,ion.y,ion.z])
        
    return vhist, spacehist


def rmse(predictions,targets):
    return np.sqrt(np.mean((predictions-targets)**2))


def rotate_matrix(axis, theta):
    return expm(cross(eye(3), axis/norm(axis)*theta))


def SpacePusher(x0,y0,z0,v0,dt):
    xx = x0 + v0[0]*dt
    yy = y0 + v0[1]*dt
    zz = z0 + v0[2]*dt
    return xx, yy, zz

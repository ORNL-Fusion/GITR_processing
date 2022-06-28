# -*- coding: utf-8 -*-
"""
author: Alyssa Hayes
"""

#%% Comparing Data from NETCDF output files and HDF5 output files

def nc_initialize():
    import numpy as np
    from netCDF4 import Dataset
    example = 'micro-trenches'
    return np, Dataset, example

def ComparePositionsData():
    np, Dataset, example = nc_initialize()
    ncFile = '../pyGITR/examples/' + example + '/output/positions.nc'
    ncData = Dataset(ncFile, "r", format="NETCDF4")
    
    nc_x = ncData.variables['x'][:]
    nc_y = ncData.variables['y'][:]
    nc_z = ncData.variables['z'][:]
    
    nc_vx = ncData.variables['vx'][:]
    nc_vy = ncData.variables['vy'][:]
    nc_vz = ncData.variables['vz'][:]
    
    nc_tT = ncData.variables['transitTime'][:]
    nc_hW = ncData.variables['hitWall'][:]
    nc_sH = ncData.variables['surfaceHit'][:]
    
    nc_w = ncData.variables['weight'][:]
    nc_q = ncData.variables['charge'][:]
    nc_hL = ncData.variables['hasLeaked'][:]
    nc_dT = ncData.variables['distTraveled'][:]
    
    return

def CompareSurfaceData():
    np, Dataset, example = nc_initialize()
    ncFile = '../pyGITR/examples/' + example + '/output/surface.nc'
    ncData = Dataset(ncFile, "r", format="NETCDF4")
    
    nc_GD = ncData.variables['grossDeposition'][:]
    nc_GE = ncData.variables['grossErosion'][:]   
    nc_aveSpyl = ncData.variables['aveSpyl'][:]
    nc_spylCounts = ncData.variables['spylCounts'][:]
    
    nc_sN = ncData.variables['surfaceNumber'][:]
    nc_sumPS = ncData.variables['sumParticlesStrike'][:]
    nc_sumWS = ncData.variables['sumWeightStrike'][:]
    
    nc_EDist = ncData.variables['surfEDist'][:]
    nc_ReflDist = ncData.variables['surfReflDist'][:]
    nc_SputDist = ncData.variables['surfSputtDist'][:]
    
    #nc_gridE = ncData.variables['gridE'][:]
    #nc_gridA = ncData.variables['gridA'][:]
    
    return

def CompareHistoryData():
    np, Dataset, example = nc_initialize()
    ncFile = '../pyGITR/examples/' + example + '/output/history.nc'
    ncData = Dataset(ncFile, "r", format="NETCDF4")

    nc_x = ncData.variables['x'][:]
    nc_y = ncData.variables['y'][:]
    nc_z = ncData.variables['z'][:]
    
    nc_v = ncData.variables['v'][:]
    nc_vx = ncData.variables['vx'][:]
    nc_vy = ncData.variables['vy'][:]
    nc_vz = ncData.variables['vz'][:]
    
    nc_q = ncData.variables['charge'][:]
    nc_w = ncData.variables['weight'][:]
    
    return

if __name__ == '__main__':
    ComparePositionsData()
    CompareSurfaceData()
    CompareHistoryData()

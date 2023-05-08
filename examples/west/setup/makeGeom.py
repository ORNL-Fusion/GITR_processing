# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 15:47:47 2023

@author: elicr
"""
import sys
import os
sys.path.insert(0, '../../../python/')

import shutil
import numpy as np
import netCDF4
import matplotlib.pyplot as plt
import io, shutil, libconf

def init():
    #set plotting style defaults
    plt.rcParams.update({'font.size':11.5})
    plt.rcParams.update({'lines.linewidth':1.2})
    plt.rcParams.update({'lines.markersize':1})

def replace_line_segment(x_priority, y_priority, x_base, y_base):
    x_final = x_base;
    y_final = y_base;

    all_close = np.array(0);
    
    for i in range(len(x_priority) - 1):
        distances = np.sqrt(np.power(x_priority[i] - x_base, 2) +
                    np.power(y_priority[i] - y_base, 2));

        condition = (distances < np.sqrt(np.power((x_priority[i] - x_priority[i + 1]), 2) +
                                         np.power((y_priority[i] - y_priority[i + 1]), 2)));
        
        close_ones = np.where(condition)[0];
        all_close = np.append(all_close,close_ones);

    all_close = np.delete(all_close,0)
    remove_indices = np.unique(all_close);
    #print('remove indices',remove_indices)
    d1 = np.sqrt(
            np.power((x_priority[0] - x_base[remove_indices[0]]), 2) +
            np.power((y_priority[0] - y_base[remove_indices[0]]),2));

    d2 = np.sqrt(
            np.power((x_priority[0] - x_base[remove_indices[-1]]), 2) +
            np.power((y_priority[0] - y_base[remove_indices[-1]]), 2));

    if (d2 < d1):
        x_priority = np.flip(x_priority, 0)
        y_priority = np.flip(y_priority, 0)

    #print('x_priority',x_priority)

    # Insert the prioritized points into the proper place in the array
    x_final = np.append(x_base[0:(remove_indices[0] - 1)], x_priority)
    x_final = np.append(x_final, x_base[(remove_indices[-1] + 1):]);
    y_final = np.append(y_base[0:(remove_indices[0] - 1)], y_priority)
    y_final = np.append(y_final, y_base[(remove_indices[-1] + 1):]);

    return x_final, y_final

def gitr_lines_from_points(r,z):

    nPoints = len(r)-1;
    lines = np.zeros([nPoints, 7]);
    lines[:, 0] = r[:-1];
    lines[:, 1] = z[:-1];
    lines[:, 2] = r[1:];
    lines[:, 3] = z[1:];

    tol = 1e12;

    for i in range(nPoints):
        if (lines[i, 3] - lines[i, 1]) == 0:
            lines[i, 4] = 0;
            lines[i, 5] = lines[i, 1];
        elif ((lines[i, 2] - lines[i, 0]) == 0):
            lines[i, 4] = np.sign(lines[i, 3] - lines[i, 1]) * tol;
            lines[i, 5] = tol;

        else:
            lines[i, 4] = (lines[i, 3] - lines[i, 1]) / (lines[i, 2] - lines[i, 0]);
            lines[i, 5] = -lines[i, 4] * lines[i, 0] + lines[i, 1];

    lines[:, 6] = np.sqrt((lines[:, 2] - lines[:, 0])**2 + (lines[:, 3] - lines[:, 1])** 2);

    return lines

def lines_to_vectors(lines, inDir, filename, plt):

    x1 = lines[:, 0]
    z1 = lines[:, 1]
    x2 = lines[:, 2]
    z2 = lines[:, 3]
    slope = lines[:, 4]

    for i in range(0,len(x1)):
        if slope[i]==0: 
            perpSlope = 1.0e12
        else:
            perpSlope = -np.sign(slope[i])/np.abs(slope[i]);

        rPerp = -inDir[i]/np.sqrt(perpSlope*perpSlope+1);
        zPerp = -inDir[i]*np.sign(perpSlope)*np.sqrt(1-rPerp*rPerp);
        plt.quiver([x1[i] + (x2[i]-x1[i])/2], [z1[i] + (z2[i]-z1[i])/2], [rPerp/10], [zPerp/10], width=0.0015, scale=5, headwidth=4)

    plt.title(filename)

    
def lines_to_gitr_geometry(filename, lines, Z, surface, inDir):

    x1 = lines[:, 0]
    z1 = lines[:, 1]
    x2 = lines[:, 2]
    z2 = lines[:, 3]
    slope = lines[:, 4]
    intercept = lines[:, 5]
    line_length = lines[:, 6]
    fileExists = os.path.exists(filename)

    if not fileExists:
        f = open(filename, "w")
        f.close()

    with io.open(filename) as f:
        config = libconf.load(f)

    config['geom'] = {}
    config['geom']['x1'] = x1.tolist()
    config['geom']['z1'] = z1.tolist()
    config['geom']['x2'] = x2.tolist()
    config['geom']['z2'] = z2.tolist()
    config['geom']['slope'] = ['%.6f' % elem for elem in slope.tolist()]
    config['geom']['intercept'] = ['%.6f' % elem for elem in intercept.tolist()]
    config['geom']['length'] = line_length.tolist()
    config['geom']['Z'] = Z.tolist()
    config['geom']['surface'] = ['%i' % elem for elem in surface.tolist()]
    config['geom']['inDir'] = ['%i' % elem for elem in inDir.tolist()]
    config['geom']['y1'] = 0.0
    config['geom']['y2'] = 0.0
    config['geom']['periodic'] = 0

    with io.open(filename, 'w') as f:
        libconf.dump(config, f)

def removeQuotes(infile='this.cfg',outfile='that.cfg'):
    with open(infile, 'r') as f, open(outfile, 'w') as fo:
        for line in f:
            fo.write(line.replace('"', ''))
            
def get_target_coordinates(solps_geometry_filename):

    # Get number of mesh elements in x and y (SOLPS coordinates), nx, ny.
    # As well as the coordinates of the corners, crx, cry, 
    # and the region number from solps_geometry_filename
    nx, ny, crx, cry, region = read_b2f_geometry(solps_geometry_filename)
    print('Number of SOLPS gridpoints in x: ', nx)
    print('Number of SOLPS gridpoints in y: ', ny)

    bottom_left = 0
    bottom_right = 1

    r_inner_target = crx[0,1:,bottom_right]
    z_inner_target = cry[0,1:,bottom_right]
    r_outer_target = crx[-1,1:,bottom_left]
    z_outer_target = cry[-1,1:,bottom_left]
    print('Number of inner and outer target points (should be ny-1): ',r_inner_target.size)
    # Therefore there should be ny-2 line segments from the solps mesh to
    # be introduced to the gitr geometry

    return r_inner_target,z_inner_target, \
           r_outer_target,z_outer_target

def read_b2f_variable(solps_geometry_filename, \
                      field_name = 'crx'):
    f = open(solps_geometry_filename, 'r')
    txt = f.readlines()
    f.close()

    field_start = 0
    field_end = 0
    found = 0

    for count, line in enumerate(txt):
        if found == 0:
            if '*cf' in line:
                words = line.split()
                if words[-1] == field_name:
                    field_start = count+1
                    found = 1;
        elif found == 1:
            if '*cf' in line:
                field_end = count
                found = 2
        elif found == 2:
            break

    field = [];
    txt_list = txt[field_start:field_end]
    for sublist in txt_list:
        split_sublist = sublist.split()
        for element in split_sublist:
            field.append(element)

    field = np.array(field)
    field = field.astype(np.float)

    return field
        
def read_b2f_geometry(solps_geometry_filename):
    nxny = read_b2f_variable(solps_geometry_filename= solps_geometry_filename, \
                            field_name='nx,ny')
    nx = int(nxny[0]+2)
    ny = int(nxny[1]+2)
    #print('nx,ny',nx,ny)
    crx = np.zeros((nx, ny, 4))
    cry = np.zeros((nx, ny, 4))

    crx_long = read_b2f_variable(solps_geometry_filename= solps_geometry_filename, \
                            field_name='crx')
    cry_long = read_b2f_variable(solps_geometry_filename= solps_geometry_filename, \
                            field_name='cry')

    for i in range(4):
        crx[:,:,i] = np.transpose(np.reshape(crx_long[i*nx*ny:(i+1)*nx*ny],(ny,nx)))
        cry[:,:,i] = np.transpose(np.reshape(cry_long[i*nx*ny:(i+1)*nx*ny],(ny,nx)))

    #print('crx shape',crx.shape)
    region = read_b2f_variable(solps_geometry_filename, \
                            field_name='region')
    #print('firstwhatever',region[0:nx*ny])
    region = np.transpose(np.reshape(region[0:nx*ny],(ny,nx)))

    return nx,ny,crx,cry,region

def intersection(a_x, a_y, b_x, b_y):
    i_a = np.array(0)
    i_b = np.array(0)

    for i in range(len(b_x)):
        for j in range(len(a_x)):
            if (b_x[i] == a_x[j] and b_y[i] == a_y[j]):
                i_a = np.append(i_a, j)
                i_b = np.append(i_b, i)

    i_a = np.delete(i_a, 0)
    i_b = np.delete(i_b, 0)

    return i_a, i_b

def remove_endline_after_comma(infile='this.cfg',outfile='that.cfg'):
    with open(infile, 'r') as f, open(outfile, 'w') as fo:
        for line in f:
            fo.write(line.replace(',\n', ','))

def remove_endline_after_comma2(infile='this.cfg',outfile='that.cfg'):
    
    with open(infile, 'r') as f, open(outfile, 'w') as fo:
        for line in f:
            fo.write(line.replace(',       ', ','))
            
def add_points_divertor(r_final,z_final,W_indicesCoarse,numAddedPoints):
    #set number of added points in a line segment to be proportional to the length of the segment
    rSurfCoarse = r_final[W_indicesCoarse]
    zSurfCoarse = z_final[W_indicesCoarse]
    dist = np.sqrt((rSurfCoarse[:-1]-rSurfCoarse[1:])**2 + (zSurfCoarse[:-1]-zSurfCoarse[1:])**2)
    totalDist = np.sum(dist)
    addedPoints = numAddedPoints*dist/totalDist
    print('addedPoints',addedPoints)
    for i,v in enumerate(addedPoints): addedPoints[i] = round(v)
    addedPoints = np.array(addedPoints,dtype='int')
    
    #populate rSurfFine and zSurfFine with the added points
    rSurfFine = rSurfCoarse[0]
    zSurfFine = zSurfCoarse[0]
    for i,v in enumerate(addedPoints):
        rBegin = rSurfCoarse[i]
        zBegin = zSurfCoarse[i]
        rEnd = rSurfCoarse[i+1]
        zEnd = zSurfCoarse[i+1]
        dr = (rEnd-rBegin)/(v+1)
        dz = (zEnd-zBegin)/(v+1)
        for j in range(1,v+1):
            rSurfFine = np.append(rSurfFine,rSurfCoarse[i]+j*dr)
            zSurfFine = np.append(zSurfFine,zSurfCoarse[i]+j*dz)
        
        rSurfFine = np.append(rSurfFine,rSurfCoarse[i+1])
        zSurfFine = np.append(zSurfFine,zSurfCoarse[i+1])
    
    r_final, z_final = replace_line_segment(rSurfFine, zSurfFine, r_final, z_final)
    W_indices = np.array(range(W_indicesCoarse[0], W_indicesCoarse[-1]+numAddedPoints))
    return r_final,z_final,W_indices

def refine_target(rSurfCoarse, zSurfCoarse, rmrsCoarse, numAddedPoints=100):
    dist = np.sqrt((rSurfCoarse[:-1]-rSurfCoarse[1:])**2 + (zSurfCoarse[:-1]-zSurfCoarse[1:])**2)
    totalDist = np.sum(dist)
    addedPoints = numAddedPoints*dist/totalDist
    print('addedPoints',addedPoints)
    for i,v in enumerate(addedPoints): addedPoints[i] = round(v)
    addedPoints = np.array(addedPoints,dtype='int')
    
    #populate rSurfFine and zSurfFine with the added points
    rSurfFine = rSurfCoarse[0]
    zSurfFine = zSurfCoarse[0]
    rmrsFine = rmrsCoarse[0]
    for i,v in enumerate(addedPoints):
        rBegin = rSurfCoarse[i]
        zBegin = zSurfCoarse[i]
        rmrsBegin = rmrsCoarse[i]
        rEnd = rSurfCoarse[i+1]
        zEnd = zSurfCoarse[i+1]
        rmrsEnd = rmrsCoarse[i+1]
        dr = (rEnd-rBegin)/(v+1)
        dz = (zEnd-zBegin)/(v+1)
        dd = (rmrsEnd-rmrsBegin)/(v+1)
        for j in range(1,v+1):
            rSurfFine = np.append(rSurfFine,rSurfCoarse[i]+j*dr)
            zSurfFine = np.append(zSurfFine,zSurfCoarse[i]+j*dz)
            rmrsFine = np.append(rmrsFine,rmrsCoarse[i]+j*dd)
        
        rSurfFine = np.append(rSurfFine,rSurfCoarse[i+1])
        zSurfFine = np.append(zSurfFine,zSurfCoarse[i+1])
        rmrsFine = np.append(rmrsFine,rmrsCoarse[i+1])
    
    return rSurfFine, zSurfFine, rmrsFine

def main(gitr_geometry_filename='gitr_geometry.cfg', \
            solps_rz = 'assets/solps_rz.txt', \
            gitr_rz = 'assets/gitr_rz.txt', \
            surfW = np.append(np.arange(99,136),np.arange(139,175)), \
            solps_mesh_extra='assets/mesh.extra', \
            profiles_filename = '../input/plasmaProfiles.nc', \
            rmrs_fine_file = 'assets/rmrs_fine.txt', \
            W_fine_file = 'assets/W_fine.txt', \
            numAddedPoints = 20, \
            solps_geom = 'assets/b2fgmtry'):
    # This program uses the solps-west mesh.extra file in combination
    # with the inner and outer (left and right) divertor target
    # coordinates which come from the solps-west-data interpolation
    # program to create a 2d geometry for GITR in which
    # the solps plasma profiles properly match the divertor target
    # geometry.
    # 4
    # This geometry is then written to a config (cfg) file for
    # use in GITR simulation
    
    if numAddedPoints >= 4 and numAddedPoints <= 10:
        sys.exit('ERROR: numAddedPoints must be greater than 10 or less than 4')

    #get geometry from solps
    solps_mesh = np.loadtxt(solps_mesh_extra)

    r = solps_mesh[:, [0,2]].transpose()
    z = solps_mesh[:, [1,3]].transpose()
    
    #save (r_ogr, z_ogr) to a file for easy visualization using viz_geom_west.m
    with open(solps_rz, 'w') as f:
        for i in range(0,len(r[0,:])):
            f.write(str(r[0,i]) +' '+ str(z[0,i]) +'\n')

    #order line segments with other code
    manual_indices = np.array(range(3, 5))
    manual_indices = np.append(manual_indices, 77)
    manual_indices = np.append(manual_indices, range(5, 22))
    manual_indices = np.append(manual_indices, 78)
    manual_indices = np.append(manual_indices, range(22, 24))
    manual_indices = np.append(manual_indices, range(79, 90))
    manual_indices = np.append(manual_indices, 91)
    manual_indices = np.append(manual_indices, range(24, 50))
    manual_indices = np.append(manual_indices, 90)
    manual_indices = np.append(manual_indices, range(50, 62))
    manual_indices = np.append(manual_indices, 0)
    manual_indices = np.append(manual_indices, 62)
    manual_indices = np.append(manual_indices, 1)
    manual_indices = np.append(manual_indices, range(63, 77))
    manual_indices = np.append(manual_indices, 2)
    manual_indices = np.append(manual_indices, range(92, 102))
    manual_indices = np.append(manual_indices, range(112, 114))
    manual_indices = np.append(manual_indices, range(102, 108))
    manual_indices = np.append(manual_indices, 110)
    manual_indices = np.append(manual_indices, range(108, 110))
    manual_indices = np.append(manual_indices, 111)
    manual_indices = np.append(manual_indices, 3)

    r_west = solps_mesh[:, [0, 2]].transpose()[0, manual_indices]
    z_west = solps_mesh[:, [1, 3]].transpose()[0, manual_indices]

    # rewriting coordinates for future visualization
    with open(solps_rz, 'w') as f:
        for i in range(0,len(r_west)):
            f.write(str(r_west[i]) +' '+ str(z_west[i]) +'\n')
        
    #get target geometries from solps
    print('manual geometry size',r_west.size)
    r_right_target,z_right_target,r_left_target,z_left_target = get_target_coordinates(solps_geom)
    
    # stitch to wall geometry
    profiles = netCDF4.Dataset(profiles_filename)
    r_right_target = profiles.variables['r_inner_target'][:]
    z_right_target = profiles.variables['z_inner_target'][:]
    r_left_target = profiles.variables['r_outer_target'][:]
    z_left_target = profiles.variables['z_outer_target'][:]
    
    r_final, z_final = replace_line_segment(r_left_target, z_left_target, r_west, z_west)
    r_final, z_final = replace_line_segment(r_right_target, z_right_target, r_final, z_final)


    # writing coordinates with targets for future visualization
    with open(gitr_rz, 'w') as f:
        for i in range(0,len(r_final)):
            f.write(str(r_final[i]) +' '+ str(z_final[i]) +'\n')
            

    plt.close()
    plt.plot(r_final, z_final, linewidth=0.1)
    plt.scatter(r_final, z_final, s=0.4)
    #plt.scatter(r_west, z_west, s=0.3)
    plt.xlabel('r')
    plt.ylabel('z')
    plt.title('Target Geometry Integrated with WEST')
    plt.show()
    # plt.savefig('final_west.png')

    ###################################################################
    # Increase Fineness of W Divertor Surface
    ###################################################################
    
    rmrsCoarse_in = profiles.variables['rmrs_inner_target'][:]
    rmrsCoarse_out = profiles.variables['rmrs_outer_target'][:]

    # 1 represents outer(left) and 2 represents inner(right)
    W_indicesCoarse1 = np.array(range(99,136))
    W_indicesCoarse2 = np.array(range(139+numAddedPoints,175+numAddedPoints))
    
    r_final,z_final,W_indices1 = add_points_divertor(r_final,z_final,W_indicesCoarse1,numAddedPoints)
    r_final,z_final,W_indices2 = add_points_divertor(r_final,z_final,W_indicesCoarse2,numAddedPoints)
    W_indicesCoarse2 = np.array(range(139,175))

    #set number of added points in a line segment to be proportional to the length of the segment
    rSurfCoarse1 = r_final[W_indicesCoarse1]
    zSurfCoarse1 = z_final[W_indicesCoarse1]
    rSurfCoarse2 = r_final[W_indicesCoarse2]
    zSurfCoarse2 = z_final[W_indicesCoarse2]

    rSurfFine1, zSurfFine1, rmrsFine1 = refine_target(rSurfCoarse1,zSurfCoarse1,rmrsCoarse_out,numAddedPoints)
    rSurfFine2, zSurfFine2, rmrsFine2 = refine_target(rSurfCoarse2,zSurfCoarse2,rmrsCoarse_in,numAddedPoints)
    
    rmrsMid1 = (rmrsFine1[:-1]+rmrsFine1[1:])/2
    rmrsMid2 = (rmrsFine2[:-1]+rmrsFine2[1:])/2
    rmrsMid = np.append(rmrsMid1,rmrsMid2)

    r_final, z_final = replace_line_segment(rSurfFine1, zSurfFine1, r_final, z_final)
    r_final, z_final = replace_line_segment(rSurfFine2, zSurfFine2, r_final, z_final)
    
    W_indicesCoarse2 = np.array(range(139+numAddedPoints,175+numAddedPoints))
    W_indices = np.array(range(W_indicesCoarse1[0], W_indicesCoarse1[-1]+numAddedPoints))
    W_indices = np.append(W_indices,np.array(range(W_indicesCoarse2[0], W_indicesCoarse2[-1]+numAddedPoints)))

    
    #define interior side of each line segment in the geometry with inDir
    inDir = np.ones(len(r_final))
    inDir[0:2] = inDir[3] = inDir[11:13] = inDir[14] = inDir[16] = inDir[18:21] = \
            inDir[34:37] = inDir[58:61] = inDir[74] = inDir[77:] = -1

    #populate lines and check that vectors point inward
    lines = gitr_lines_from_points(r_final, z_final)
    lines_to_vectors(lines, inDir, 'inDir',plt)
    plt.plot(r_final,z_final)
    plt.savefig('plots/inDir.png')

    #give the divertor target segments, targ_indices, a material and an interactive surface
    Z = np.zeros(len(r_final)+1)
    surfaces = np.zeros(len(r_final)+1)

    
    Z[W_indices] = 74;
    surfaces[W_indices] = 1;
    
    #save (r_final, z_final) and W indices to a file for pulling into the particle source
    with open(gitr_rz, 'w') as f:
        for i in range(0,len(r_final)):
            f.write(str(r_final[i]) +' '+ str(z_final[i]) +'\n')
    
    with open(W_fine_file, 'w') as f:
        for i in range(0,len(W_indices)):
            f.write(str(W_indices[i])+'\n')

    with open(rmrs_fine_file, 'w') as f:
        for i in range(0,len(rmrsMid)):
            f.write(str(rmrsMid[i]) +'\n')
    
    #populate geometry input file to GITR
    lines_to_gitr_geometry(gitr_geometry_filename+'0', lines, Z, surfaces, inDir)
    removeQuotes(gitr_geometry_filename+'0', gitr_geometry_filename)
    
    #gitr.remove_endline_after_comma(infile=gitr_geometry_filename+"0", outfile=gitr_geometry_filename+"00")
    #gitr.remove_endline_after_comma2(infile=gitr_geometry_filename+"00", outfile=gitr_geometry_filename)
    
    
if __name__ == "__main__":
    main()
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 15:24:52 2023

@author: elicr
"""

# This code was used to get the r and z values into a txt file 
# so that they could be used for post processing.

import matplotlib.pyplot as plt
# sys.path.append('/home/tqd/code/netcdf4-python')
import netCDF4
import numpy as np
# import cv2
import io, libconf
# from mpl_toolkits.mplot3d import Axes3D
# from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import math
import os

def make_gitr_geometry_from_solps_west(gitr_geometry_filename='gitrGeometry.cfg', \
                                  solps_mesh_extra='/assets/mesh.extra', \
                                  solps_geom = '/assets/b2fgmtry'):
    # This program uses the solps-west mesh.extra file in combination
    # with the inner and outer (left and right) divertor target
    # coordinates which come from the solps-west-data interpolation
    # program to create a 2d geometry for GITR in which
    # the solps plasma profiles properly match the divertor target
    # geometry.
    #
    # This geometry is then written to a config (cfg) file for
    # use in GITR simulation

    #get geometry from solps
    solps_mesh = np.loadtxt(solps_mesh_extra)

    r = solps_mesh[:, [0,2]].transpose()
    z = solps_mesh[:, [1,3]].transpose()

    #order line segments
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

    plt.plot(r_west, z_west)
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.title('Raw WEST Geometry from SOLPS')
    #plt.savefig('mesh_extra_west.pdf')

    #plt.scatter(r,z,s=0.4)
    #plt.savefig('mesh_extra_west_scatter.png')


    #get target geometries from solps
    #print('manual geometry size',r_west.size)
    r_left_target,z_left_target,r_right_target,z_right_target = get_target_coordinates(solps_geom)
    plt.plot(r_left_target, z_left_target)
    plt.plot(r_right_target, z_right_target)
    plt.title('Raw WEST Targets from SOLPS')
    #plt.savefig('targets_west.pdf')

    #integrate target geometry into base geometry
    #uncomment print statements here and in replace_line_segments_west
    #to help solve errors integrating targets into the base geometry
    #print('START r_west size: ', r_west.size)
    #print('ADD r_inner_target size: ', r_left_target.size)
    r_final, z_final = replace_line_segment_west(r_left_target, z_left_target, r_west, z_west)
    #print('CHECK r_final size after replacing inner target: ', r_final.size)
    #print('ADD r_outer_target size: ', r_right_target.size)
    r_final, z_final = replace_line_segment_west(r_right_target, z_right_target, r_final, z_final)
    #print('CHECK r_final size after replacing outer target: ', r_final.size)

    plt.close()
    plt.plot(r_final, z_final, linewidth=0.1)
    plt.scatter(r_final, z_final, s=0.4)
    #plt.scatter(r_west, z_west, s=0.3)
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.title('Target Geometry Integrated with WEST')
    #plt.savefig('final_west.pdf')


    #define interior side of each line segment in the geometry with inDir
    inDir = np.ones(len(r_final))
    inDir[0:2] = inDir[3] = inDir[11:13] = inDir[14] = inDir[16] = inDir[18:21] = \
            inDir[34:37] = inDir[58:61] = inDir[74] = inDir[77:] = -1

    #populate lines and check that vectors point inward
    lines = gitr_lines_from_points(r_final, z_final)
    lines_to_vectors(lines, inDir, 'vectors_west.pdf')


    Z = np.zeros(len(r_final)+1)
    surfaces = np.zeros(len(r_final)+1)

    i_a, i_b = intersection(r_final, z_final, r_left_target, z_left_target)
    Z[i_b] = 74;
    surfaces[i_b] = 1;

    i_a, i_b = intersection(r_final, z_final, r_right_target, z_right_target)
    Z[i_b] = 74;
    surfaces[i_b] = 1;

    lines_to_gitr_geometry(gitr_geometry_filename, lines, Z, surfaces, inDir)

    removeQuotes(infile=gitr_geometry_filename, outfile=gitr_geometry_filename+"0")

    remove_endline_after_comma(infile=gitr_geometry_filename+"0", outfile=gitr_geometry_filename+"00")
    remove_endline_after_comma2(infile=gitr_geometry_filename+"00", outfile=gitr_geometry_filename)
    return r_final,z_final

def make_gitr_geometry_from_solps_west(gitr_geometry_filename='gitr_geometry.cfg', \
                                  solps_mesh_extra='/assets/mesh.extra', \
                                  solps_geom = '/assets/b2fgmtry'):
    # This program uses the solps-west mesh.extra file in combination
    # with the inner and outer (left and right) divertor target
    # coordinates which come from the solps-west-data interpolation
    # program to create a 2d geometry for GITR in which
    # the solps plasma profiles properly match the divertor target
    # geometry.
    #
    # This geometry is then written to a config (cfg) file for
    # use in GITR simulation

    #get geometry from solps
    solps_mesh = np.loadtxt(solps_mesh_extra)

    r = solps_mesh[:, [0,2]].transpose()
    z = solps_mesh[:, [1,3]].transpose()

    #order line segments
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

    plt.plot(r_west, z_west)
    plt.xlabel('r')
    plt.ylabel('z')
    plt.title('Raw WEST Geometry from SOLPS')
    plt.savefig('mesh_extra_west.png')

    #plt.scatter(r,z,s=0.4)
    #plt.savefig('mesh_extra_west_scatter.png')


    #get target geometries from solps
    #print('manual geometry size',r_west.size)
    r_left_target,z_left_target,r_right_target,z_right_target = get_target_coordinates(solps_geom)
    plt.plot(r_left_target, z_left_target)
    plt.plot(r_right_target, z_right_target)
    plt.title('Raw WEST Targets from SOLPS')
    plt.savefig('targets_west.png')

    #integrate target geometry into base geometry
    #uncomment print statements here and in replace_line_segments_west
    #to help solve errors integrating targets into the base geometry
    #print('START r_west size: ', r_west.size)
    #print('ADD r_inner_target size: ', r_left_target.size)
    r_final, z_final = replace_line_segment_west(r_left_target, z_left_target, r_west, z_west)
    #print('CHECK r_final size after replacing inner target: ', r_final.size)
    #print('ADD r_outer_target size: ', r_right_target.size)
    r_final, z_final = replace_line_segment_west(r_right_target, z_right_target, r_final, z_final)
    #print('CHECK r_final size after replacing outer target: ', r_final.size)

    plt.close()
    plt.plot(r_final, z_final, linewidth=0.1)
    plt.scatter(r_final, z_final, s=0.4)
    #plt.scatter(r_west, z_west, s=0.3)
    plt.xlabel('r')
    plt.ylabel('z')
    plt.title('Target Geometry Integrated with WEST')
    plt.savefig('final_west.png')


    #define interior side of each line segment in the geometry with inDir
    inDir = np.ones(len(r_final))
    inDir[0:2] = inDir[3] = inDir[11:13] = inDir[14] = inDir[16] = inDir[18:21] = \
            inDir[34:37] = inDir[58:61] = inDir[74] = inDir[77:] = -1

    #populate lines and check that vectors point inward
    lines = gitr_lines_from_points(r_final, z_final)
    lines_to_vectors(lines, inDir, 'vectors_west.png',plt)


    Z = np.zeros(len(r_final)+1)
    surfaces = np.zeros(len(r_final)+1)

    i_a, i_b = intersection(r_final, z_final, r_left_target, z_left_target)
    Z[i_b] = 74;
    surfaces[i_b] = 1;

    i_a, i_b = intersection(r_final, z_final, r_right_target, z_right_target)
    Z[i_b] = 74;
    surfaces[i_b] = 1;

    lines_to_gitr_geometry(gitr_geometry_filename, lines, Z, surfaces, inDir)

    removeQuotes(infile=gitr_geometry_filename, outfile=gitr_geometry_filename+"0")

    remove_endline_after_comma(infile=gitr_geometry_filename+"0", outfile=gitr_geometry_filename+"00")
    remove_endline_after_comma2(infile=gitr_geometry_filename+"00", outfile=gitr_geometry_filename)

def replace_line_segment_west(x_priority, y_priority, x_base, y_base):
    x_left_bound = min(x_priority[0], x_priority[-1])
    x_right_bound = max(x_priority[0], x_priority[-1])
    y_bottom_bound = min(y_priority[0], y_priority[-1])
    y_top_bound = max(y_priority[0], y_priority[-1])

    remove_indices = np.empty([0])

    for i in range(0, len(x_base)):
        if ((x_base[i]>x_left_bound) and (x_base[i]<x_right_bound)\
                and (y_base[i]>y_bottom_bound) and (y_base[i]<y_top_bound)):
                    remove_indices = np.append(remove_indices, i)

    remove_indices = remove_indices.astype(int)
    #print('SUBTRACT remove_indices: ', remove_indices.size)

    x_final = np.append(x_base[0:(remove_indices[0])], x_priority)
    x_final = np.append(x_final, x_base[(remove_indices[-1] + 1):]);
    y_final = np.append(y_base[0:(remove_indices[0])], y_priority)
    y_final = np.append(y_final, y_base[(remove_indices[-1] + 1):]);

    return x_final, y_final

def get_target_coordinates(solps_geometry_filename='/Users/tyounkin/Dissertation/ITER/mq3/solps/b2fgmtry'):

    # Get number of mesh elements in x and y (SOLPS coordinates), nx, ny.
    # As well as the coordinates of the corners, crx, cry, 
    # and the region number from solps_geometry_filename
    nx, ny, crx, cry, region = read_b2f_geometry(solps_geometry_filename)
    print('Number of SOLPS gridpoints in x: ', nx)
    print('Number of SOLPS gridpoints in y: ', ny)

    geom_shape = crx.shape
    bottom_left = 0
    bottom_right = 1
    top_left = 2;
    top_right = 3;

    r_inner_target = crx[0,1:,bottom_right]
    z_inner_target = cry[0,1:,bottom_right]
    r_outer_target = crx[-1,1:,bottom_left]
    z_outer_target = cry[-1,1:,bottom_left]
    print('Number of inner and outer target points (should be ny-1): ',r_inner_target.size)
    # Therefore there should be ny-2 line segments from the solps mesh to
    # be introduced to the gitr geometry

    return r_inner_target,z_inner_target, \
           r_outer_target,z_outer_target

def read_b2f_geometry(solps_geometry_filename='/Users/tyounkin/Dissertation/ITER/mq3/solps/b2fgmtry'):
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

def read_b2f_variable(solps_geometry_filename='/Users/tyounkin/Dissertation/ITER/mq3/solps/b2fgmtry', \
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

def gitr_lines_from_points(r,z):

    nPoints = len(r)-1;
    lines = np.zeros([nPoints, 7]);
    lines[:, 0] = r[:-1];
    lines[:, 1] = z[:-1];
    lines[:, 2] = r[1:];
    lines[:, 3] = z[1:];

    tol = 1e12;
    tol_small = 1e-12;

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
    intercept = lines[:, 5]
    line_length = lines[:, 6]

    for i in range(0,len(x1)):
        if slope[i]==0: 
            perpSlope = 1.0e12
        else:
            perpSlope = -np.sign(slope[i])/np.abs(slope[i]);

        rPerp = -inDir[i]/np.sqrt(perpSlope*perpSlope+1);
        zPerp = -inDir[i]*np.sign(perpSlope)*np.sqrt(1-rPerp*rPerp);
        plt.quiver([x1[i] + (x2[i]-x1[i])/2], [z1[i] + (z2[i]-z1[i])/2], [rPerp/10], [zPerp/10], width=0.0015, scale=5, headwidth=4)

    plt.title(filename)
    # plt.savefig('plots/'+filename+'.pdf')

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

def remove_endline_after_comma(infile='this.cfg',outfile='that.cfg'):
    with open(infile, 'r') as f, open(outfile, 'w') as fo:
        for line in f:
            fo.write(line.replace(',\n', ','))

def remove_endline_after_comma2(infile='this.cfg',outfile='that.cfg'):
    with open(infile, 'r') as f, open(outfile, 'w') as fo:
        for line in f:
            fo.write(line.replace(',       ', ','))


# This program uses the solps-west mesh.extra file in combination
# with the inner and outer (left and right) divertor target
# coordinates which come from the solps-west-data interpolation
# program to create a 2d geometry for GITR in which
# the solps plasma profiles properly match the divertor target
# geometry.
#
# This geometry is then written to a config (cfg) file for
# use in GITR simulation

gitr_geometry_filename='gitrGeometry.cfg'
solps_mesh_extra='assets/mesh.extra'
solps_geom = 'assets/b2fgmtry'

#get geometry from solps
solps_mesh = np.loadtxt(solps_mesh_extra)

r = solps_mesh[:, [0,2]].transpose()
z = solps_mesh[:, [1,3]].transpose()

#order line segments
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

plt.plot(r_west, z_west)
plt.xlabel('r [m]')
plt.ylabel('z [m]')
plt.title('Raw WEST Geometry from SOLPS')
#plt.savefig('mesh_extra_west.pdf')

#plt.scatter(r,z,s=0.4)
#plt.savefig('mesh_extra_west_scatter.png')


#get target geometries from solps
#print('manual geometry size',r_west.size)
r_left_target,z_left_target,r_right_target,z_right_target = get_target_coordinates(solps_geom)
plt.plot(r_left_target, z_left_target)
plt.plot(r_right_target, z_right_target)
plt.title('Raw WEST Targets from SOLPS')
#plt.savefig('targets_west.pdf')

#integrate target geometry into base geometry
#uncomment print statements here and in replace_line_segments_west
#to help solve errors integrating targets into the base geometry
#print('START r_west size: ', r_west.size)
#print('ADD r_inner_target size: ', r_left_target.size)
r_final, z_final = replace_line_segment_west(r_left_target, z_left_target, r_west, z_west)
#print('CHECK r_final size after replacing inner target: ', r_final.size)
#print('ADD r_outer_target size: ', r_right_target.size)
r_final, z_final = replace_line_segment_west(r_right_target, z_right_target, r_final, z_final)
#print('CHECK r_final size after replacing outer target: ', r_final.size)

plt.close()
plt.plot(r_final, z_final, linewidth=0.1)
plt.scatter(r_final, z_final, s=0.4)
#plt.scatter(r_west, z_west, s=0.3)
plt.xlabel('r [m]')
plt.ylabel('z [m]')
plt.title('Target Geometry Integrated with WEST')
#plt.savefig('final_west.pdf')


#define interior side of each line segment in the geometry with inDir
inDir = np.ones(len(r_final))
inDir[0:2] = inDir[3] = inDir[11:13] = inDir[14] = inDir[16] = inDir[18:21] = \
        inDir[34:37] = inDir[58:61] = inDir[74] = inDir[77:] = -1

#populate lines and check that vectors point inward
lines = gitr_lines_from_points(r_final, z_final)
lines_to_vectors(lines, inDir, 'vectors_west.pdf',plt)


Z = np.zeros(len(r_final)+1)
surfaces = np.zeros(len(r_final)+1)

i_a, i_b = intersection(r_final, z_final, r_left_target, z_left_target)
Z[i_b] = 74;
surfaces[i_b] = 1;

i_a, i_b = intersection(r_final, z_final, r_right_target, z_right_target)
Z[i_b] = 74;
surfaces[i_b] = 1;

lines_to_gitr_geometry(gitr_geometry_filename, lines, Z, surfaces, inDir)

removeQuotes(infile=gitr_geometry_filename, outfile=gitr_geometry_filename+"0")

remove_endline_after_comma(infile=gitr_geometry_filename+"0", outfile=gitr_geometry_filename+"00")
remove_endline_after_comma2(infile=gitr_geometry_filename+"00", outfile=gitr_geometry_filename)

# Create gitr_rz txt file for post processing

gitr_rz = 'assets/gitr_rz.txt'
with open(gitr_rz, 'w') as f:
    for i in range(0,len(r_final)):
        f.write(str(r_final[i]) +' '+ str(z_final[i]) +'\n')

















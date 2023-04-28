import sys
import os
sys.path.insert(0, '../../../python/')

import solps
import numpy as np
import matplotlib.pyplot as plt
import io, libconf
import netCDF4

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
    plt.savefig('plots/'+filename+'.pdf')
    
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

def V4e_v004(gitr_geometry_filename='gitrGeometry4v.cfg', \
                                    solps_geomfile = 'assets/sas-vw_v004.ogr', \
                                    solps_targfile = 'assets/b2fgmtry', \
                                    solps_rz = 'assets/solps_rz.txt', \
                                    gitr_rz = 'assets/gitr_rz.txt', \
                                    profiles_filename = '../input/plasmaProfiles.nc', \
                                    numAddedPoints = 100):
    
    # This program uses the solps geometry .ogr file to create a 2d geometry for GITR
    # in which the solps plasma profiles properly match the solps divertor target.
    # This geometry is then written to a config (cfg) file for use in GITR simulation.
    
    #read in ogr r,z wall geometry
    with open(solps_geomfile) as f: solps_geom = f.readlines()[1:]
    solps_geom[-1] += '\n' #need this for index counting in r,z extraction
    
    # uses np.
    r_ogr = z_ogr = np.empty(0)
    for row in solps_geom:
        rvalue = float(row.split(' ')[0])
        r_ogr = np.append(r_ogr, rvalue)
    
        zvalue = float(row.split(' ')[1][0:-1])
        z_ogr = np.append(z_ogr, zvalue)
    
    r_ogr = np.append(r_ogr, r_ogr[0])
    z_ogr = np.append(z_ogr, z_ogr[0])
    
    #save (r_ogr, z_ogr) to a file for easy visualization using viz_geom_sasvw.m
    with open(solps_rz, 'w') as f:
        for i in range(0,len(r_ogr)):
            f.write(str(r_ogr[i]) +' '+ str(z_ogr[i]) +'\n')
    
    #order line segments as determined visually using viz_geom_sasvw.m
    manual_indices = np.zeros(int(114/2), dtype=int)
    i=1
    for j in range(1,114):
        if j%2 == 0:
            manual_indices[i] = j
            i+=1
    
    manual_indices = np.append(manual_indices, range(113,114))
    
    r_wall = r_ogr[manual_indices]/1000 #mm->m
    z_wall = z_ogr[manual_indices]/1000 #mm->m
    
    #get target geometry from b2fgmtry and stitch to wall geometry
    r_right_target, z_right_target, r_left_target, z_left_target = solps.get_target_coordinates(solps_targfile)
    profiles = netCDF4.Dataset(profiles_filename)
    r_right_target = profiles.variables['r_inner_target'][:]
    z_right_target = profiles.variables['z_inner_target'][:]
    r_final, z_final = replace_line_segment(r_left_target, z_left_target, r_wall, z_wall)
    r_final, z_final = replace_line_segment(r_right_target, z_right_target, r_wall, z_wall)
    
    ###################################################################
    # Increase Fineness of W Divertor Surface
    ###################################################################
    W_indicesCoarse = np.array(range(24,38))
    
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
    
    #plot correctly-ordered line segments
    plt.rcParams['font.size'] = 14
    print('plotting correctly-ordered line segments to solps_wall.pdf')
    plt.close()
    plt.plot(r_final, z_final, linewidth=0.1, label='V1e_v004')
    plt.scatter(r_final, z_final, s=0.1)
    plt.axis('scaled')
    plt.xlabel('r [mm]')
    plt.ylabel('z [mm]')
    plt.title('DIII-D SAS-VW4 Geometry')
    plt.savefig('plots/solps_final.pdf')
    
    #define interior side of each line segment in the geometry with inDir
    inDir = np.ones(len(r_final))
    inDir[2:9] = inDir[10] = inDir[16:24] = -1
    inDir[51+numAddedPoints] = inDir[53+numAddedPoints] = inDir[58+numAddedPoints:60+numAddedPoints] = -1
    inDir[64+numAddedPoints:74+numAddedPoints] = -1
    
    #populate lines and check that vectors point inward
    lines = gitr_lines_from_points(r_final, z_final)
    lines_to_vectors(lines, inDir, 'inDir', plt)
    
    #give the divertor target segments, targ_indices, a material and an interactive surface
    Z = np.zeros(len(r_final))
    surfaces = np.zeros(len(r_final))
    
    W_indices = np.array(range(24,38+numAddedPoints))
    
    plt.close()
    plt.plot(r_right_target, z_right_target, '-k', label='Carbon', linewidth=0.5)
    plt.plot(r_final[W_indices], z_final[W_indices], 'violet', label='Tungsten', linewidth=0.6)
    plt.scatter(r_final[W_indices], z_final[W_indices], color='violet', s=8)
    plt.legend()
    plt.xticks(fontsize=11)
    plt.xlabel('r [m]')
    plt.ylabel('z [m]')
    plt.title('Upper Outer SAS-WV4 Divertor in DIII-D')
    plt.savefig('plots/W wall ID')
    
    Z[W_indices] = 74;
    surfaces[W_indices] = 1;
    
    #populate geometry input file to GITR
    lines_to_gitr_geometry(gitr_geometry_filename+'0', lines, Z, surfaces, inDir)
    removeQuotes(gitr_geometry_filename+'0', gitr_geometry_filename)
    
    #gitr.remove_endline_after_comma(infile=gitr_geometry_filename+"0", outfile=gitr_geometry_filename+"00")
    #gitr.remove_endline_after_comma2(infile=gitr_geometry_filename+"00", outfile=gitr_geometry_filename)
    
    #save (r_final, z_final) to a file for easy visualization in outputs
    with open(gitr_rz, 'w') as f:
        for i in range(0,len(r_final)):
            f.write(str(r_final[i]) +' '+ str(z_final[i]) +'\n')
    
    
    print('r_min:', min(r_final), '\nr_max:', max(r_final), '\nz_min:', min(z_final), '\nz_max:', max(z_final))
    print('created gitrGeometry4v.cfg')
    return r_final, z_final, r_final[W_indices], z_final[W_indices]

if __name__ == "__main__":
    V4e_v004()
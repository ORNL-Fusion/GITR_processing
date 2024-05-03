import sys
import os
sys.path.insert(0, '../../../python/')

import numpy as np
import matplotlib.pyplot as plt
import io, libconf, netCDF4

def init():
    #set plotting style defaults
    plt.rcParams.update({'font.size':10})
    plt.rcParams.update({'lines.linewidth':2.5})
    plt.rcParams.update({'lines.markersize':5})

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
    x_final = np.append(x_base[0:(remove_indices[0])], x_priority)
    x_final = np.append(x_final, x_base[(remove_indices[-1] + 1):]);
    y_final = np.append(y_base[0:(remove_indices[0])], y_priority)
    y_final = np.append(y_final, y_base[(remove_indices[-1] + 1):]);
    
    # removes boundary points from base x,y if they fall inside priority x,y
    # x_final = np.append(x_base[0:(remove_indices[0] - 1)], x_priority)
    # x_final = np.append(x_final, x_base[(remove_indices[-1] + 1):]);
    # y_final = np.append(y_base[0:(remove_indices[0] - 1)], y_priority)
    # y_final = np.append(y_final, y_base[(remove_indices[-1] + 1):]);

    return x_final, y_final

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

def lines_to_vectors(lines, inDir, plt, plot_variables):

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
        plt.title('inDir')
        #if plot_variables: plt.savefig('plots/geom/inDir.png')
    return
    
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

def refine_target(rSurfCoarse, zSurfCoarse, rmrsCoarse, numAddedPoints=100):
    dist = np.sqrt((rSurfCoarse[:-1]-rSurfCoarse[1:])**2 + (zSurfCoarse[:-1]-zSurfCoarse[1:])**2)
    totalDist = np.sum(dist)
    addedPoints = numAddedPoints*dist/totalDist
    #print('addedPoints:',addedPoints)
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
        
    rSurfFine = np.delete(rSurfFine,-1)
    zSurfFine = np.delete(zSurfFine,-1)
    
    return rSurfFine, zSurfFine, rmrsFine, sum(addedPoints)

def main(gitr_geometry_filename='gitrGeometry.cfg', \
             solps_geomfile = 'assets/sas-vw_v005_mod.ogr', \
             solps_targfile = 'assets/b2fgmtry', \
             profiles_file = '../input/plasmaProfiles.nc', \
             W_indices_profiles = np.arange(11,22), \
             tile_shift_indices = [1,9], \
             solps_rz = 'assets/solps_rz.txt', \
             gitr_rz = 'assets/gitr_rz.txt', \
             rmrs_fine_file = 'assets/rmrs_fine.txt', \
             W_fine_file = 'assets/W_fine.txt', \
             numAddedPoints = 100, \
             plot_variables = 0):
    
    # This program uses the solps geometry .ogr file to create a 2d geometry for GITR
    # in which the solps plasma profiles properly match the solps divertor target.
    # This geometry is then written to a config (cfg) file for use in GITR simulation.
    
    init()
    
    #read in ogr r,z wall geometry
    with open(solps_geomfile) as f: solps_geom = f.readlines()[1:]
    solps_geom[-1] += '\n' #need this for index counting in r,z extraction
    
    #read r,z from .ogr into np arrays then close the loop
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
    
    #order line segments as determined visually using vizGeom.m
    stop = len(r_ogr)
    manual_indices = np.zeros(int((stop+1)/2), dtype=int)
    i=1
    for j in range(1,stop):
        if j%2 == 0:
            manual_indices[i] = j
            i+=1
        
    r_wall = r_ogr[manual_indices]/1000 #mm->m
    z_wall = z_ogr[manual_indices]/1000 #mm->m
    
    #get target geometry from b2fgmtry and stitch to wall geometry
    profiles = netCDF4.Dataset(profiles_file)
    r_right_target = profiles.variables['r_inner_target'][:]
    z_right_target = profiles.variables['z_inner_target'][:]
    r_left_target = profiles.variables['r_outer_target'][:]
    z_left_target = profiles.variables['z_outer_target'][:]
    
    r_final, z_final = replace_line_segment(r_left_target, z_left_target, r_wall, z_wall)
    r_final, z_final = replace_line_segment(r_right_target, z_right_target, r_final, z_final)
    
    if plot_variables:
        #plot correctly-ordered line segments
        plt.close()
        plt.plot(r_final, z_final)
        plt.plot(r_left_target, z_left_target, 'green')
        plt.plot(r_right_target, z_right_target, 'red')
        plt.scatter(r_final, z_final, s=0.1)
        plt.axis('scaled')
        plt.xlabel('r [m]')
        plt.ylabel('z [m]')
        plt.title('DIII-D SAS-VW Geometry')
        plt.savefig('plots/geom/solps_final.pdf')
        plt.show(block=False)
    
    ###################################################################
    # Increase Fineness of W Divertor Surface
    ###################################################################
    
    rmrsCoarse = profiles.variables['rmrs_inner_target'][W_indices_profiles]
    
    W_indicesCoarse = np.array(range(69,80)) #80?
    print('length Coarse W_indices + numAddedPoints =',\
          len(W_indicesCoarse), '+', numAddedPoints, '=',len(W_indicesCoarse)+numAddedPoints)
    
    # find what W_indices (gitrGeometry) and W_indices_profiles (plasmaProfiles) 
    # should be along the plasmaProfiles.nc targets
    if plot_variables:
        plt.close()
        if tile_shift_indices != []:
            for i in tile_shift_indices:
                plt.axhline(y=z_right_target[W_indices_profiles][i], color='k', linestyle='dotted')
        plt.plot(r_right_target, z_right_target, '-k', label='Carbon', linewidth=5)
        plt.plot(r_final[W_indicesCoarse], z_final[W_indicesCoarse], 'violet', label='W Indices Relative to Full Wall', linewidth=3)
        plt.scatter(r_final[W_indicesCoarse], z_final[W_indicesCoarse], marker='_', color='violet')
        plt.plot(r_right_target[W_indices_profiles], z_right_target[W_indices_profiles], 'c', label='W Indices Relative to InterpValues.nc', linewidth=1)
        plt.legend()
        plt.axis('scaled')
        plt.xlabel('r [m]')
        plt.ylabel('z [m]')
        plt.title('Upper Outer SAS-VW Divertor in DIII-D \n makeGeom')
        plt.savefig('plots/geom/makeGeomCoarse.png')
        plt.show(block=False)
    
    print('\n')
    print('Vertex between Legs 1 and 2:,', \
          r_right_target[W_indices_profiles][tile_shift_indices[0]], z_right_target[W_indices_profiles][tile_shift_indices[0]])
    print('Vertex between Legs 2 and 3:,', \
          r_right_target[W_indices_profiles][tile_shift_indices[1]], z_right_target[W_indices_profiles][tile_shift_indices[1]],'\n')
    
    #set number of added points in a line segment to be proportional to the length of the segment
    rSurfCoarse = r_final[W_indicesCoarse]
    zSurfCoarse = z_final[W_indicesCoarse]

    rSurfFine, zSurfFine, rmrsFine, numAddedPoints = refine_target(rSurfCoarse,zSurfCoarse,rmrsCoarse,numAddedPoints)
    
    rmrsMid = (rmrsFine[:-1]+rmrsFine[1:])/2
    r_final, z_final = replace_line_segment(rSurfFine, zSurfFine, r_final, z_final)
    W_indices = np.array(range(W_indicesCoarse[0], W_indicesCoarse[-1]+numAddedPoints+1)) #+1 may need to be added because range function is exclusive
    print('length Fine W_indices:',len(W_indices))
    
    
    #test to check that the refined rmrsMid fall between the matching coarse rmrs midpoint values
    print('\n')
    print('TEST')
    rmrsMidCoarse = profiles.variables['rmrs_inner_target_midpoints'][W_indices_profiles]
    #print(rmrsMidCoarse)
    #print(rmrsMid)
    rmrsTestCoarse = np.ones(len(rmrsMidCoarse))
    rmrsTestFine = np.ones(len(rmrsMid))
    plt.close()
    tile_shift_indices = [1,9]
    if tile_shift_indices != []:
        for i in tile_shift_indices:
            plt.axvline(x=rmrsCoarse[i], color='k', linestyle='dotted')
    plt.scatter(rmrsMidCoarse, rmrsTestCoarse, s=5, color='orange', label='Coarse')
    plt.scatter(rmrsMid, rmrsTestFine, s=1, color='cyan', label='Fine')
    plt.title('Ones plotted rmrs midpoint values')
    plt.xlabel('rmrs a.k.a. D-Dsep [m]')
    plt.legend()
    plt.show(block=False)
    plt.close()
    
    rmrsCoarse = np.append(rmrsCoarse[:4],rmrsCoarse[5:]) #remove after debugging some pre-PSI nonsense
    tile_shift_indices = [1,8]
    if tile_shift_indices != []:
        for i in tile_shift_indices:
            plt.axvline(x=rmrsCoarse[i], color='k', linestyle='dotted')
    
    print('coarse coords', len(rmrsCoarse))
    print('coarse mid', len(rmrsMidCoarse))
    rmrsTestCoarse = np.ones(len(rmrsCoarse))
    rmrsTestFine = np.ones(len(rmrsFine))
    plt.scatter(rmrsCoarse, rmrsTestCoarse, s=20, color='orange', label='Coarse Coords')
    plt.scatter(rmrsFine, rmrsTestFine, s=1, color='cyan', label='Fine Coords')
    rmrsTestCoarse = np.ones(len(rmrsMidCoarse))
    rmrsTestFine = np.ones(len(rmrsMid))
    plt.scatter(rmrsMidCoarse[1:-1], rmrsTestCoarse[1:-1], s=5, color='magenta', label='Coarse Mids')
    #plt.scatter(rmrsMid, rmrsTestFine, s=1, color='magenta', label='Fine')
    
    plt.title('Ones plotted rmrs coord values')
    plt.xlabel('rmrs a.k.a. D-Dsep [m]')
    plt.legend()
    plt.show(block=True)
    plt.close()
    print('TEST')
    print('\n')
    
    
    #find strikepoint
    strikepoint_index = np.where(rmrsFine==0)[0]
    print('Strikepoint Coords:', rSurfFine[strikepoint_index], zSurfFine[strikepoint_index])
    
    if plot_variables:
        plt.close()
        plt.plot(r_right_target, z_right_target, '-k', label='Carbon', linewidth=0.5)
        plt.plot(r_final[W_indices], z_final[W_indices], 'violet', label='Tungsten', linewidth=0.5)
        plt.scatter(rSurfCoarse, zSurfCoarse, marker='.', s=20, color='green')
        plt.scatter(r_final[W_indices], z_final[W_indices], marker='.', s=10, color='violet')
        plt.scatter(rSurfFine[strikepoint_index], zSurfFine[strikepoint_index], label='Strikepoint', marker='x', color='k', s=150, zorder=5)
        plt.legend()
        plt.axis('scaled')
        plt.xlabel('r [m]')
        plt.ylabel('z [m]')
        plt.title('Cross Section of SAS-VW Divertor')
        plt.savefig('plots/geom/makeGeom.png')
        plt.show(block=True)
    
    #print('length of 3rd leg:',np.sqrt((rSurfCoarse[-1]-rSurfCoarse[-2])**2+(zSurfCoarse[-1]-zSurfCoarse[-2])**2))
    
    if plot_variables:
        #plot correctly-ordered line segments
        plt.close()
        plt.plot(r_final, z_final, 'k', label='Carbon')
        #plt.scatter(r_final, z_final, s=1)
        plt.plot(r_final[W_indices], z_final[W_indices], 'orchid', label='Tungsten')
        plt.axis('scaled')
        plt.xlabel('r [m]')
        plt.ylabel('z [m]')
        plt.title('Cross Section of DIII-D Geometry')
        plt.legend()
        plt.savefig('plots/geom/gitr_final.png')
        plt.show(block=False)
    
    #define interior side of each line segment in the geometry with inDir
    inDir = np.ones(len(r_final))
    inDir[:35] = inDir[38:49] = inDir[52] = inDir[59:68] = -1
    inDir[190] = inDir[193] = inDir[201] = inDir[210] = inDir[214] = inDir[223] = inDir[-5:] = -1
    
    #populate lines and check that vectors point inward
    lines = gitr_lines_from_points(r_final, z_final)
    lines_to_vectors(lines, inDir, plt, plot_variables)
    
    #give the divertor target segments, targ_indices, a material and an interactive surface
    Z = np.zeros(len(r_final))
    surfaces = np.zeros(len(r_final))
    
    Z[W_indices] = 74;
    surfaces[W_indices] = 1;
    
    #populate geometry input file to GITR
    lines_to_gitr_geometry(gitr_geometry_filename+'0', lines, Z, surfaces, inDir)
    removeQuotes(gitr_geometry_filename+'0', gitr_geometry_filename)
    
    #gitr.remove_endline_after_comma(infile=gitr_geometry_filename+"0", outfile=gitr_geometry_filename+"00")
    #gitr.remove_endline_after_comma2(infile=gitr_geometry_filename+"00", outfile=gitr_geometry_filename)
    
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
    
    print('r_min:', min(r_final), '\nr_max:', max(r_final), '\nz_min:', min(z_final), '\nz_max:', max(z_final))
    print('Created gitrGeometry.cfg')
    
    return

if __name__ == "__main__":
    main(plot_variables = 1)




import sys
# adds path to gitr.py and solps.py
sys.path.insert(0, '../../../python/')

import gitr
import solps
import numpy as np
import matplotlib.pyplot as plt

# This function uses the original geometry .ogr file to create a 2d geometry for GITR
# in which the solps divertor target replaces the .ogr divertor target geometry.
# This geometry is then written to a config (cfg) file for use in GITR simulation.

# start of the function and labeling of files
def V6e_v002(gitr_geometry_filename='gitrGeometry.cfg', \
                                  solps_geomfile = 'assets/geom-SASV6/SAS-V6e_v002.ogr', \
                                  solps_targfile = 'assets/b2fgmtry', \
                                  solps_rz = 'assets/solps_rz.txt', \
                                  gitr_rz = 'assets/gitr_rz.txt', \
                                  surf_coarse = 'assets/surf_coarse.txt', \
                                  surf_ind = 'assets/surf_ind.txt', \
                                  numAddedPoints = 100):

    ###################################################################
    # Turn your file data into coordinate points (r and z for 2-D or r, z, and t for 3-D)
    # that can be plotted as your desired shape. Assume the points are not in order
    ###################################################################
    
    # read data in ogr r,z wall geometry from original gemotry file
    with open(solps_geomfile) as f: solps_geom = f.readlines()[1:]
    solps_geom[-1] += '\n' # need this for index counting in r,z extraction

    # creates empty matricies for r_ogr and z_ogr
    r_ogr = z_ogr = np.empty(0)
    for row in solps_geom:
        rvalue = float(row.split(' ')[0])
        r_ogr = np.append(r_ogr, rvalue)

        zvalue = float(row.split(' ')[1][0:-1])
        z_ogr = np.append(z_ogr, zvalue)

    # creates a closed geometery by putting the first data point at the end
    r_ogr = np.append(r_ogr, r_ogr[0])
    z_ogr = np.append(z_ogr, z_ogr[0])

    # save (r_ogr, z_ogr) to a file for easy visualization using viz_geom_sasvw.m
    with open(solps_rz, 'w') as f:
        for i in range(0,len(r_ogr)):
            f.write(str(r_ogr[i]) +' '+ str(z_ogr[i]) +'\n')

    ###################################################################
    # Reordering Coordinate points since .ogr file is not in order
    ###################################################################
    
    # order line segments r and z as determined visually using viz_geom_sasvw.m(make into python)
    # so that they're in spatial order
    manual_indices = np.zeros(int(122/2), dtype=int)
    i=1
    for j in range(1,122):
        if j%2 == 0:
            manual_indices[i] = j
            i+=1

    manual_indices = np.append(manual_indices, range(121,160))

    # must convert points into meters once they are ordered
    r_wall = r_ogr[manual_indices]/1000 #mm->m
    z_wall = z_ogr[manual_indices]/1000 #mm->m

    #get target geometry from b2fgmtry and stitch to wall geometry
    r_right_target, z_right_target, r_left_target, z_left_target = solps.get_target_coordinates(solps_targfile)

    r_final, z_final = gitr.replace_line_segment(r_left_target, z_left_target, r_wall, z_wall)
    r_final, z_final = gitr.replace_line_segment(r_right_target, z_right_target, r_wall, z_wall)
    r_final_coarse, z_final_coarse = r_final, z_final
    
    ###################################################################
    # Increase Fineness of W Divertor Surface
    ###################################################################
    W_indicesCoarse = np.array(range(30,45))
    
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
    
    r_final, z_final = gitr.replace_line_segment(rSurfFine, zSurfFine, r_final, z_final)
    W_indices = np.array(range(W_indicesCoarse[0], W_indicesCoarse[-1]+numAddedPoints))

    #plot correctly-ordered line segments
    print('plotting correctly-ordered line segments to solps_wall.pdf')
    plt.close()
    plt.plot(r_final, z_final, linewidth=0.1, label='V6e_v002')
    plt.scatter(r_final, z_final, s=0.1)
    plt.axis('scaled')
    plt.xlabel('r [mm]')
    plt.ylabel('z [mm]')
    plt.title('DIII-D SAS-VW Geometry')
    plt.savefig('plots/solps_final.pdf')

    #define interior side of each line segment in the geometry with inDir
    inDir = np.ones(len(r_final))
    inDir[2:9] = inDir[10] = inDir[17:30] = inDir[61+numAddedPoints] = inDir[63+numAddedPoints] = \
        inDir[68+numAddedPoints:70+numAddedPoints] = inDir[74+numAddedPoints:] = -1

    #populate lines and check that vectors point inward
    lines = gitr.gitr_lines_from_points(r_final, z_final)
    gitr.lines_to_vectors(lines, inDir, 'inDir', plt)

    plt.close()
    fs = 14
    plt.plot(r_right_target, z_right_target, '-k', label='Carbon', linewidth=0.5)
    plt.plot(r_final[W_indices], z_final[W_indices], 'violet', label='Tungsten', linewidth=0.6)
    plt.scatter(r_final[W_indices], z_final[W_indices], color='violet', marker='_', s=10)
    plt.legend()
    plt.xlabel('r [m]',fontsize=fs)
    plt.ylabel('z [m]',fontsize=fs)
    plt.xticks(fontsize=fs-3)
    plt.yticks(fontsize=fs)
    plt.title('Upper Outer SAS-VW Divertor in DIII-D',fontsize=fs)
    plt.savefig('plots/W wall ID')

    #give the divertor target segments, targ_indices, a material and an interactive surface
    Z = np.zeros(len(r_final))
    surfaces = np.zeros(len(r_final))
    Z[W_indices] = 74;
    surfaces[W_indices] = 1;

    #populate geometry input file to GITR
    gitr.lines_to_gitr_geometry(gitr_geometry_filename+'0', lines, Z, surfaces, inDir)
    gitr.removeQuotes(gitr_geometry_filename+'0', gitr_geometry_filename)

    # rewrites .cfg files to avoid clutter (deprocated, clean up of gitr.py need before deletion)
    #gitr.remove_endline_after_comma(infile=gitr_geometry_filename+"0", outfile=gitr_geometry_filename+"00")
    #gitr.remove_endline_after_comma2(infile=gitr_geometry_filename+"00", outfile=gitr_geometry_filename)

    #save (r_final, z_final) to a file for easy visualization in outputs
    with open(gitr_rz, 'w') as f:
        for i in range(0,len(r_final)):
            f.write(str(r_final[i]) +' '+ str(z_final[i]) +'\n')

    with open(surf_coarse, 'w') as f:
        for i in range(0,len(W_indicesCoarse)):
            f.write(str(W_indicesCoarse[i])+'\n')

    with open(surf_ind, 'w') as f:
        for i in range(0,len(W_indices)):
            f.write(str(W_indices[i])+'\n')

    print('r_min:', min(r_final), '\nr_max:', max(r_final), '\nz_min:', min(z_final), '\nz_max:', max(z_final))
    print('created gitrGeometry.cfg')
    return r_final, z_final, r_final[W_indices], z_final[W_indices], r_final_coarse, z_final_coarse, addedPoints
if __name__ == "__main__":
    V6e_v002()

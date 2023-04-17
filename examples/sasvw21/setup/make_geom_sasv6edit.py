import sys
sys.path.insert(0, '../../../python/')

import numpy as np
import matplotlib.pyplot as plt
import solps
import gitr

# Plotting function
def plot(r,z,sizescat,sizelin,name,savename,scatter,line,legend,color,fontsize):
    if line == 1:
        plt.plot(r, z, linewidth=sizelin, label=name, color=color)
    else:
        return
    if scatter == 1:    
        plt.scatter(r, z, s=sizescat,color=color)
    else:
        return
    if legend == 1:
        plt.legend()
    plt.axis('scaled')
    plt.xlabel('r [m]',fontsize=fontsize)
    plt.ylabel('z [m]',fontsize=fontsize)
    plt.title(name,fontsize=fontsize)
    savename0 = 'plots/'+savename+'.pdf'
    plt.savefig(savename0)
    
def V6e_v002(gitr_geometry_filename='gitrGeometry.cfg', \
                                  solps_geomfile = 'assets/geom-SASV6/SAS-V6e_v002.ogr', \
                                  solps_targfile = 'assets/b2fgmtry', \
                                  solps_rz = 'assets/geom-SASV6/solps_rz.txt', \
                                  gitr_rz = 'assets/geom-SASV6/gitr_rz.txt'):

    # This program uses the solps geometry .ogr file to create a 2d geometry for GITR
    # in which the solps plasma profiles properly match the divertor target geometry.
    # This geometry is then written to a config (cfg) file for use in GITR simulation.

    #read in ogr r,z wall geometry
    with open(solps_geomfile) as f: solps_geom = f.readlines()[1:]
    solps_geom[-1] += '\n' #need this for index counting in r,z extraction

    #use np.empty instead of np.zeros
    r_ogr = z_ogr = np.zeros(1)
    r_ogr = np.delete(r_ogr,0)
    z_ogr = np.delete(z_ogr,0)
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
    manual_indices = np.zeros(int(122/2), dtype=int)
    i=1
    for j in range(1,122):
        if j%2 == 0:
            manual_indices[i] = j
            i+=1

    manual_indices = np.append(manual_indices, range(121,160))

    r_wall = r_ogr[manual_indices]/1000 #mm->m
    z_wall = z_ogr[manual_indices]/1000 #mm->m

    #get target geometry from b2fgmtry and stitch to wall geometry
    r_right_target, z_right_target, r_left_target, z_left_target = solps.get_target_coordinates(solps_targfile)

    r_final, z_final = gitr.replace_line_segment(r_left_target, z_left_target, r_wall, z_wall)
    r_final, z_final = gitr.replace_line_segment(r_right_target, z_right_target, r_wall, z_wall)

    #plot correctly-ordered line segments
    print('plotting correctly-ordered line segments to solps_wall.pdf')
    plt.close()
    plot(r_final, z_final,0.1,0.1,'DIII-D SAS-VW Geometry','W wall ID',1,1,0,'b',10)

    #define interior side of each line segment in the geometry with inDir
    inDir = np.ones(len(r_final))
    inDir[2:9] = inDir[10] = inDir[17:30] = inDir[61] = inDir[63] = inDir[68:70] = inDir[74:] = -1

    #populate lines and check that vectors point inward
    lines = gitr.gitr_lines_from_points(r_final, z_final)
    gitr.lines_to_vectors(lines, inDir, 'inDir', plt)
    
    #give the divertor target segments, targ_indices, a material and an interactive surface
    Z = np.zeros(len(r_final))
    surfaces = np.zeros(len(r_final))

    W_indices = np.array(range(30,45))

    plt.close()
    fs =14
    plot(r_right_target, z_right_target,1,0.5,'Target','W wall ID',0,1,1,'black',14)
    plot(r_final[W_indices], z_final[W_indices],8,0.6,'W','W wall ID',1,1,1,'purple',14)
    plt.xticks(fontsize=fs-3)
    plt.yticks(fontsize=fs)

    Z[W_indices] = 74;
    surfaces[W_indices] = 1;

    #populate geometry input file to GITR
    gitr.lines_to_gitr_geometry(gitr_geometry_filename+'0', lines, Z, surfaces, inDir)
    gitr.removeQuotes(gitr_geometry_filename+'0', gitr_geometry_filename)

    #gitr.remove_endline_after_comma(infile=gitr_geometry_filename+"0", outfile=gitr_geometry_filename+"00")
    #gitr.remove_endline_after_comma2(infile=gitr_geometry_filename+"00", outfile=gitr_geometry_filename)

    #save (r_final, z_final) to a file for easy visualization in outputs
    with open(gitr_rz, 'w') as f:
        for i in range(0,len(r_final)):
            f.write(str(r_final[i]) +' '+ str(z_final[i]) +'\n')


    print('r_min:', min(r_final), '\nr_max:', max(r_final), '\nz_min:', min(z_final), '\nz_max:', max(z_final))
    print('created gitrGeometry.cfg')
    return r_final, z_final, r_final[W_indices], z_final[W_indices]

if __name__ == "__main__":
    V6e_v002()

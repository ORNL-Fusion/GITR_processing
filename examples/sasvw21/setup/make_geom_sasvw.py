import sys
sys.path.insert(0, '../../../python/')

import gitr as g
import numpy as np
import matplotlib.pyplot as plt

def make_gitr_geometry_from_solps_sasvw(gitr_geometry_filename='gitrGeometry.cfg', \
                                  solps_geomfile = 'assets/geom-SASV/SAS-V6e_v002.ogr', \
                                  solps_rz = 'assets/geom-SASV/solps_rz.txt'):

    # This program uses the solps geometry .ogr file to create a 2d geometry for GITR
    # in which the solps plasma profiles properly match the divertor target geometry.
    # This geometry is then written to a config (cfg) file for use in GITR simulation.

    #read in ogr r,z wall geometry
    with open(solps_geomfile) as f: solps_geom = f.readlines()[1:]
    solps_geom[-1] += '\n' #need this for index counting in r,z extraction

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

    r = r_ogr[manual_indices]
    z = z_ogr[manual_indices]
    
    #plot correctly-ordered line segments
    print('plotting correctly-ordered line segments to solps_geom.pdf')
    plt.close()
    plt.plot(r, z, linewidth=0.1, label='V6e_v002')
    plt.scatter(r, z, s=0.1)
    plt.axis('scaled')
    plt.xlabel('r [mm]')
    plt.ylabel('z [mm]')
    plt.title('DIII-D SAS-VW Geometry')
    plt.legend()
    plt.savefig('plots/solps_geom.pdf')
    """
    #define interior side of each line segment in the geometry with inDir
    inDir = np.ones(len(r))
    inDir[0:12] = inDir[15:18] = inDir[19:41] = inDir[47] = inDir[57:61] = inDir[101:103] = inDir[106:] = -1

    #populate lines and check that vectors point inward
    lines = g.gitr_lines_from_points(r,z)
    g.lines_to_vectors(lines, inDir, 'inDir', plt)

    #give the divertor target segments, targ_indices, a material and an interactive surface
    Z = np.zeros(len(r))
    surfaces = np.zeros(len(r))

    targ_indices = np.array(range(58,82))
    Z[targ_indices] = 74;
    surfaces[targ_indices] = 1;

    #populate geometry input file to GITR
    g.lines_to_gitr_geometry(gitr_geometry_filename, lines, Z, surfaces, inDir)

    #g.removeQuotes(infile=gitr_geometry_filename, outfile=gitr_geometry_filename+"0")

    #g.remove_endline_after_comma(infile=gitr_geometry_filename+"0", outfile=gitr_geometry_filename+"00")
    #g.remove_endline_after_comma2(infile=gitr_geometry_filename+"00", outfile=gitr_geometry_filename)

    print('created gitrGeometry.cfg')
    return r/1000, z/1000 #unit correction to meters
    """
if __name__ == "__main__":
    make_gitr_geometry_from_solps_sasvw()

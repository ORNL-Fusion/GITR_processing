import sys
sys.path.insert(0, '../../../python/')

import gitr as g
import numpy as np
import matplotlib.pyplot as plt

def make_gitr_geometry_from_solps_sasvw(gitr_geometry_filename='gitrGeometry.cfg', \
                                  solps_geomfile = 'assets/vvfile'):
    # This program uses the solps geometry .ogr file to create a 2d geometry for GITR
    # in which the solps plasma profiles properly match the divertor target geometry.
    #
    # This geometry is then written to a config (cfg) file for use in GITR simulation. 

    #get geometry from solps
    solps_geom = np.loadtxt(solps_geomfile)

    #order line segments as determined visually using viz_geom_sasvw.m
    manual_indices = np.array(range(0,5))
    manual_indices = np.append(manual_indices, 6)
    manual_indices = np.append(manual_indices, 5)
    manual_indices = np.append(manual_indices, range(7,110))
    manual_indices = np.append(manual_indices, range(111,117))

    r = solps_geom.transpose()[0,manual_indices]
    z = solps_geom.transpose()[1,manual_indices]
    r1 = r[:-1]
    r2 = r[1:]
    z1 = z[:-1]
    z2 = z[1:]

    #plot correctly-ordered line segments
    plt.plot(r,z,linewidth=0.1)
    plt.scatter(r, z, s=0.4)
    plt.xlabel('r [mm]')
    plt.ylabel('z [mm]')
    plt.title('DIII-D SAS-VW Geometry')
    plt.savefig('solps_geom.pdf')
    plt.close()

    #define in-direction, inDir

    """
    #define interior side of each line segment in the geometry with inDir
    inDir = np.ones(len(r_final))
    inDir[0:2] = inDir[3] = inDir[11:13] = inDir[14] = inDir[16] = inDir[18:21] = \
            inDir[34:37] = inDir[58:61] = inDir[74] = inDir[77:] = -1

    #populate lines and check that vectors point inward
    lines = gitr_lines_from_points_west(r_final, z_final)
    lines_to_vectors_west(lines, inDir, 'vectors_west.pdf')


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
    """

if __name__ == "__main__":
    make_gitr_geometry_from_solps_sasvw()

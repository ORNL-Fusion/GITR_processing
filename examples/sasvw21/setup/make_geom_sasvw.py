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

    #plot correctly-ordered line segments
    plt.plot(r,z,linewidth=0.1)
    plt.scatter(r, z, s=0.1)
    plt.axis('scaled')
    plt.xlabel('r [mm]')
    plt.ylabel('z [mm]')
    plt.title('DIII-D SAS-VW Geometry')
    plt.savefig('plots/solps_geom.pdf')

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

if __name__ == "__main__":
    make_gitr_geometry_from_solps_sasvw()

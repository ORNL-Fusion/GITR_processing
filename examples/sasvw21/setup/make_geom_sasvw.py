import sys
sys.path.insert(0, '../../../python/')

import gitr as g
import numpy as np

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

    print(len(r), len(r1), len(r2))

    import matplotlib.pyplot as plt
    plt.plot(r,z,linewidth=0.05)
    plt.savefig('trash.pdf')

    """
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
    r_left_target,z_left_target,r_right_target,z_right_target = solps.get_target_coordinates(solps_geom)
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

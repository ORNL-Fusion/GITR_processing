#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# Geometry and meshing for micro-trenches experiments on DiMES
# Reference paper:
@author: guterl
"""
import gmsh
import math
from pyGITR.Geom import GeomSetup
import gmsh
from pyGITR.gmsh_helper import SetGroups
# Clean up gmsh session if script already executed
# (Re)initialize gsmh session

try:
    gmsh.finalize()
except:
    pass
finally:
    gmsh.initialize()

# #Define DiMES cap surface
xDiMES=0.0
yDiMES=0.0
zDiMES=0.0
rDiMES=0.025



#Define micro-trenches on DiMES cap surface
Lmesh = 1e-4
dxTrench1,dyTrench1,dzTrench1=30e-6,30e-6,-2e-6 #dimension
xTrench1,yTrench1,zTrench1=0.0,0,0 # position
xTrench1,yTrench1=xTrench1-dxTrench1/2,yTrench1-dyTrench1/2 # position
axTrench1,ayTrench1,azTrench1=0,0,1 #ax rotation
alphaTrench1=0#math.pi/4 #angle of rotation

# Define large rectanlge around
xc_rect=0.0
yc_rect=0.0
zc_rect=0.0
dx_rect=1.1*dxTrench1
dy_rect=1.1*dyTrench1

xc_rect = xc_rect - dx_rect/2
yc_rect = yc_rect - dy_rect/2
Lxbc = dx_rect/2
Lybc = dy_rect/2

# dxTrench2,dyTrench2,dzTrench2=30e-6,30e-6,-2e-6 # position
# xTrench2,yTrench2,zTrench2=0.003,0,0 #dimension
# xTrench2,yTrench2=xTrench2-dxTrench2/2,yTrench2-dyTrench2/2 #dimension

# axTrench2,ayTrench2,azTrench2=0,0,1 #ax rotation
# alphaTrench2=math.pi/4 #angle of rotation

# Mesh parameters
NpTrench = 2 #number of points on edges of the micro-trenches
NpCap = 150 #number of points on edges of Dimes cap

# Tags
TagDiMES0=1
TagDiMES=20
TagTrench1=100
TagTrench2=50
TaglDot=30
TaglDot0=900
TagRect=10

# Create DiMES surface and micro-trenches
gmsh.model.occ.addDisk(xDiMES,yDiMES,zDiMES,rDiMES,rDiMES,TagDiMES0)
gmsh.model.occ.addRectangle(xc_rect,yc_rect,zc_rect,dx_rect,dy_rect,TagRect)
print(gmsh.model.occ.getEntities(2)) # look at tags of new 2D elements
gmsh.model.occ.addBox(xTrench1, yTrench1, zTrench1, dxTrench1, dyTrench1, dzTrench1,TagTrench1)
gmsh.model.occ.rotate([(3, TagTrench1)], xTrench1, yTrench1, zTrench1, axTrench1, ayTrench1, azTrench1, alphaTrench1)
print(gmsh.model.occ.getEntities(2)) # look at tags of new 2D elements
#gmsh.model.occ.addBox(xTrench2, yTrench2, zTrench2, dxTrench2, dyTrench2, dzTrench2,TagTrench2)
#print(gmsh.model.occ.getEntities(2)) # look at tags of new 2D elements

gmsh.model.occ.cut([(2,TagRect)],[(2,16)],removeTool=False)
gmsh.model.occ.cut([(2,TaglDot)],[(3,TagTrench1)],TaglDot,removeTool=False)

# Remove the 3D objects.  2D faces of the objectS are not deleted unless recursive is set to True
gmsh.model.occ.remove([(3,TagTrench1)])
gmsh.model.occ.remove([(3,TagTrench2)])

# Remove 2D face on the top of the micro-trenches to open the micro-trenches.
gmsh.model.occ.remove([(2,16)])
gmsh.model.occ.remove([(2,112)])

# Set mesh size for the micro-trenches using vertices of the micro-trench if not using field methods (see below with box and ball)
gmsh.model.occ.mesh.setSize([(0, i) for i in [3,5,7,9]], abs(dzTrench1)/NpTrench)

# Synchronize necessary before mesh setup and generation
gmsh.model.occ.synchronize()

# Setup mesh
# Mesh size field for trench #1
# Because trench has been rotated, we have to define a wider box
gmsh.model.mesh.field.add("Box", 1)
gmsh.model.mesh.field.setNumber(1, "VIn", abs(dzTrench1)/NpTrench)
gmsh.model.mesh.field.setNumber(1, "VOut", 1)
gmsh.model.mesh.field.setNumber(1, "XMin", xTrench1-dxTrench1*2)
gmsh.model.mesh.field.setNumber(1, "XMax", xTrench1+dxTrench1*2)
gmsh.model.mesh.field.setNumber(1, "YMin", yTrench1-dyTrench1*2)
gmsh.model.mesh.field.setNumber(1, "YMax", yTrench1+dyTrench1*2)
gmsh.model.mesh.field.setNumber(1, "ZMin", -1)
gmsh.model.mesh.field.setNumber(1, "ZMax", 1)
gmsh.model.mesh.field.setNumber(1, "Thickness",abs(dzTrench1)/NpTrench*10)

# # Mesh size field for trench #2
# gmsh.model.mesh.field.add("Box", 12)
# gmsh.model.mesh.field.setNumber(12, "VIn", abs(dzTrench2)/NpTrench)
# gmsh.model.mesh.field.setNumber(12, "VOut", 1)
# gmsh.model.mesh.field.setNumber(12, "XMin", xTrench2)
# gmsh.model.mesh.field.setNumber(12, "XMax", xTrench2+dxTrench2)
# gmsh.model.mesh.field.setNumber(12, "YMin", yTrench2)
# gmsh.model.mesh.field.setNumber(12, "YMax", yTrench2+dyTrench2)
# gmsh.model.mesh.field.setNumber(12, "ZMin", -1)
# gmsh.model.mesh.field.setNumber(12, "ZMax", 1)
# gmsh.model.mesh.field.setNumber(12, "Thickness",abs(dzTrench2)/NpTrench*10)

# # Mesh size field for large dot
gmsh.model.mesh.field.add("Ball", 2)
gmsh.model.mesh.field.setNumber(2, "VIn", 0.0005)
gmsh.model.mesh.field.setNumber(2, "VOut", 0.005)
gmsh.model.mesh.field.setNumber(2, "XCenter", xlDot)
gmsh.model.mesh.field.setNumber(2, "YCenter", ylDot)
gmsh.model.mesh.field.setNumber(2, "ZCenter", zlDot)
gmsh.model.mesh.field.setNumber(2, "Radius", rlDot)
gmsh.model.mesh.field.setNumber(2, "Thickness", rlDot*0.1)

# Set mesh size fieldslist
gmsh.model.mesh.field.add("Min", 3)
gmsh.model.mesh.field.setNumbers(3, "FieldsList", [1,12,2])

# Turn off background(global) mesh size rules and set fieldslist as mesh size rules.
# gmsh.option.setNumber("Mesh.MeshSizeFromCurvature",0)
# gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
# gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
# gmsh.model.mesh.field.setAsBackgroundMesh(3)


# We can constraint the min and max element sizes to stay within reasonnable values
gmsh.option.setNumber("Mesh.MeshSizeMin", Lmesh)
gmsh.option.setNumber("Mesh.MeshSizeMax", Lmesh)

# Generate 2D mesh
mesh = gmsh.model.mesh.generate(2)

# Launch the GUI to see the results:
#gmsh.fltk.run()

# Write mesh into a meshio format
gmsh.write("single_micro-trench.msh")

# close gmsh session
gmsh.finalize()
# Load geometry
g = GeomSetup('single_micro-trench.msh', Verbose=True)

# Show existing groups
g.ShowGroups()

# Plot mesh for some groups
g.Plot()

# Zoom on the plot (the zoom function for 3D axis in matplotlib 3.3.3 is not working so it has to be done manually)
g.SetAxisLim(-Lxbc*2, Lxbc*2)

# Show normals for the boundbox
g.ShowNormals(L=0.02)

# Set properties for each group
g.SetElemAttr([],'surface', 1)


# set inDir. Empty list = all elements
g.SetElemAttr([], 'inDir',-1)

# Set Z for material
g.SetElemAttr([], 'Z', 0)

# Set potential for biasing
g.SetElemAttr([],'potential',0)
g.SetAttr('lambda',0.0003)
g.SetAttr('zref',0)
g.SetAttr('periodic_bc_x',1)
g.SetAttr('periodic_bc_y',1)
g.SetAttr('periodic_bc_x0',-Lxbc)
g.SetAttr('periodic_bc_x1',+Lxbc)
g.SetAttr('periodic_bc_y0',-Lybc)
g.SetAttr('periodic_bc_y1',+Lybc)
# Plot geometry showing values of Z with color
g.Plot(ElemAttr='Z', Alpha=0.1)


# Write the geometry file
g.WriteGeomFile(Folder='input',OverWrite=True)

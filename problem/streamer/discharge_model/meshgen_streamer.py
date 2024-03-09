# Copied from streamer_mesh2.py. Supports embedding the points in the side of the domain

import gmsh
import numpy as np
from gmsh import model as gm

#############################################
case_name = 'chen_geom'
#############################################

gmsh.initialize()
gm.add(case_name)

############## CREATE PLANE GEOMETRY ##############

gm.occ.importShapes('mesh_square.STEP')

gm.occ.synchronize()
# # print(gm.getEntities())
# gmsh.fltk.run()
# exit()

# Exraction of the bbox parameters
# xmin, ymin, zmin, xmax, ymax, zmax = gm.occ.getBoundingBox(2, 1)
# size_scale_factor = .1/(xmax-xmin)
# print(xmin)
# print(xmax)
# print(ymin)
# print(ymax)
# exit()

# Scale and rotate mesh -> This is to align the imported CAD geometry with the desired CSYS and scaling for the problem.
# dimtags = gm.occ.getEntities()
# gmsh.model.occ.dilate(dimtags, 0,0,0, size_scale_factor, size_scale_factor, size_scale_factor)

# gmsh.model.occ.rotate(dimtags, 0,0,0, 1, 0, 0, np.pi/2)
xmin, ymin, zmin, xmax, ymax, zmax = gm.occ.getBoundingBox(2, 1)    # Need to over-write the BBox values so that we can generate an accurate length scale later

internal_tags = []
internal_tags.append(gmsh.model.occ.addPoint(0, 100, 0, tag=10))
for tag in internal_tags:
    gmsh.model.occ.fragment([(0,tag)], gmsh.model.occ.getEntities())

gm.occ.synchronize()    # This is important!
# gmsh.fltk.run()
# exit()

# Uncomment for refinement
# Uniform refinement in X from wall
gm.mesh.field.add("Distance", 4)
gm.mesh.field.setNumbers(4, "CurvesList", [6])
gm.mesh.field.setNumber(4, "Sampling", 300)
gmsh.model.mesh.field.add("MathEval", 12)
gmsh.model.mesh.field.setString(12, "F", ".1+.25*(F4-6)*(Atan(1000*(F4-6))/Pi + 0.5) - Atan(1000)/Pi + 0.5")

# High refinement right at wall
gmsh.model.mesh.field.add("MathEval", 14)
gmsh.model.mesh.field.setString(14, "F", ".02+.08*F4^2")

# Lower boundary
gm.mesh.field.add("Distance", 5)
gm.mesh.field.setNumbers(5, "PointsList", [3])
gmsh.model.mesh.field.add("MathEval", 13)
gmsh.model.mesh.field.setString(13, "F", ".02+.01*F5^2")


# Use the minimum of all the fields as the background mesh field:
gm.mesh.field.add("Min", 20)
# gm.mesh.field.setNumbers(20, "FieldsList", [12,13])
gm.mesh.field.setNumbers(20, "FieldsList", [12,14])

gm.mesh.field.setAsBackgroundMesh(20)
gmsh.option.setNumber("Mesh.MeshSizeMax", 60)
gm.mesh.setSmoothing(2, 1, 10)

gm.mesh.generate(2)
gmsh.write('streamer_tmp.msh3')
gmsh.fltk.run()
gmsh.finalize()
exit()

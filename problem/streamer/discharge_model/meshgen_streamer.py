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

# refined_pts = np.linspace()
internal_tags = []
# internal_tags.append(gmsh.model.occ.addPoint(0, 98, 0, tag=10))
# internal_tags.append(gmsh.model.occ.addPoint(0, 96, 0, tag=11))
# internal_tags.append(gmsh.model.occ.addPoint(0, 94, 0, tag=12))
# internal_tags.append(gmsh.model.occ.addPoint(0, 92, 0, tag=13))
# internal_tags.append(gmsh.model.occ.addPoint(0, 90, 0, tag=14))
# internal_tags.append(gmsh.model.occ.addPoint(0, 88, 0, tag=15))
# internal_tags.append(gmsh.model.occ.addPoint(0, 86, 0, tag=16))
# internal_tags.append(gmsh.model.occ.addPoint(0, 84, 0, tag=17))
# internal_tags.append(gmsh.model.occ.addPoint(0, 82, 0, tag=18))
# internal_tags.append(gmsh.model.occ.addPoint(0, 80, 0, tag=19))
# internal_tags.append(gmsh.model.occ.addPoint(0, 78, 0, tag=20))
# gmsh.model.mesh.embed(0, internal_tags,1, 2)       DON'T USE THIS

internal_tags.append(gmsh.model.occ.addPoint(0, 100, 0, tag=10))
# internal_tags.append(gmsh.model.occ.addPoint(0, 40, 0, tag=11))
for tag in internal_tags:
    gmsh.model.occ.fragment([(0,tag)], gmsh.model.occ.getEntities())

gm.occ.synchronize()    # This is important!
# gm.mesh.setTransfiniteCurve(6, 100)

# Uncomment for refinement
gm.mesh.field.add("Distance", 4)
# gm.mesh.field.setNumbers(4, "PointsList", [12, 13, 14, 15, 16, 11])
gm.mesh.field.setNumbers(4, "CurvesList", [6])
gm.mesh.field.setNumber(4, "Sampling", 300)
gmsh.model.mesh.field.add("MathEval", 12)
# gmsh.model.mesh.field.setString(12, "F", ".045+.65*(F4/8)^1.5")
gmsh.model.mesh.field.setString(12, "F", ".06+.45*(F4/8)^1.5")

# Use the minimum of all the fields as the background mesh field:
gm.mesh.field.add("Min", 20)
gm.mesh.field.setNumbers(20, "FieldsList", [12])

gm.mesh.field.setAsBackgroundMesh(20)
gmsh.option.setNumber("Mesh.MeshSizeMax", 60)
# gmsh.option.setNumber("Mesh.MeshSizeMin", 60)
gm.mesh.setSmoothing(2, 1, 10)

gm.mesh.generate(2)
gmsh.write('streamer_tmp.msh3')
gmsh.fltk.run()
gmsh.finalize()
exit()
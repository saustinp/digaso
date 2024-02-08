

import gmsh
import numpy as np
from gmsh import model as gm

#############################################
case_name = 'chen_geom'
#############################################

gmsh.initialize()
gm.add(case_name)

############## CREATE PLANE GEOMETRY ##############

gm.occ.importShapes('double_electrode-2.STEP')

gm.occ.synchronize()
# # print(gm.getEntities())
# gmsh.fltk.run()
# exit()

# Exraction of the bbox parameters
xmin, ymin, zmin, xmax, ymax, zmax = gm.occ.getBoundingBox(2, 1)
size_scale_factor = 100/(xmax-xmin)
# print(xmin)
# print(xmax)
# print(ymin)
# print(ymax)
# exit()

# Scale and rotate mesh -> This is to align the imported CAD geometry with the desired CSYS and scaling for the problem.
dimtags = gm.occ.getEntities()
gmsh.model.occ.dilate(dimtags, 0,0,0, size_scale_factor, size_scale_factor, size_scale_factor)

# gmsh.model.occ.rotate(dimtags, 0,0,0, 1, 0, 0, np.pi/2)
xmin, ymin, zmin, xmax, ymax, zmax = gm.occ.getBoundingBox(2, 1)    # Need to over-write the BBox values so that we can generate an accurate length scale later

# internal_tags = []
# internal_tags.append(gmsh.model.occ.addPoint(0, 98, 0, tag=10))      # 10 set arbitrarily
# internal_tags.append(gmsh.model.occ.addPoint(0, 96, 0, tag=11))      # 10 set arbitrarily
# internal_tags.append(gmsh.model.occ.addPoint(0, 94, 0, tag=12))      # 10 set arbitrarily
# internal_tags.append(gmsh.model.occ.addPoint(0, 92, 0, tag=13))      # 10 set arbitrarily
# internal_tags.append(gmsh.model.occ.addPoint(0, 90, 0, tag=14))      # 10 set arbitrarily
# internal_tags.append(gmsh.model.occ.addPoint(0, 88, 0, tag=15))      # 10 set arbitrarily
# internal_tags.append(gmsh.model.occ.addPoint(0, 86, 0, tag=16))      # 10 set arbitrarily
# internal_tags.append(gmsh.model.occ.addPoint(0, 84, 0, tag=17))      # 10 set arbitrarily
# internal_tags.append(gmsh.model.occ.addPoint(0, 82, 0, tag=18))      # 10 set arbitrarily
# internal_tags.append(gmsh.model.occ.addPoint(0, 80, 0, tag=19))      # 10 set arbitrarily
# internal_tags.append(gmsh.model.occ.addPoint(0, 78, 0, tag=20))      # 10 set arbitrarily
# # gmsh.model.mesh.embed(0, internal_tags,1, 2)       DON'T USE THIS
# for tag in internal_tags:
#     gmsh.model.occ.fragment([(0,tag)], gmsh.model.occ.getEntities())

gm.occ.synchronize()    # This is important!
# gmsh.fltk.run()
# gmsh.finalize()
# exit()

# Refinement
gm.mesh.field.add("Distance", 4)
# gm.mesh.field.setNumbers(4, "PointsList", [2,3])
gm.mesh.field.setNumbers(4, "CurvesList", [   2  ])
gm.mesh.field.setNumber(4, "Sampling", 300)
gmsh.model.mesh.field.add("MathEval", 12)
gmsh.model.mesh.field.setString(12, "F", "(.1+.5*(F4/10)^1.5)*.8")

gm.mesh.field.add("Distance", 5)
gm.mesh.field.setNumbers(5, "PointsList", [2,3])
# gm.mesh.field.setNumbers(4, "CurvesList", [   2  ])
gm.mesh.field.setNumber(4, "Sampling", 300)
gmsh.model.mesh.field.add("MathEval", 13)
gmsh.model.mesh.field.setString(13, "F", "(.03+.5*(F5/10)^1.5)*.8")

# Use the minimum of all the fields as the background mesh field:
gm.mesh.field.add("Min", 20)
gm.mesh.field.setNumbers(20, "FieldsList", [12,13])

gm.mesh.field.setAsBackgroundMesh(20)
gmsh.option.setNumber("Mesh.MeshSizeMax", 10)
# gmsh.option.setNumber("Mesh.MeshSizeMin", 60)
gm.mesh.setSmoothing(2, 1, 10)

gm.mesh.generate(2)
gmsh.write('streamer_tmp.msh3')
gmsh.fltk.run()
gmsh.finalize()
exit()

# %%

import vtk
from src.hierarchicalMesh import HierarchicalMesh
from projectionStructure import ProjectionStructure

anchor = HierarchicalMesh(None, None, None)

filename = "../meshes/mid_mesh.stl"
# the zero is just for color assignment
mesh = ProjectionStructure(filename,0)
# this way the mesh is interpreted as projection structure and a new papermesh is generated
hm = HierarchicalMesh(None,mesh,None)

# this way the mesh is assigned as papermesh with no structures inside the hm
'''
hm = HierarchicalMesh(None,None,None)
hm.papermesh = mesh.getActor().GetMapper().GetInput()
hm.mesh = mesh.getActor()
hm.setName(hm.writePapermeshStlAndOff("Temp"))
'''

anchor.add(hm)

filename = "../meshes/inner_mesh.stl"
mesh = ProjectionStructure(filename,0)
hm = HierarchicalMesh(None,mesh,None)

anchor.add(hm)

filename = "../meshes/outer_mesh.stl"
mesh = ProjectionStructure(filename,0)
hm = HierarchicalMesh(None,mesh,None)

anchor.add(hm)

anchor.toString()

print(anchor.toList())
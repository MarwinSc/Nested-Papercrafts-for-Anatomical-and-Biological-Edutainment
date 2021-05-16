import vtk
from hierarchicalMesh import HierarchicalMesh




anchor = HierarchicalMesh(None, None, "Anchor")


stlReader = vtk.vtkSTLReader()
stlReader.SetFileName("../meshes/mid_mesh.stl")
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(stlReader.GetOutputPort())
actor = vtk.vtkActor()
actor.SetMapper(mapper)

anchor.add(actor, "../meshes/mid_mesh.stl")

stlReader.SetFileName("../meshes/inner_mesh.stl")
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(stlReader.GetOutputPort())
actor = vtk.vtkActor()
actor.SetMapper(mapper)

anchor.add(actor, "../meshes/inner_mesh.stl")

stlReader.SetFileName("../meshes/outer_mesh.stl")
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(stlReader.GetOutputPort())
actor = vtk.vtkActor()
actor.SetMapper(mapper)

anchor.add(actor, "../meshes/outer_mesh.stl")

stlReader.SetFileName("../meshes/inner_mesh.stl")
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(stlReader.GetOutputPort())
actor = vtk.vtkActor()
actor.SetMapper(mapper)

anchor.add(actor, "../meshes/inner_mesh.stl")

print(anchor.toList())
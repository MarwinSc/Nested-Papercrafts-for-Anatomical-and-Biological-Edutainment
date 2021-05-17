import vtk
from PyQt5 import QtWidgets
import sys

from hierarchicalMesh import HierarchicalMesh
from viewpointselection.app import MainWindow
import viewpointselection.entropy

anchor = HierarchicalMesh(None, None, "Anchor")


stlReader = vtk.vtkSTLReader()
stlReader.SetFileName("../meshes/mid_mesh.stl")
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(stlReader.GetOutputPort())
actor = vtk.vtkActor()
actor.SetMapper(mapper)
actor.GetProperty().SetOpacity(0.5)
actor.GetProperty().SetColor(0, 1, 0)

anchor.add(actor, "../meshes/mid_mesh.stl")

stlReader.SetFileName("../meshes/inner_mesh.stl")
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(stlReader.GetOutputPort())
actor = vtk.vtkActor()
actor.SetMapper(mapper)
actor.GetProperty().SetOpacity(0.5)
actor.GetProperty().SetColor(1, 0, 0)

anchor.add(actor, "../meshes/inner_mesh.stl")

stlReader.SetFileName("../meshes/outer_mesh.stl")
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(stlReader.GetOutputPort())
actor = vtk.vtkActor()
actor.SetMapper(mapper)
actor.GetProperty().SetOpacity(0.5)
actor.GetProperty().SetColor(0, 0, 1)

anchor.add(actor, "../meshes/outer_mesh.stl")

stlReader.SetFileName("../meshes/inner_mesh.stl")
mapper = vtk.vtkPolyDataMapper()
mapper.SetInputConnection(stlReader.GetOutputPort())
actor = vtk.vtkActor()
actor.SetMapper(mapper)
actor.GetProperty().SetOpacity(0.5)
actor.GetProperty().SetColor(1, 1, 0)

anchor.add(actor, "../meshes/inner_mesh.stl")

print(anchor.toList())

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    ec = MainWindow(anchor.toList())


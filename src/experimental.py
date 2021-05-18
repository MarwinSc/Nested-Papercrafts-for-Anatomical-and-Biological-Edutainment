import vtk
from PyQt5 import QtWidgets
import sys

from hierarchicalMesh import HierarchicalMesh
from viewpointselection.app import ViewPointComputation


def build_hierarchy():
    anchor = HierarchicalMesh(None, None, "Anchor")
    stlReader = vtk.vtkSTLReader()
    stlReader.SetFileName("../meshes/mid_mesh.stl")
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(stlReader.GetOutputPort())
    actor1 = vtk.vtkActor()
    actor1.SetMapper(mapper)
    actor1.GetProperty().SetOpacity(0.2)
    actor1.GetProperty().SetColor(0.0, 1.0, 0.0)

    anchor.add(actor1, "../meshes/mid_mesh.stl")

    stlReader.SetFileName("../meshes/inner_mesh.stl")
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(stlReader.GetOutputPort())
    actor2 = vtk.vtkActor()
    actor2.SetMapper(mapper)
    actor2.GetProperty().SetOpacity(0.2)
    actor2.GetProperty().SetColor(1.0, 0.0, 0.0)

    anchor.add(actor2, "../meshes/inner_mesh.stl")

    stlReader.SetFileName("../meshes/outer_mesh.stl")
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(stlReader.GetOutputPort())
    actor3 = vtk.vtkActor()
    actor3.SetMapper(mapper)
    actor3.GetProperty().SetOpacity(0.2)
    actor3.GetProperty().SetColor(0.0, 0.0, 1.0)

    anchor.add(actor3, "../meshes/outer_mesh.stl")

    stlReader.SetFileName("../meshes/inner_mesh.stl")
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(stlReader.GetOutputPort())
    actor4 = vtk.vtkActor()
    actor4.SetMapper(mapper)
    actor4.GetProperty().SetOpacity(0.2)
    actor4.GetProperty().SetColor(1.0, 1.0, 0.0)

    anchor.add(actor4, "../meshes/inner_mesh.stl")
    return anchor


if __name__ == "__main__":
    anchor = build_hierarchy()
    print(anchor.toList())
    app = QtWidgets.QApplication(sys.argv)
    vpc = ViewPointComputation(anchor.toList())
    #vpc = ViewPointComputation()
    app.exec_()
    vpc.close()

    cut_plane = vpc.planes[0]
    plane = vtk.vtkPlane()
    plane.SetOrigin(cut_plane.GetPosition())
    plane.SetNormal(cut_plane.GetOrientation())
    anchor.cut_with_plane(plane)

    # anchor.convexify -> anchor.simplify -> anchor.recusrive_difference
    # for plane in vpc.planes():
        # anchor.cut_with_plane (?)
        # is_stable = vpc.stability([anchor.all_previously_generated_objects()], "../viewpointselection/data/urdf/plane.urdf")
        # if is_stable:
            # break
    # anchor.unfold_all_previously_generated_objects()
    print(vpc.stability(anchor.generated_by_plane_cut, "../viewpointselection/data/urdf/plane.urdf"))

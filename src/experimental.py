import vtk
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy


def read_and_render_model(path, ambient, render):
    reader = vtk.vtkOBJReader()
    reader.SetFileName(path)
    reader.Update()
    model_mapper = vtk.vtkPolyDataMapper()
    model_mapper.SetInputConnection(reader.GetOutputPort())
    model_actor = vtk.vtkActor()
    model_actor.SetMapper(model_mapper)
    model_actor.GetProperty().SetAmbient(ambient)

    prop = model_actor.GetProperty()

    model_actor.GetProperty().EdgeVisibilityOn()
    model_mapper = vtk.vtkPolyDataMapper()
    model_mapper.SetInputConnection(reader.GetOutputPort())
    render.AddActor(model_actor)


def read_and_render_IDs(path, path_mirror, render, scale=0.1):
    reader = vtk.vtkOBJReader()
    reader.SetFileName(path)
    reader.Update()

    polydata = reader.GetOutput()
    points = vtk_to_numpy(polydata.GetPoints().GetData()).astype(float)
    indices = np.unique(points, axis=0, return_index=True)[1]
    points = [points[index] for index in sorted(indices)]

    reader.SetFileName(path_mirror)
    reader.Update()

    polydata = reader.GetOutput()
    mirrorPoints = vtk_to_numpy(polydata.GetPoints().GetData()).astype(float)
    indices = np.unique(mirrorPoints, axis=0, return_index=True)[1]
    mirrorPoints = [mirrorPoints[index] for index in sorted(indices)]

    id = 1
    for i in range(0, len(points), 4):
        label = vtk.vtkVectorText()
        label.SetText(str(id))
        id += 1
        pos = (points[i] + points[i + 1] + points[i + 2] + points[i + 3]) / 4

        mirrorPos = (mirrorPoints[i] + mirrorPoints[i + 1] + mirrorPoints[i + 2] + mirrorPoints[i + 3]) / 4

        lblMapper = vtk.vtkPolyDataMapper()
        lblMapper.SetInputConnection(label.GetOutputPort())

        # Set up an actor for the node label
        lblActor = vtk.vtkFollower()
        lblActor.SetMapper(lblMapper)
        lblActor.SetScale(scale, scale, scale)
        lblActor.SetPosition(pos[0] - scale / 2, pos[1] - scale / 2, 0.0)
        lblActor.GetProperty().SetColor(255, 0, 0)
        render.AddActor(lblActor)

        lblMapper = vtk.vtkPolyDataMapper()
        lblMapper.SetInputConnection(label.GetOutputPort())
        # Set up an actor for the node label
        lblActor = vtk.vtkFollower()
        lblActor.SetMapper(lblMapper)
        lblActor.SetScale(scale, scale, scale)
        lblActor.SetPosition(mirrorPos[0] - scale / 2, mirrorPos[1] - scale / 2, 0.0)
        lblActor.GetProperty().SetColor(0, 255, 0)
        render.AddActor(lblActor)

def main():
    renderer = vtk.vtkRenderer()

    read_and_render_model(r"C:\repositories\AnatomicalEdutainer\meshes\model.obj", 0.8, renderer)
    read_and_render_model(r"C:\repositories\AnatomicalEdutainer\meshes\gluetabs.obj", 0.3, renderer)
    read_and_render_model(r"C:\repositories\AnatomicalEdutainer\meshes\mirrorgt.obj", 0.2, renderer)

    read_and_render_IDs(r"C:\repositories\AnatomicalEdutainer\meshes\gluetabs.obj", r"C:\repositories\AnatomicalEdutainer\meshes\mirrorgt.obj", renderer)

    renderer.SetBackground(vtk.vtkNamedColors().GetColor3d('White'))
    window = vtk.vtkRenderWindow()
    window.AddRenderer(renderer)
    window.SetSize(720, 480)

    interactiveRenderer = vtk.vtkRenderWindowInteractor()
    interactiveRenderer.SetRenderWindow(window)
    style = vtk.vtkInteractorStyleTrackballCamera()
    interactiveRenderer.SetInteractorStyle(style)
    interactiveRenderer.Initialize()
    interactiveRenderer.Start()


if __name__ == "__main__":
    main()

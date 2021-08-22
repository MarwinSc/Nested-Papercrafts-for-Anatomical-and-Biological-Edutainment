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


def read_and_render_IDs(path, ambient, render):
    reader = vtk.vtkOBJReader()
    reader.SetFileName(path)
    reader.Update()

    polydata = reader.GetOutput()
    faces = vtk_to_numpy(polydata.GetPolys().GetData())
    points = vtk_to_numpy(polydata.GetPoints().GetData()).astype(float)

    id = 1
    for i in range(0, len(faces) - 2, 5):
        txt = vtk.vtkTextActor()
        txt.SetInput(str(id))
        id += 1
        txtprop = txt.GetTextProperty()
        txtprop.SetFontFamilyToArial()
        txtprop.SetFontSize(18)
        txt.GetProperty().SetColor(ambient, 0, 0)
        pos = (points[faces[i]] + points[faces[i + 1]] + points[faces[i + 2]] + points[faces[i + 3]] + points[
            faces[i + 4]] + points[faces[i + 5]]) / 6

        print(pos)
        txt.SetPosition(pos[0], pos[1])
        render.AddActor(txt)


def main():
    renderer = vtk.vtkRenderer()

    read_and_render_model(r"C:\repositories\AnatomicalEdutainer\meshes\model.obj", 0.8, renderer)
    read_and_render_model(r"C:\repositories\AnatomicalEdutainer\meshes\gluetabs.obj", 0.3, renderer)
    read_and_render_model(r"C:\repositories\AnatomicalEdutainer\meshes\mirrorgt.obj", 0.2, renderer)

    read_and_render_IDs(r"C:\repositories\AnatomicalEdutainer\meshes\gluetabs.obj", 1.0, renderer)

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

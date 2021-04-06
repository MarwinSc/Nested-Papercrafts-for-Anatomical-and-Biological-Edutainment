import vtkmodules.all as vtk
import numpy as np
from vtkmodules.numpy_interface.dataset_adapter import numpy_support

def getbufferRenIntWin(camera = vtk.vtkCamera(),width=2000,height=2000):
    ren = vtk.vtkRenderer()
    ren.SetBackground(255.0, 255.0, 255.0)
    ren.SetActiveCamera(camera)
    renWin = vtk.vtkRenderWindow()
    renWin.SetSize(width,height)
    renWin.AddRenderer(ren)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    wti = vtk.vtkWindowToImageFilter()
    wti.SetInput(renWin)
    #wti.SetScale(1000,1000)
    wti.SetInputBufferTypeToRGB()
    wti.ReadFrontBufferOff()
    return ren,iren,renWin,wti

def printSinglePoint(self, x, y, ren, color):
    newPoints = vtk.vtkPoints()
    vertices = vtk.vtkCellArray()

    p = [x, 0.0, y]

    id = newPoints.InsertNextPoint(p)
    vertices.InsertNextCell(1)
    vertices.InsertCellPoint(id)

    # Create a polydata object
    point = vtk.vtkPolyData()

    point.SetPoints(newPoints)
    point.SetVerts(vertices)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(point)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(color)
    actor.GetProperty().SetPointSize(2)
    ren.AddActor(actor)

def NpToVtk(img,dx,dy,dz):
    resultImg = vtk.vtkImageData()

    vtkResult = numpy_support.numpy_to_vtk(img.reshape(dy * dx, dz))

    resultImg.SetSpacing(1., 1., 1.)
    resultImg.SetOrigin(0., 0., 0.)
    resultImg.SetDimensions(dx, dy, 1)
    resultImg.AllocateScalars(numpy_support.get_vtk_array_type(img.dtype), dz)
    resultImg.GetPointData().SetScalars(vtkResult)
    return resultImg

def writeImage(img, path):
    castFilter = vtk.vtkImageCast()
    castFilter.SetInputData(img)
    castFilter.SetOutputScalarTypeToUnsignedChar()
    castFilter.Update()

    writer = vtk.vtkPNGWriter()
    writer.SetFileName(path)
    writer.SetInputConnection(castFilter.GetOutputPort())
    writer.Write()
    return castFilter.GetOutput()
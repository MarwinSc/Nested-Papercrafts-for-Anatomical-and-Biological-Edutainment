import vtkmodules.all as vtk

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
    wti.SetInputBufferTypeToRGB()
    wti.ReadFrontBufferOff()
    return ren,iren,renWin,wti

def printSinglePoint(self,p):
    newPoints = vtk.vtkPoints()
    vertices = vtk.vtkCellArray()

    id = newPoints.InsertNextPoint(p)
    vertices.InsertNextCell(1)
    vertices.InsertCellPoint(id)

    # Create a polydata object
    point = vtk.vtkPolyData()

    # Set the points and vertices we created as the geometry and topology of the polydata
    point.SetPoints(newPoints)
    point.SetVerts(vertices)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(point)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor([255.0, 0.0, 0.0])
    actor.GetProperty().SetPointSize(10)
    #self.ren.AddActor(actor)
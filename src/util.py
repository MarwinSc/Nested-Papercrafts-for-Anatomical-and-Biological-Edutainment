import vtkmodules.all as vtk
import numpy as np
import os
import trimesh
import meshio
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

def stlToOff(meshpath):
    dirname = os.path.dirname(__file__)
    filename = os.path.join(dirname, meshpath)
    mesh = trimesh.load(filename)
    filename = os.path.join(dirname, "../out/3D/mesh.off")
    trimesh.exchange.export.export_mesh(mesh, filename, file_type="off")

def meshioIO(inPath,outPath):
    mesh = meshio.read(
        inPath,  # string, os.PathLike, or a buffer/open file
        #file_format="stl",  # optional if filename is a path; inferred from extension
    )
    meshio.write(
        outPath,  # str, os.PathLike, or buffer/ open file
        mesh,
        # file_format="vtk",  # optional if first argument is a path; inferred from extension
    )

def writeStl(mesh,name):
    dirname = os.path.dirname(__file__)
    filename = os.path.join(dirname, "../out/3D/"+name+".stl")
    stlWriter = vtk.vtkSTLWriter()
    stlWriter.SetFileName(filename)
    stlWriter.SetInputData(mesh)
    stlWriter.Write()

def writeObj(mesh,name):
    dirname = os.path.dirname(__file__)
    filename = os.path.join(dirname, "../out/3D/"+name+".obj")
    objWriter = vtk.vtkOBJWriter()
    objWriter.SetFileName(filename)
    objWriter.SetInputData(mesh)
    objWriter.Write()

def shrinkWrap(mesh,shrink_mesh, edgeSmoothing = True):
    smoother = vtk.vtkSmoothPolyDataFilter()
    smoother.SetInputData(0, mesh)
    smoother.SetInputData(1, shrink_mesh)

    #smoother.SetNumberOfIterations(2)
    #smoother.SetRelaxationFactor(0.1)
    if edgeSmoothing:
        smoother.FeatureEdgeSmoothingOn()
    smoother.Update()

    clean = vtk.vtkCleanPolyData()
    clean.SetInputData(smoother.GetOutput())
    clean.Update()

    return clean.GetOutput()

def offsetMesh(mesh, factor = 5.0):

    clean = vtk.vtkCleanPolyData()
    normals = vtk.vtkPolyDataNormals()
    clean.SetInputData(mesh)
    normals.SetInputConnection(clean.GetOutputPort())
    normals.SplittingOff()

    offsetted = vtk.vtkWarpVector()
    offsetted.SetInputConnection(normals.GetOutputPort())
    offsetted.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS,vtk.vtkDataSetAttributes.NORMALS)
    offsetted.SetScaleFactor(factor)
    offsetted.Update()

    clean = vtk.vtkCleanPolyData()
    clean.SetInputData(offsetted.GetOutput())
    clean.Update()

    return clean.GetOutput()

def subdivideMesh(mesh, iterations = 1):
    subdivider = vtk.vtkLinearSubdivisionFilter()
    subdivider.SetNumberOfSubdivisions(iterations)
    subdivider.SetInputData(mesh)
    subdivider.Update()

    clean = vtk.vtkCleanPolyData()
    clean.SetInputData(subdivider.GetOutput())
    clean.Update()

    return clean.GetOutput()

def smoothMesh(mesh,secondMesh = None,iterations = 15,relaxation = 0.1):
    smooth = vtk.vtkSmoothPolyDataFilter()
    if secondMesh:
        smooth.SetInputData(0,mesh)
        smooth.SetInputData(1,secondMesh)
    else:
        smooth.SetInputData(mesh)
    smooth.SetNumberOfIterations(iterations)
    smooth.SetRelaxationFactor(relaxation)
    smooth.FeatureEdgeSmoothingOff()
    smooth.BoundarySmoothingOn()
    smooth.Update()

    clean = vtk.vtkCleanPolyData()
    clean.SetInputData(smooth.GetOutput())
    clean.Update()

    return clean.GetOutput()

def cutAwayBlackAreaOfImage(img):
    mask = np.where(img > [0.0, 0.0, 0.0])

    width = mask[0].max() - mask[0].min()
    height = mask[1].max() - mask[1].min()
    if width > height:
        img = img[mask[0].min():mask[0].max(), mask[1].min(): mask[1].min() + width, :]
    else:
        img = img[mask[0].min():mask[0].min() + height, mask[1].min(): mask[1].max(), :]
    return img
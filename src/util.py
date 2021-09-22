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

def VtkToNp(img):
    width = img.GetExtent()[1] + 1
    height = img.GetExtent()[3] + 1
    img = numpy_support.vtk_to_numpy(img.GetPointData().GetScalars())[:, 0:3]
    img.astype(float)
    img = np.reshape(np.ravel(img), (width, height, 3))
    return img

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

def readStl(name):
    stlReader = vtk.vtkSTLReader()
    stlReader.SetFileName(name)
    stlReader.Update()
    return stlReader.GetOutput()

def writeStl(mesh,name):
    dirname = os.path.dirname(__file__)
    filename = os.path.join(dirname, "../out/3D/"+name+".stl")
    stlWriter = vtk.vtkSTLWriter()
    stlWriter.SetFileName(filename)
    mesh = cleanMesh(mesh)
    stlWriter.SetInputData(mesh)
    stlWriter.Write()

def writeObj(mesh,name):
    dirname = os.path.dirname(__file__)
    filename = os.path.join(dirname, "../out/3D/"+name+".obj")
    objWriter = vtk.vtkOBJWriter()
    objWriter.SetFileName(filename)
    objWriter.SetInputData(mesh)
    objWriter.Write()

def readObj(path):
    importer = vtk.vtkOBJReader()
    importer.SetFileName(path)
    importer.Update()
    return importer.GetOutput()

def shrinkWrap(mesh,shrink_mesh, edgeSmoothing = True):
    smoother = vtk.vtkSmoothPolyDataFilter()
    smoother.SetInputData(0, mesh)
    smoother.SetInputData(1, shrink_mesh)

    #smoother.SetNumberOfIterations(2)
    #smoother.SetRelaxationFactor(0.1)
    if edgeSmoothing:
        smoother.FeatureEdgeSmoothingOn()
    smoother.Update()

    #clean = vtk.vtkCleanPolyData()
    #clean.SetInputData(smoother.GetOutput())
    #clean.Update()

    return smoother.GetOutput()

def offsetMesh(mesh, factor = 5.0):

    clean = vtk.vtkCleanPolyData()
    normals = vtk.vtkPolyDataNormals()
    clean.SetInputData(mesh)
    normals.SetInputConnection(clean.GetOutputPort())
    normals.SetInputData(mesh)
    normals.SplittingOff()

    offsetted = vtk.vtkWarpVector()
    offsetted.SetInputConnection(normals.GetOutputPort())
    offsetted.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS,vtk.vtkDataSetAttributes.NORMALS)
    offsetted.SetScaleFactor(factor)
    offsetted.Update()

    #clean = vtk.vtkCleanPolyData()
    #clean.SetInputData(offsetted.GetOutput())
    #clean.Update()

    return offsetted.GetOutput()

def subdivideMesh(mesh, iterations = 1):
    subdivider = vtk.vtkLinearSubdivisionFilter()
    subdivider.SetNumberOfSubdivisions(iterations)
    subdivider.SetInputData(mesh)
    subdivider.Update()

    #clean = vtk.vtkCleanPolyData()
    #clean.SetInputData(subdivider.GetOutput())
    #clean.Update()

    return subdivider.GetOutput()

def smoothMesh(mesh,secondMesh = None,iterations = 15,relaxation = 0.1):
    smooth = vtk.vtkSmoothPolyDataFilter()
    if secondMesh:
        smooth.SetInputData(0,mesh)
        smooth.SetInputData(1,secondMesh)
    else:
        smooth.SetInputData(mesh)
    smooth.SetNumberOfIterations(iterations)
    smooth.SetRelaxationFactor(relaxation)
    smooth.FeatureEdgeSmoothingOn()
    smooth.BoundarySmoothingOff()
    smooth.Update()

    #clean = vtk.vtkCleanPolyData()
    #clean.SetInputData(smooth.GetOutput())
    #clean.Update()

    return smooth.GetOutput()

def cleanMesh(mesh,tolerance = 0.0):
    clean = vtk.vtkCleanPolyData()
    clean.SetInputData(mesh)
    clean.SetTolerance(tolerance)
    clean.Update()
    return clean.GetOutput()


def cleanedMeshToMeshWithDoubleVertices(mesh):

    newGeometry = vtk.vtkPolyData()
    newPoints = vtk.vtkPoints()
    newCells = vtk.vtkCellArray()

    print(mesh.GetNumberOfCells())

    for i in range(mesh.GetNumberOfCells()):

        points = mesh.GetCell(i).GetPoints()

        triCell = vtk.vtkTriangle()

        pointId = newPoints.InsertNextPoint(points.GetPoint(0))
        triCell.GetPointIds().SetId(0, pointId)

        pointId = newPoints.InsertNextPoint(points.GetPoint(1))
        triCell.GetPointIds().SetId(1, pointId)

        pointId = newPoints.InsertNextPoint(points.GetPoint(2))
        triCell.GetPointIds().SetId(2, pointId)

        newCells.InsertNextCell(triCell)

    newGeometry.SetPoints(newPoints)
    newGeometry.SetPolys(newCells)

    return newGeometry

def projectMeshToBoundsAlongCubeNormals(mesh):
    bounds = mesh.GetBounds()
    width = bounds[1] - bounds[0]
    depth = bounds[3] - bounds[2]
    height = bounds[5] - bounds[4]

    newGeometry = vtk.vtkPolyData()
    newPoints = vtk.vtkPoints()
    newCells = vtk.vtkCellArray()

    hull = vtk.vtkHull()
    hull.SetInputData(mesh)
    hull.AddCubeFacePlanes()
    hull.Update()

    centersFilter = vtk.vtkCellCenters()
    centersFilter.SetInputData(mesh)
    centersFilter.VertexCellsOn()
    centersFilter.Update()

    cellLocator = vtk.vtkCellLocator()
    cellLocator.SetDataSet(hull.GetOutput())
    cellLocator.BuildLocator()

    for i in range(centersFilter.GetOutput().GetNumberOfPoints()):
        p = [0.0, 0.0, 0.0]
        centersFilter.GetOutput().GetPoint(i, p)
        triCell = vtk.vtkTriangle()

        '''

        #for each point
        for j in range(3):

            absCoords = [abs(coord) for coord in points.GetPoint(j)]
            indexOfHighest = absCoords.index(max(absCoords))

            p1 = list(points.GetPoint(j))
            p2 = [p1[0] + (normals.GetTuple(i)[0] * 100.0), p1[1] + (normals.GetTuple(i)[1] * 100.0), p1[2] + (normals.GetTuple(i)[2] * 100.0)]

            tolerance = 0.1

            t = vtk.mutable(0)
            x = [0.0, 0.0, 0.0]
            pcoords = [0.0, 0.0, 0.0]
            subId = vtk.mutable(0)
            cellLocator.IntersectWithLine(p1, p2, tolerance, t, x, pcoords, subId)

            pointId = newPoints.InsertNextPoint(tuple(x))
            triCell.GetPointIds().SetId(j, pointId)
        '''
        newCells.InsertNextCell(triCell)

        newGeometry.SetPoints(newPoints)
        newGeometry.SetPolys(newCells)

def projectMeshToBounds(mesh):

    bounds = mesh.GetBounds()

    width = bounds[1]-bounds[0]
    depth = bounds[3]-bounds[2]
    height = bounds[5]-bounds[4]
    threshold = [width/8,depth/8,height/8]
    newGeometry = vtk.vtkPolyData()
    newPoints = vtk.vtkPoints()
    newCells = vtk.vtkCellArray()

    hull = vtk.vtkHull()
    hull.SetInputData(mesh)
    hull.AddCubeFacePlanes()
    hull.Update()

    normalsFilter = vtk.vtkPolyDataNormals()
    normalsFilter.SetInputData(mesh)
    normalsFilter.ComputePointNormalsOn()
    normalsFilter.ComputeCellNormalsOff()
    normalsFilter.SplittingOff()
    normalsFilter.Update()
    array = normalsFilter.GetOutput()
    normals = array.GetCellData().GetArray("Normals")

    cellLocator = vtk.vtkCellLocator()
    cellLocator.SetDataSet(hull.GetOutput())
    cellLocator.BuildLocator()

    for i in range(mesh.GetNumberOfCells()):
        points = mesh.GetCell(i).GetPoints()
        triCell = vtk.vtkTriangle()

        #for each point
        for j in range(3):

            absCoords = [abs(coord) for coord in points.GetPoint(j)]
            indexOfHighest = absCoords.index(max(absCoords))

            p1 = list(points.GetPoint(j))
            p2 = [p1[0] + (normals.GetTuple(i)[0] * 100.0), p1[1] + (normals.GetTuple(i)[1] * 100.0), p1[2] + (normals.GetTuple(i)[2] * 100.0)]

            tolerance = 0.1

            t = vtk.mutable(0)
            x = [0.0, 0.0, 0.0]
            pcoords = [0.0, 0.0, 0.0]
            subId = vtk.mutable(0)
            cellLocator.IntersectWithLine(p1, p2, tolerance, t, x, pcoords, subId)

            '''
            for y in range(3):
                if abs(points.GetPoint(j)[y]) > threshold[y]:
                    if np.sign(points.GetPoint(j)[y]) == -1:
                        point[y] = bounds[(y + 1) * 2 - 2]
                    elif np.sign(points.GetPoint(j)[y]) == 1:
                        point[y] = bounds[(y + 1) * 2 - 1]
            '''

            pointId = newPoints.InsertNextPoint(tuple(x))
            triCell.GetPointIds().SetId(j, pointId)

        newCells.InsertNextCell(triCell)

    newGeometry.SetPoints(newPoints)
    newGeometry.SetPolys(newCells)

    newGeometry = smoothMesh(newGeometry, hull.GetOutput())

    

    return newGeometry


def cutAwayBlackAreaOfImage(img):
    mask = np.where(img > [0.0, 0.0, 0.0])

    width = mask[0].max() - mask[0].min()
    height = mask[1].max() - mask[1].min()
    if width > height:
        img = img[mask[0].min():mask[0].max(), mask[1].min(): mask[1].min() + width, :]
    else:
        img = img[mask[0].min():mask[0].min() + height, mask[1].min(): mask[1].max(), :]
    return img

def appendMeshes(meshes):
    append = vtk.vtkAppendPolyData()
    for a in meshes:
        append.AddInputData(a)
    append.Update()
    clean = vtk.vtkCleanPolyData()
    clean.SetInputData(append.GetOutput())
    clean.Update()
    return clean.GetOutput()
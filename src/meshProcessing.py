import vtkmodules.all as vtk
from vtkmodules.numpy_interface.dataset_adapter import numpy_support
import numpy as np
import os
import util
import meshInteraction
import trimesh
import meshio
from mu3d.mu3dpy.mu3d import Graph

#Class responsible for mesh related processing steps and I/O.
class MeshProcessing():

    dirname = os.path.dirname(__file__)

    dedicatedPaperMeshes = []
    tempPaperActor = vtk.vtkActor()

    meshInteractor = meshInteraction.MeshInteraction(dedicatedPaperMeshes)

    # create initial Mesh to subdivide and wrap
    def createPapermesh(self,actorList, graph, boolReturn = False):

        polysAppended = self.appendAllMeshes(actorList)

        hull = vtk.vtkHull()
        hull.SetInputData(polysAppended)
        hull.AddCubeFacePlanes()
        # hull.AddRecursiveSpherePlanes(1)
        hull.Update()

        triangleFilter = vtk.vtkTriangleFilter()
        triangleFilter.SetInputConnection(hull.GetOutputPort())
        triangleFilter.Update()

        #paperMesh = triangleFilter.GetOutput()

        # for i in range(2):
        #temp = self.__createPaperMeshHelper(paperMesh, hull, False, actorList)

        subdivider = vtk.vtkLinearSubdivisionFilter()
        subdivider.SetNumberOfSubdivisions(1)
        subdivider.SetInputConnection(triangleFilter.GetOutputPort())

        smoother = vtk.vtkSmoothPolyDataFilter()
        smoother.SetInputConnection(0, subdivider.GetOutputPort())
        smoother.SetInputData(1, polysAppended)
        # smoother.SetNumberOfIterations(2)
        # smoother.SetRelaxationFactor(0.1)
        smoother.FeatureEdgeSmoothingOn()

        clean = vtk.vtkCleanPolyData()
        normals = vtk.vtkPolyDataNormals()
        clean.SetInputConnection(smoother.GetOutputPort())
        normals.SetInputConnection(clean.GetOutputPort())
        normals.SplittingOff()
        offsetted = vtk.vtkWarpVector()
        offsetted.SetInputConnection(normals.GetOutputPort())
        offsetted.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, vtk.vtkDataSetAttributes.NORMALS)
        offsetted.SetScaleFactor(5.0)
        offsetted.Update()

        paperMesh = offsetted.GetOutput()

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(paperMesh)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        #ren, iren, renWin, wti = util.getbufferRenIntWin(vtk.vtkCamera(),2000,2000)
        #ren.AddActor(actor)
        #renWin.Render()

        # ----write for blender
        #objWriter = vtk.vtkOBJExporter()
        # objWriter.SetRenderWindow(ren.GetRenderWindow())
        #objWriter.SetInput(renWin)

        filename = os.path.join(self.dirname, "../out/3D/papermesh.stl")
        #objWriter.SetFilePrefix(filename)
        #objWriter.Write()

        stlWriter = vtk.vtkSTLWriter()
        stlWriter.SetFileName(filename)
        stlWriter.SetInputData(actor.GetMapper().GetInput())
        stlWriter.Write()

        mesh = meshio.read(
            filename,  # string, os.PathLike, or a buffer/open file
            file_format="stl",  # optional if filename is a path; inferred from extension
        )

        filename = os.path.join(self.dirname, "../out/3D/papermesh.off")

        meshio.write(
            filename,  # str, os.PathLike, or buffer/ open file
            mesh,
            # file_format="vtk",  # optional if first argument is a path; inferred from extension
        )

        graph.load(filename)
        if not graph.unfold(10000, 0):
            raise Exception("unfolding failed")
        else:
            filename = os.path.join(self.dirname, "../out/3D/unfolded/model.obj")
            gluetabs_filename = os.path.join(self.dirname, "../out/3D/unfolded/gluetabs.obj")

            graph.save(filename, gluetabs_filename)

        if not boolReturn:
            self.tempPaperActor = actor

        return actor

    def appendAllMeshes(self, actorList):

        append = vtk.vtkAppendPolyData()
        clean = vtk.vtkCleanPolyData()

        for a in actorList:
            append.AddInputData(a.GetMapper().GetInput())

        clean.SetInputConnection(append.GetOutputPort())
        clean.Update()

        return clean.GetOutput()

    def __createPaperMeshHelper(self, mesh, hull, isDedicatedMesh, actorList):

        subdivider = vtk.vtkLinearSubdivisionFilter()
        smoother = vtk.vtkSmoothPolyDataFilter()
        clean = vtk.vtkCleanPolyData()
        normals = vtk.vtkPolyDataNormals()
        offsetted = vtk.vtkWarpVector()

        triangleFilter = vtk.vtkTriangleFilter()
        triangleFilter.SetInputConnection(hull.GetOutputPort())

        ##deprecated
        if isDedicatedMesh:
            subdivider.SetNumberOfSubdivisions(1)
            subdivider.SetInputConnection(triangleFilter.GetOutputPort())
            smoother.SetInputConnection(0, subdivider.GetOutputPort())
            smoother.SetInputData(1, mesh)
        else:
            subdivider.SetNumberOfSubdivisions(1)
            subdivider.SetInputConnection(triangleFilter.GetOutputPort())

            smoother.SetInputConnection(0, subdivider.GetOutputPort())
            smoother.SetInputData(1, self.appendAllMeshes(actorList))

        clean.SetInputConnection(smoother.GetOutputPort())
        normals.SetInputConnection(clean.GetOutputPort())
        normals.SplittingOff()

        ##deprecated
        if (isDedicatedMesh):
            offsetted.SetInputConnection(normals.GetOutputPort())
            offsetted.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS,vtk.vtkDataSetAttributes.NORMALS)
            offsetted.SetScaleFactor(5.0)
            return offsetted
        else:
            offsetted.SetInputConnection(normals.GetOutputPort())
            offsetted.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS,vtk.vtkDataSetAttributes.NORMALS)
            offsetted.SetScaleFactor(15.0)
            return offsetted

    # for each loaded structure shrinkwrap the global paper mesh onto them
    def createDedicatedMeshes(self,inflateStrucList,actorList):

        for a in range(len(actorList)):

            if inflateStrucList[a]:

                hull = vtk.vtkHull()
                hull.SetInputData(actorList[a].GetMapper().GetInput())
                hull.AddCubeFacePlanes()
                hull.Update()

                smoother = vtk.vtkSmoothPolyDataFilter()

                smoother.SetInputData(0, self.tempPaperActor.GetMapper().GetInput())
                #!!!!
                meshToShrinkOnto = self.__createPaperMeshHelper(actorList[a].GetMapper().GetInput(),hull,True,actorList)
                smoother.SetInputConnection(1, meshToShrinkOnto.GetOutputPort())

                clean = vtk.vtkCleanPolyData()
                clean.SetInputConnection(smoother.GetOutputPort())

                offsetted = vtk.vtkWarpVector()
                offsetted.SetInputConnection(clean.GetOutputPort())
                offsetted.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, vtk.vtkDataSetAttributes.NORMALS)
                offsetted.SetScaleFactor(5.0)
                offsetted.Update()

                poly = offsetted.GetOutput()

            else:

                poly = self.tempPaperActor.GetMapper().GetInput()

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputData(poly)

            actor = vtk.vtkActor()
            actor.SetMapper(mapper)

            self.dedicatedPaperMeshes.append(actor)


    def importUnfoldedMesh(self, boolReturn = False):

        #ren, iren, renWin, wti = util.getbufferRenIntWin()

        importer = vtk.vtkOBJReader()
        #importer = vtk.vtkOBJImporter()

        #filename = os.path.join(self.dirname, "../out/3D/unfolded/model.obj")

        #mesh = meshio.read(filename, "obj")
        #meshio.write(filename, mesh)

        print("import geht nicht")

        filename = os.path.join(self.dirname, "../out/3D/unfolded/model.obj")
        importer.SetFileName(filename)
        #filename = os.path.join(self.dirname, "../out/3D/unfolded/mesh.mtl")
        #importer.SetFileNameMTL(filename)
        #filename = os.path.join(self.dirname, ".../out/3D/unfolded/")
        #importer.SetTexturePath(filename)
        #importer.SetRenderWindow(renWin)

        print("was will er eig")

        importer.Update()
        mesh = importer.GetOutput()

        print("import geht")

        textureCoordinates = mesh.GetPointData().GetTCoords()
        absmax = 0.0
        normmax = 0.0
        newTCoords = vtk.vtkFloatArray()
        newTCoords.SetNumberOfComponents(2)

        print("get tcoords geht")

        for i in range(int(textureCoordinates.GetNumberOfTuples()/3)):
            uvs = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
            textureCoordinates.GetTuple(mesh.GetCell(i).GetPointId(0), uvs[0])
            textureCoordinates.GetTuple(mesh.GetCell(i).GetPointId(1), uvs[1])
            textureCoordinates.GetTuple(mesh.GetCell(i).GetPointId(2), uvs[2])
            abstemp=np.amax(abs(np.array(uvs)))
            normtemp=np.amax((np.array(uvs)))
            if abstemp > absmax: absmax = abstemp
            if normtemp > normmax: normmax = normtemp

        print("1 for loop geht")

        for i in range(int(textureCoordinates.GetNumberOfTuples()/3)):
            u = (textureCoordinates.GetTuple2((mesh.GetCell(i).GetPointId(0)))[0]+absmax)/(absmax+normmax)
            v = (textureCoordinates.GetTuple2((mesh.GetCell(i).GetPointId(0)))[1]+absmax)/(absmax+normmax)
            newTCoords.InsertNextTuple2(u,v)

            u = (textureCoordinates.GetTuple2((mesh.GetCell(i).GetPointId(1)))[0]+absmax)/(absmax+normmax)
            v = (textureCoordinates.GetTuple2((mesh.GetCell(i).GetPointId(1)))[1]+absmax)/(absmax+normmax)
            newTCoords.InsertNextTuple2(u,v)

            u = (textureCoordinates.GetTuple2((mesh.GetCell(i).GetPointId(2)))[0]+absmax)/(absmax+normmax)
            v = (textureCoordinates.GetTuple2((mesh.GetCell(i).GetPointId(2)))[1]+absmax)/(absmax+normmax)
            newTCoords.InsertNextTuple2(u,v)

        mesh.GetPointData().SetTCoords(newTCoords)

        print("TCoords gehen")

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(mesh)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        #actors = ren.GetActors()
        #actor = actors.GetLastActor()

        #transform = vtk.vtkTransform()
        #transform.RotateX(90.0)
        #rotate = vtk.vtkTransformFilter()
        #rotate.SetInputData(actor.GetMapper().GetInput())
        #rotate.SetTransform(transform)
        #rotate.Update()

        #actor.GetMapper().SetInputConnection(rotate.GetOutputPort())
        #actor.Modified()

        actor.GetProperty().SetColor([0.5,0.5,0.5])
        actor.GetProperty().BackfaceCullingOn()

        if not boolReturn:
            self.tempPaperActor = actor
        return actor

    #Enormous rendering method that just produces a long texture image, of concatenated renderings of the triangles from the papermesh.
    def projectPerTriangle(self,inflateStruc,actorList):
        meshNr = 0

        for a in self.dedicatedPaperMeshes:
            paper = a.GetMapper().GetInput()
            textureCoordinates = paper.GetPointData().GetTCoords()

            centersFilter = vtk.vtkCellCenters()
            centersFilter.SetInputData(paper)
            centersFilter.VertexCellsOn()
            centersFilter.Update()

            normals = vtk.vtkPolyDataNormals()
            normals.SetInputData(paper)

            normals.ComputePointNormalsOff()
            normals.ComputeCellNormalsOn()
            normals.SplittingOff()
            normals.FlipNormalsOn()
            normals.Update()

            array = normals.GetOutput()
            normalDataDouble = array.GetCellData().GetArray("Normals")

            # -----------------

            camera = vtk.vtkCamera()
            #        camera.SetViewUp(0, 1, 0)
            camera.ParallelProjectionOn()
            camera.Zoom(0.01)
            camera.SetClippingRange(0.0001, 500.01)
            if inflateStruc:
                if not inflateStruc[meshNr]:
                    camera.SetClippingRange(0.0001, 60.01)

            depthPeeling = True
            occlusion = 0.1
            numberOfPeels = 10

            buffer = vtk.vtkRenderer()
            buffer.SetBackground(255.0, 255.0, 255.0)
            buffer.SetActiveCamera(camera)
            buffer.SetLayer(0)

            bufferPaper = vtk.vtkRenderer()
            buffer.SetBackground(255.0, 255.0, 255.0)
            buffer.SetActiveCamera(camera)
            bufferPaper.SetLayer(1)

            bufferPoints = vtk.vtkRenderer()
            bufferPoints.SetBackground(255.0, 255.0, 255.0)
            bufferPoints.SetActiveCamera(camera)
            bufferPoints.SetLayer(2)

            bufferWin = vtk.vtkRenderWindow()
            bufferWin.SetNumberOfLayers(3)
            w = 500
            h = 500
            bufferWin.SetSize(w, h)
            bufferWin.AddRenderer(buffer)
            bufferWin.AddRenderer(bufferPaper)
            bufferWin.AddRenderer(bufferPoints)

            bufferIren = vtk.vtkRenderWindowInteractor()
            bufferIren.SetRenderWindow(bufferWin)
            bufferWin.SetOffScreenRendering(True)

            if depthPeeling:
                buffer.SetUseDepthPeeling(True)
                buffer.SetOcclusionRatio(occlusion)
                buffer.SetMaximumNumberOfPeels(numberOfPeels)
            else:
                buffer.SetUseDepthPeeling(False)

            # ?
            #bufferPaper.AddActor(a)

            paperMeshColor = [0.0, 0.0, 0.0, 255.0]
            transparent = [0.0, 0.0, 0.0, 0.0]
            cellData = vtk.vtkUnsignedCharArray()
            cellData.SetNumberOfComponents(4)
            cellData.SetNumberOfTuples(centersFilter.GetOutput().GetNumberOfPoints())
            for i in range(centersFilter.GetOutput().GetNumberOfPoints()):
                cellData.InsertTuple(i, transparent)

            paper.GetCellData().SetScalars(cellData)
            paper.GetCellData().Modified()
            paper.Modified()

            uvArray = vtk.vtkDoubleArray()
            uvArray.SetNumberOfComponents(2)
            newGeometry = vtk.vtkPolyData()
            newPoints = vtk.vtkPoints()
            newCells = vtk.vtkCellArray()

            # helper for placement
            lengths = []
            img = np.array([[],[],[]])

            print("zeile 370")

            for i in range(centersFilter.GetOutput().GetNumberOfPoints()):
                p = [0.0, 0.0, 0.0]
                centersFilter.GetOutput().GetPoint(i, p)

                print("zeile 376")

                p2 = normalDataDouble.GetTuple3(i)

                position = [p[0] + (p2[0] * -5), p[1] + (p2[1] * -5), p[2] + (p2[2] * -5)]

                # p2 = [p[0] + (p2[0]*10), p[1] + (p2[1]*10), p[2] + (p2[2]*10)]
                print("zeile 383")

                camera.SetPosition(position)
                camera.SetFocalPoint(p)
                #            camera.OrthogonalizeViewUp()

                points = paper.GetCell(i).GetPoints()
                print("zeile 390")

                uvs = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
                textureCoordinates.GetTuple(paper.GetCell(i).GetPointId(0), uvs[0])
                textureCoordinates.GetTuple(paper.GetCell(i).GetPointId(1), uvs[1])
                textureCoordinates.GetTuple(paper.GetCell(i).GetPointId(2), uvs[2])
                print("zeile 394")

                uvs = np.array(uvs)
                uvs = list(zip(uvs[:, 1], uvs[:, 0]))

                # set one cell opaque
                #cellData.SetTuple(i, transparent)
                #paper.GetCellData().SetScalars(cellData)
                #paper.GetCellData().Modified()
                #paper.Modified()

                mapper = vtk.vtkPolyDataMapper()
                mapper.SetInputData(paper)
                actor = vtk.vtkActor()
                actor.SetMapper(mapper)
                bufferPaper.AddActor(actor)

                pointActor = vtk.vtkActor()
                bufferPoints.RemoveAllViewProps()

                print("zeile 410")

                if True:
                    colors = vtk.vtkNamedColors()
                    Colors = vtk.vtkUnsignedCharArray()
                    Colors.SetNumberOfComponents(3)
                    Colors.SetName("Colors")
                    Colors.InsertNextTuple3(*colors.GetColor3ub('Red'))
                    Colors.InsertNextTuple3(*colors.GetColor3ub('Lime'))
                    Colors.InsertNextTuple3(*colors.GetColor3ub('Blue'))

                    # Create the topology of the point (a vertex)
                    new_points = vtk.vtkPoints()
                    vertices = vtk.vtkCellArray()
                    scalars = vtk.vtkDoubleArray()
                    scalars.SetNumberOfComponents(3)

                    id = new_points.InsertNextPoint(points.GetPoint(0))
                    vertices.InsertNextCell(1)
                    vertices.InsertCellPoint(id)

                    id = new_points.InsertNextPoint(points.GetPoint(1))
                    vertices.InsertNextCell(1)
                    vertices.InsertCellPoint(id)

                    id = new_points.InsertNextPoint(points.GetPoint(2))
                    vertices.InsertNextCell(1)
                    vertices.InsertCellPoint(id)

                    # Create a polydata object
                    point = vtk.vtkPolyData()

                    # Set the points and vertices we created as the geometry and topology of the polydata
                    point.SetPoints(new_points)
                    point.SetVerts(vertices)
                    point.GetPointData().SetScalars(Colors)

                    mapper = vtk.vtkPolyDataMapper()
                    mapper.SetInputData(point)

                    pointActor.SetMapper(mapper)
                    pointActor.GetProperty().SetPointSize(2)
                    bufferPoints.AddActor(pointActor)

                print("zeile 455")


                # render frame
                triangle = self.renderHelper(camera, buffer, bufferPaper, bufferPoints, bufferWin, i, actorList[meshNr])

                print("zeile 456")

                blue = triangle[:, :, 2]
                red = triangle[:, :, 0]
                green = triangle[:, :, 1]
                maskBlue = np.logical_and(np.logical_and(blue > 250, red < 1), green < 1)
                maskRed = np.logical_and(np.logical_and(red > 250, blue < 1), green < 1)
                maskGreen = np.logical_and(np.logical_and(green > 250, red < 1), blue < 1)

                uvArray.InsertNextTuple2(np.where(maskRed)[1][0] + img.shape[1],np.where(maskRed)[0][0])
                uvArray.InsertNextTuple2(np.where(maskGreen)[1][0] + img.shape[1],np.where(maskGreen)[0][0])
                uvArray.InsertNextTuple2(np.where(maskBlue)[1][0] + img.shape[1],np.where(maskBlue)[0][0])

                if len(img[0]) == 0:
                    img = triangle
                    black = np.zeros((300, img.shape[1], img.shape[2]))
                    img = np.vstack((img, black))
                else:
                    black = np.zeros((img.shape[0] - triangle.shape[0], triangle.shape[1], triangle.shape[2]))
                    triangle = np.vstack((triangle, black))
                    img = np.hstack((img, triangle))

                triCell = vtk.vtkTriangle()

                pointId = newPoints.InsertNextPoint(points.GetPoint(0))
                triCell.GetPointIds().SetId(0,pointId)

                pointId = newPoints.InsertNextPoint(points.GetPoint(1))
                triCell.GetPointIds().SetId(1,pointId)

                pointId = newPoints.InsertNextPoint(points.GetPoint(2))
                triCell.GetPointIds().SetId(2,pointId)

                newCells.InsertNextCell(triCell)
                #cellData.SetTuple(i, paperMeshColor)
                #paper.GetCellData().SetScalars(cellData)
                #paper.GetCellData().Modified()
                #paper.Modified()

                bufferPaper.RemoveAllViewProps()
                buffer.RemoveAllViewProps()

            mask = np.where(img > 0)
            width = mask[0].max() - mask[0].min()
            height = mask[1].max() - mask[1].min()
            if width > height:
                img = img[mask[0].min():mask[0].max(), mask[1].min(): mask[1].min() + width, :]
            else:
                img = img[mask[0].min():mask[0].min() + height, mask[1].min(): mask[1].max(), :]

            resultImg = vtk.vtkImageData()

            dy, dx, dz = img.shape

            vtkResult = numpy_support.numpy_to_vtk(img.reshape(dy * dx, dz))

            resultImg.SetSpacing(1., 1., 1.)
            resultImg.SetOrigin(0., 0., 0.)
            resultImg.SetDimensions(dx, dy, 1)
            resultImg.AllocateScalars(numpy_support.get_vtk_array_type(img.dtype), dz)
            resultImg.GetPointData().SetScalars(vtkResult)

            castFilter = vtk.vtkImageCast()
            castFilter.SetInputData(resultImg)
            castFilter.SetOutputScalarTypeToUnsignedChar()
            castFilter.Update()

            filename = os.path.join(self.dirname, "../out/2D/texture/texture{}.png".format(meshNr))
            writer = vtk.vtkPNGWriter()
            writer.SetFileName(filename)
            writer.SetInputConnection(castFilter.GetOutputPort())
            writer.Write()

            for i in range(uvArray.GetNumberOfTuples()):
                uvArray.SetTuple2(i,uvArray.GetTuple2(i)[0]/img.shape[1],uvArray.GetTuple2(i)[1]/img.shape[0])

            newGeometry.SetPoints(newPoints)
            newGeometry.SetPolys(newCells)
            newGeometry.GetPointData().SetTCoords(uvArray)

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputData(newGeometry)
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)

            texture = vtk.vtkTexture()
            texture.SetInputConnection(castFilter.GetOutputPort())
            actor.SetTexture(texture)

            self.dedicatedPaperMeshes[meshNr] = actor

            print("rendern geht")

            self.createUnfoldedPaperMesh(meshNr)

            print("unfold geht")

            meshNr += 1

        return self.dedicatedPaperMeshes

    def renderHelper(self, camera, buffer, bufferPaper, bufferPoints, bufferWin, count, actor):

        buffer.AddActor(actor)

        buffer.SetActiveCamera(camera)
        bufferPaper.SetActiveCamera(camera)
        bufferPoints.SetActiveCamera(camera)

        bufferWin.Render()

        wti = vtk.vtkWindowToImageFilter()
        wti.SetInput(bufferWin)
        wti.SetInputBufferTypeToRGB()
        wti.Update()

        return self.cropRenderedTriangle(wti.GetOutput(),count)

    def cropRenderedTriangle(self, image, count = 1):
        img = numpy_support.vtk_to_numpy(image.GetPointData().GetScalars())[:, 0:3]

        img.astype(float)
        print(img.shape)
        img = np.reshape(np.ravel(img), (500, 500, 3))


        blue = img[:,:,2]
        red = img[:,:,0]
        green = img[:,:,1]
        maskBlue = np.logical_and(np.logical_and(blue > 250, red < 1),green<1)
        maskRed = np.logical_and(np.logical_and(red > 250, blue < 1), green < 1)
        maskGreen = np.logical_and(np.logical_and(green > 250, red < 1),blue<1)
        mask = np.where((maskBlue+maskGreen+maskRed))

        if len(mask[0]) > 0 and len(mask[1]) > 0:
            result = img[mask[0].min()-2:mask[0].max()+2, mask[1].min()-2 : mask[1].max()+2, :]
        else:
            result = [[[1.0,0.0,0.0],[1.0,0.0,0.0],[1.0,0.0,0.0]],[[1.0,0.0,0.0],[1.0,0.0,0.0],[1.0,0.0,0.0]],[[1.0,0.0,0.0],[1.0,0.0,0.0],[1.0,0.0,0.0]]]
            result = np.array(result)

        return result

    def printSinglePoint(self,x,y,ren,color):
        newPoints = vtk.vtkPoints()
        vertices = vtk.vtkCellArray()

        p = [x,0.0,y]

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
        actor.GetProperty().SetColor(color)
        actor.GetProperty().SetPointSize(2)
        ren.AddActor(actor)

    #Method that maps the created texture to an mesh which is created according to the unfolded uv layout
    def createUnfoldedPaperMesh(self, idx):

        mesh = self.tempPaperActor.GetMapper().GetInput()

        width = 200

        textureCoordinates = mesh.GetPointData().GetTCoords()

        points = vtk.vtkPoints()
        cells = vtk.vtkCellArray()

        widthMarker = [0.0,0.0]
        heightMarker = [0.0,0.0]

        #creating a mesh according to the uv layout
        for i in range(mesh.GetNumberOfCells()):
            uvs = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
            textureCoordinates.GetTuple(mesh.GetCell(i).GetPointId(0), uvs[0])
            textureCoordinates.GetTuple(mesh.GetCell(i).GetPointId(1), uvs[1])
            textureCoordinates.GetTuple(mesh.GetCell(i).GetPointId(2), uvs[2])

            cells.InsertNextCell(3)

            for j in range(len(uvs)):
                x = (uvs[j][0] * width)-(width/2.0)
                y = (uvs[j][1] * width)-(width/2.0)

                point = [x,-1.0,y]

                #x += width/2.0
                #y += width/2.0

                if x>widthMarker[1]: widthMarker[1] = x
                elif x<widthMarker[0]: widthMarker[0] = x
                if y>heightMarker[1]: heightMarker[1] = y
                elif y<heightMarker[0]: heightMarker[0] = y

                id = points.InsertNextPoint(point)
                cells.InsertCellPoint(id)

        textureCoordinates = self.dedicatedPaperMeshes[idx].GetMapper().GetInput().GetPointData().GetTCoords()

        unfoldedPaper = vtk.vtkPolyData()

        unfoldedPaper.SetPoints(points)
        unfoldedPaper.SetPolys(cells)
        unfoldedPaper.GetPointData().SetTCoords(textureCoordinates)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(unfoldedPaper)

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        #----------------------

        filename = os.path.join(self.dirname, "../out/2D/texture/texture{}.png".format(idx))
        readerFac = vtk.vtkImageReader2Factory()
        imageReader = readerFac.CreateImageReader2(filename)
        imageReader.SetFileName(filename)

        texture = vtk.vtkTexture()
        texture.SetInputConnection(imageReader.GetOutputPort())

        actor.SetTexture(texture)

        camera = vtk.vtkCamera()
        camera.SetPosition(0,-600,0)
        camera.SetFocalPoint(0,0,0)
        camera.SetViewUp(0,0,1)
        ren, iren, renWin, wti = util.getbufferRenIntWin(camera,width=3000,height=3000)
        #print(wti.GetViewport())

        self.printSinglePoint(widthMarker[0],heightMarker[0],ren,[0.0,0.0,255.0])
        self.printSinglePoint(widthMarker[1],heightMarker[1],ren,[0.0,255.0,0.0])

        ren.AddActor(actor)
        renWin.Render()

        print("zeile 725")

        wti.Update()

        print()
        imgWidth = wti.GetOutput().GetExtent()[1] + 1
        imgHeight = wti.GetOutput().GetExtent()[3] + 1
        print(imgWidth, " ", imgHeight)

        img = self.vtkToNpHelper(wti.GetOutput(), imgWidth, imgHeight)

        #castFilter = vtk.vtkImageCast()
        #castFilter.SetInputConnection(wti.GetOutputPort())
        #castFilter.SetOutputScalarTypeToUnsignedChar()
        #castFilter.Update()

        print("zeile 752")

        resultImg = vtk.vtkImageData()

        dy, dx, dz = img.shape

        vtkResult = numpy_support.numpy_to_vtk(img.reshape(dy * dx, dz))

        resultImg.SetSpacing(1., 1., 1.)
        resultImg.SetOrigin(0., 0., 0.)
        resultImg.SetDimensions(dx, dy, 1)
        resultImg.AllocateScalars(numpy_support.get_vtk_array_type(img.dtype), dz)
        resultImg.GetPointData().SetScalars(vtkResult)

        print("zeile 766")

        castFilter = vtk.vtkImageCast()
        castFilter.SetInputData(resultImg)
        castFilter.SetOutputScalarTypeToUnsignedChar()
        castFilter.Update()

        print("zeile 773")

        writer = vtk.vtkPNGWriter()
        filename = os.path.join(self.dirname, "../out/2D/unfolding{}.png".format(idx))

        writer.SetFileName(filename)
        writer.SetInputConnection(castFilter.GetOutputPort())
        writer.Write()

    def vtkToNpHelper(self, img, width, height):

        img = numpy_support.vtk_to_numpy(img.GetPointData().GetScalars())[:, 0:3]
        print(img.shape)
        img.astype(float)

        print("zeile 740")

        print(img.shape)
        img = np.reshape(np.ravel(img), (width, height, 3))

        print("zeile 742")

        mask = np.where(img != 255)
        width = mask[0].max() - mask[0].min()
        height = mask[1].max() - mask[1].min()
        if width > height:
            img = img[mask[0].min():mask[0].max(), mask[1].min(): mask[1].min() + width, :]
        else:
            img = img[mask[0].min():mask[0].min() + height, mask[1].min(): mask[1].max(), :]

        return img

    def booleanTrimesh(self,ren,graph,actorlist,nestedlist):

        cubeactor = self.createPapermesh(nestedlist,graph,True)
        #binder.onUnfold()
        #cubeactor = self.importUnfoldedMesh(True)

        poly = self.tempPaperActor.GetMapper().GetInput()

        clip = vtk.vtkClipPolyData()
        clip.GenerateClipScalarsOn()
        clip.GenerateClippedOutputOn()
        clip.SetInputData(poly)

        plane = vtk.vtkPlane()
        plane.SetOrigin(cubeactor.GetMapper().GetInput().GetCenter())
        plane.SetNormal(0.0, 0.0, 1.0)

        clip.SetClipFunction(plane)
        clip.SetValue(0)
        clip.Update()

        fillHoles = vtk.vtkFillHolesFilter()
        fillHoles.SetInputData(clip.GetOutput())
        fillHoles.SetHoleSize(10000.0)
        fillHoles.Update()

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(fillHoles.GetOutput())
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

    #----------------------
        clip2 = vtk.vtkClipPolyData()
        clip2.GenerateClipScalarsOn()
        clip2.GenerateClippedOutputOn()
        clip2.SetInputData(poly)

        plane2 = vtk.vtkPlane()
        plane2.SetOrigin(cubeactor.GetMapper().GetInput().GetCenter())
        plane2.SetNormal(0.0, 0.0, -1.0)

        clip2.SetClipFunction(plane2)
        clip2.SetValue(0)
        clip2.Update()

        fillHoles2 = vtk.vtkFillHolesFilter()
        fillHoles2.SetInputData(clip2.GetOutput())
        fillHoles2.SetHoleSize(10000.0)
        fillHoles2.Update()

        mapper2 = vtk.vtkPolyDataMapper()
        mapper2.SetInputData(fillHoles2.GetOutput())
        actor2 = vtk.vtkActor()
        actor2.SetMapper(mapper2)
    #-------------------------------

        filename = os.path.join(self.dirname, "../out/3D/boolMesh.stl")
        stlWriter = vtk.vtkSTLWriter()
        stlWriter.SetFileName(filename)
        stlWriter.SetInputData(actor.GetMapper().GetInput())
        stlWriter.Write()

        filename = os.path.join(self.dirname, "../out/3D/boolMesh2.stl")
        stlWriter = vtk.vtkSTLWriter()
        stlWriter.SetFileName(filename)
        stlWriter.SetInputData(actor2.GetMapper().GetInput())
        stlWriter.Write()

        filename = os.path.join(self.dirname, "../out/3D/cutout.stl")
        stlWriter = vtk.vtkSTLWriter()
        stlWriter.SetFileName(filename)
        stlWriter.SetInputData(cubeactor.GetMapper().GetInput())
        stlWriter.Write()

        filename = os.path.join(self.dirname, "../out/3D/boolMesh.stl")
        mesh = trimesh.load(filename)

        filename = os.path.join(self.dirname, "../out/3D/boolMesh2.stl")
        mesh2 = trimesh.load(filename)

        filename = os.path.join(self.dirname, "../out/3D/cutout.stl")
        cutout = trimesh.load(filename)

        print("imports working")

        boolUpActor, boolDownActor = self.trimeshBooleanBugWorkaround(mesh,mesh2,cutout)

        writer = vtk.vtkSTLWriter()
        writer.SetInputData(boolUpActor.GetMapper().GetInput())
        filename = os.path.join(self.dirname, "../out/3D/upper.stl")
        writer.SetFileName(filename)
        writer.Write()
#
#        binder.onUnfold()
#        unfoldedUpActor = self.importUnfoldedMesh(True)
#
#        print("unfolded the Upper Actor")
#
        writer = vtk.vtkSTLWriter()
        writer.SetInputData(boolDownActor.GetMapper().GetInput())
        filename = os.path.join(self.dirname, "../out/3D/lower.stl")
        writer.SetFileName(filename)
        writer.Write()
#
#        binder.onUnfold()
#        unfoldedDownActor = self.importUnfoldedMesh(True)
#
#        print("unfolded the Lower Actor")

        self.unfoldTest()

        self.dedicatedPaperMeshes = [boolUpActor,boolDownActor]

        ren.RemoveAllViewProps()
        ren.AddActor(boolUpActor)
        ren.AddActor(boolDownActor)

    def unfoldTest(self):

        filename = os.path.join(self.dirname, "../out/3D/upper.stl")
        #        filename = os.path.join(self.dirname, "../out/3D/papermesh_cleaned.off")

        mesh = meshio.read(
            filename,  # string, os.PathLike, or a buffer/open file
            file_format="stl",  # optional if filename is a path; inferred from extension
        )

        filename = os.path.join(self.dirname, "../out/3D/papermesh.off")

        meshio.write(
            filename,  # str, os.PathLike, or buffer/ open file
            mesh,
            # file_format="vtk",  # optional if first argument is a path; inferred from extension
        )

        #       print(filename)
        graph = Graph()

        graph.load(filename)
        if not graph.unfold(5000, 0):
            print("failed to unfold :(")
        else:
            print("succesfully unfolded :)")
            filename = os.path.join(self.dirname, "../out/3D/unfolded/upper.obj")
            gluetabs_filename = os.path.join(self.dirname, "../out/3D/unfolded/gluetabs_upper.obj")

            graph.save(filename, gluetabs_filename)

    def trimeshBooleanBugWorkaround(self,mesh,mesh2,cutout):

        volumeSum = 0
        boolDownActor = vtk.vtkActor()
        boolUpActor = vtk.vtkActor()

        for i in range(10):
            # ---------------------------

            # meshes = [mesh,cutout]
            # meshes = []
            # meshes.append(mesh)
            # meshes.append(cutout)
            # boolMesh = trimesh.boolean.difference(meshes,engine='scad')
            boolMesh = mesh.difference(cutout, engine='blender')

            print("after first bool")

            filename = os.path.join(self.dirname, "../out/3D/bool.stl")
            boolMesh.export(filename)

            reader = vtk.vtkSTLReader()
            reader.SetFileName(filename)
            reader.Update()

            # -------------------

            # meshes2 = [mesh2,cutout]
            # meshes2 = []
            # meshes2.append(mesh2)
            # meshes2.append(cutout)
            # boolMesh2 = trimesh.boolean.difference(meshes2,engine='blender')
            boolMesh2 = mesh2.difference(cutout, engine='blender')

            filename = os.path.join(self.dirname, "../out/3D/bool2.stl")
            boolMesh2.export(filename)

            reader2 = vtk.vtkSTLReader()
            reader2.SetFileName(filename)
            reader2.Update()

            mass = vtk.vtkMassProperties()
            mass.SetInputData(reader.GetOutput())
            mass.Update()
            volume = mass.GetVolume()

            mass2 = vtk.vtkMassProperties()
            mass2.SetInputData(reader2.GetOutput())
            mass2.Update()
            volume2 = mass2.GetVolume()

            if volumeSum < (volume + volume2):

                print(volumeSum)

                volumeSum = (volume + volume2)

                #-------------------

                #translation = vtk.vtkTransform()
                #translation.Translate(0.0, 0.0, 100.0)

                #transformFilter = vtk.vtkTransformFilter()
                #transformFilter.SetInputData(reader.GetOutput())
                #transformFilter.SetTransform(translation)
                #transformFilter.Update()
                # boolUpActor.GetProperty().SetOpacity(0.5)

                triangle = vtk.vtkTriangleFilter()
                triangle.SetInputData(reader.GetOutput())
                triangle.Update()

                boolUpMapper = vtk.vtkPolyDataMapper()
                boolUpMapper.SetInputData(triangle.GetOutput())
                boolUpMapper.Update()
                boolUpActor.SetMapper(boolUpMapper)

                #--------------

                #translation2 = vtk.vtkTransform()
                #translation2.Translate(0.0, 0.0, -100.0)

                #transformFilter2 = vtk.vtkTransformFilter()
                #transformFilter2.SetInputData(reader2.GetOutput())
                #transformFilter2.SetTransform(translation2)
                #transformFilter2.Update()
                # boolDownActor.GetProperty().SetOpacity(0.5)

                triangle2 = vtk.vtkTriangleFilter()
                triangle2.SetInputData(reader2.GetOutput())
                triangle2.Update()

                boolDownMapper = vtk.vtkPolyDataMapper()
                boolDownMapper.SetInputData(triangle2.GetOutput())
                boolDownMapper.Update()
                boolDownActor.SetMapper(boolDownMapper)

        return boolUpActor,boolDownActor
import vtkmodules.all as vtk
import numpy as np
import os
import meshInteraction
import trimesh
import meshio
import projector
from PyQt5 import QtWidgets
from mu3d import Graph


##Class responsible for mesh related processing steps and I/O.
class MeshProcessing():

    dirname = os.path.dirname(__file__)

    dedicatedPaperMeshes = []
    tempPaperActor = vtk.vtkActor()

    projector = projector.Projector()
    meshInteractor = meshInteraction.MeshInteraction(dedicatedPaperMeshes)

    ##Create initial Mesh to subdivide and wrap
    def createPapermesh(self,actorList, boolReturn = False):

        polysAppended = self.appendAllMeshes(actorList)

        hull = vtk.vtkHull()
        hull.SetInputData(polysAppended)
        hull.AddCubeFacePlanes()
        # hull.AddRecursiveSpherePlanes(1)
        hull.Update()

        triangleFilter = vtk.vtkTriangleFilter()
        triangleFilter.SetInputConnection(hull.GetOutputPort())
        triangleFilter.Update()

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

        paperMesh = self.calcMeshNormals(paperMesh)
        #paperMesh.GetPointData().SetNormals(normals)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(paperMesh)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        if not boolReturn:
            self.tempPaperActor = actor
        return actor

    ##Calls the mu3d framework to unfold the mesh given with hardcodet filename, could be changed todo.
    def mu3dUnfoldPaperMesh(self, actor, graph, iterations):
        filename = os.path.join(self.dirname, "../out/3D/papermesh.stl")
        # objWriter.SetFilePrefix(filename)
        # objWriter.Write()

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
        if not graph.unfold(iterations, 0):
            msgBox = QtWidgets.QMessageBox()
            msgBox.setText("failed to unfold :( in {} iterations".format(iterations))
            msgBox.exec()
            return False
        else:
            print("succesfully unfolded :) in {} iterations".format(iterations))

            filename = os.path.join(self.dirname, "../out/3D/unfolded/model.obj")
            gluetabs_filename = os.path.join(self.dirname, "../out/3D/unfolded/gluetabs.obj")

            graph.save(filename, gluetabs_filename)

            importer = vtk.vtkOBJReader()
            filename = os.path.join(self.dirname, "../out/3D/unfolded/model.obj")
            importer.SetFileName(filename)
            importer.Update()

            mesh = importer.GetOutput()
            mesh = self.normalizeUV(mesh)

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputData(mesh)
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)

            actor.GetProperty().SetColor([0.5, 0.5, 0.5])
            actor.GetProperty().BackfaceCullingOn()
            actor.GetProperty().SetOpacity(0.5)

            self.tempPaperActor = actor

            return True


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

    ##For each loaded structure shrinkwrap the global paper mesh onto them
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

    ##Imports a previous unfolded mesh, with hardcodet filename, could be changed todo.
    def importUnfoldedMesh(self, boolReturn = False):

        importer = vtk.vtkOBJReader()

        #hardcodet file path todo
        filename = os.path.join(self.dirname, "../out/3D/unfolded/model.obj")
        importer.SetFileName(filename)

        importer.Update()
        mesh = importer.GetOutput()

        mesh = self.normalizeUV(mesh)

        mesh = self.calcMeshNormals(mesh)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(mesh)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        actor.GetProperty().SetColor([0.5,0.5,0.5])
        actor.GetProperty().BackfaceCullingOn()

        if not boolReturn:
            self.tempPaperActor = actor
        return actor

    def normalizeUV(self,mesh):
        textureCoordinates = mesh.GetPointData().GetTCoords()
        absmax = 0.0
        normmax = 0.0
        newTCoords = vtk.vtkFloatArray()
        newTCoords.SetNumberOfComponents(2)

        # Find maxima to normalize uvs.
        for i in range(int(textureCoordinates.GetNumberOfTuples() / 3)):
            uvs = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
            textureCoordinates.GetTuple(mesh.GetCell(i).GetPointId(0), uvs[0])
            textureCoordinates.GetTuple(mesh.GetCell(i).GetPointId(1), uvs[1])
            textureCoordinates.GetTuple(mesh.GetCell(i).GetPointId(2), uvs[2])
            abstemp = np.amax(abs(np.array(uvs)))
            normtemp = np.amax((np.array(uvs)))
            if abstemp > absmax: absmax = abstemp
            if normtemp > normmax: normmax = normtemp

        # Normalize the uvs.
        for i in range(int(textureCoordinates.GetNumberOfTuples() / 3)):
            u = (textureCoordinates.GetTuple2((mesh.GetCell(i).GetPointId(0)))[0] + absmax) / (absmax + normmax)
            v = (textureCoordinates.GetTuple2((mesh.GetCell(i).GetPointId(0)))[1] + absmax) / (absmax + normmax)
            newTCoords.InsertNextTuple2(u, v)

            u = (textureCoordinates.GetTuple2((mesh.GetCell(i).GetPointId(1)))[0] + absmax) / (absmax + normmax)
            v = (textureCoordinates.GetTuple2((mesh.GetCell(i).GetPointId(1)))[1] + absmax) / (absmax + normmax)
            newTCoords.InsertNextTuple2(u, v)

            u = (textureCoordinates.GetTuple2((mesh.GetCell(i).GetPointId(2)))[0] + absmax) / (absmax + normmax)
            v = (textureCoordinates.GetTuple2((mesh.GetCell(i).GetPointId(2)))[1] + absmax) / (absmax + normmax)
            newTCoords.InsertNextTuple2(u, v)

        mesh.GetPointData().SetTCoords(newTCoords)
        return mesh


    def calcMeshNormals(self,polydata):
        #normals = vtk.vtkTriangleMeshPointNormals()

        normals = vtk.vtkPolyDataNormals()
        normals.ComputePointNormalsOff()
        normals.ComputeCellNormalsOn()
        normals.ConsistencyOn()
        #normals.AutoOrientNormalsOn()

        normals.SetInputData(polydata)
        normals.Update()
        #return normals.GetOutput().GetPointData().GetNormals()
        return normals.GetOutput()

    def project(self,inflateStruc,actorList, resolution):
        meshes = []
        idx = 0
        for a in self.dedicatedPaperMeshes:
            mesh = (self.projector.projectPerTriangle(a,inflateStruc,actorList,idx, resolution))
            meshes.append(mesh)
            self.projector.createUnfoldedPaperMesh(mesh, self.tempPaperActor, idx)
            idx+=1
        return meshes
    
    def booleanTrimesh(self,ren,graph,actorlist,nestedlist):

        cubeactor = self.createPapermesh(nestedlist,True)
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

        filledMesh = fillHoles.GetOutput()
        filledMesh = self.calcMeshNormals(filledMesh)
        #filledMesh.GetPointData().SetNormals(normals2)
        filledMesh.GetPoints().GetData().Modified()

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(filledMesh)
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

        filledMesh2 = fillHoles2.GetOutput()
        filledMesh2 = self.calcMeshNormals(filledMesh2)
        #filledMesh2.GetPointData().SetNormals(normals2)
        filledMesh2.GetPoints().GetData().Modified()

        mapper2 = vtk.vtkPolyDataMapper()
        mapper2.SetInputData(filledMesh2)
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

        filename = os.path.join(self.dirname, "../out/3D/boolMeshTEST.stl")
        mesh = trimesh.load(filename)
        filename = os.path.join(self.dirname, "../out/3D/mesh.off")
        trimesh.exchange.export.export_mesh(mesh, filename, file_type="off")

        filename = os.path.join(self.dirname, "../out/3D/boolMesh2TEST.stl")
        mesh2 = trimesh.load(filename)

        filename = os.path.join(self.dirname, "../out/3D/cutout.stl")
        cutout = trimesh.load(filename)
        filename = os.path.join(self.dirname, "../out/3D/cutout.off")
        trimesh.exchange.export.export_mesh(cutout, filename, file_type="off")

        print("imports working")

        boolUpActor, boolDownActor = self.trimeshBooleanBugWorkaround(mesh,mesh2,cutout)

        writer = vtk.vtkSTLWriter()
        writer.SetInputData(self.decimateMesh(boolUpActor.GetMapper().GetInput()))
        filename = os.path.join(self.dirname, "../out/3D/upper.stl")
        writer.SetFileName(filename)
        writer.Write()

        #TEST
        writer = vtk.vtkSTLWriter()
        writer.SetInputData(self.decimateMesh(boolUpActor.GetMapper().GetInput()))
        filename = os.path.join(self.dirname, "../out/3D/upper_decimated.stl")
        writer.SetFileName(filename)
        writer.Write()

#
#        binder.onUnfold()
#        unfoldedUpActor = self.importUnfoldedMesh(True)
#
#        print("unfolded the Upper Actor")
#
        writer = vtk.vtkSTLWriter()
        writer.SetInputData(self.decimateMesh(boolDownActor.GetMapper().GetInput()))
        filename = os.path.join(self.dirname, "../out/3D/lower.stl")
        writer.SetFileName(filename)
        writer.Write()
#
#        binder.onUnfold()
#        unfoldedDownActor = self.importUnfoldedMesh(True)
#
#        print("unfolded the Lower Actor")

        self.unfoldTest()

        importer = vtk.vtkOBJReader()
        filename = os.path.join(self.dirname,  "../out/3D/unfolded/upper.obj")
        importer.SetFileName(filename)
        importer.Update()
        upper = importer.GetOutput()

        print("unfolded upper")

        self.unfoldTest(name="lower")

        importer = vtk.vtkOBJReader()
        filename = os.path.join(self.dirname,  "../out/3D/unfolded/lower.obj")
        importer.SetFileName(filename)
        importer.Update()
        lower = importer.GetOutput()

        print("unfolded lower")

        self.dedicatedPaperMeshes = [upper,lower]

        ren.RemoveAllViewProps()
        ren.AddActor(boolUpActor)
        ren.AddActor(boolDownActor)

    def unfoldTest(self, name = "upper"):

        filename = os.path.join(self.dirname, "../out/3D/" + name + ".stl")
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
        if not graph.unfold(50000, 0):
            print("failed to unfold :(")
        else:
            print("succesfully unfolded :)")
            filename = os.path.join(self.dirname, "../out/3D/unfolded/" + name + ".obj")
            gluetabs_filename = os.path.join(self.dirname, "../out/3D/unfolded/gluetabs_" + name + ".obj")

            graph.save(filename, gluetabs_filename)

    def decimateMesh(self,mesh):

        subdivider = vtk.vtkLoopSubdivisionFilter()
        subdivider.SetInputData(mesh)
        subdivider.SetNumberOfSubdivisions(1)
        subdivider.Update()

        decimater = vtk.vtkDecimatePro()
        decimater.SetInputData(subdivider.GetOutput())
        decimater.SetTargetReduction(0.1)
        decimater.PreserveTopologyOn()
        decimater.Update()

        quadrClustering = vtk.vtkQuadricClustering()
        quadrClustering.SetInputData(decimater.GetOutput())
        quadrClustering.SetNumberOfXDivisions(32)
        quadrClustering.SetNumberOfYDivisions(32)
        quadrClustering.SetNumberOfZDivisions(4)
        quadrClustering.Update()
        return quadrClustering.GetOutput()

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
            #boolMesh = trimesh.repair.fix_normals(boolMesh)
            #boolMesh = trimesh.repair.fix_winding(boolMesh)
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
            boolMesh2 = mesh2.difference(cutout, engine='scad')
            #boolMesh2 = trimesh.repair.fix_normals(boolMesh2)
            #boolMesh2 = trimesh.repair.fix_winding(boolMesh2)

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

                upMesh = triangle.GetOutput()
                upMesh = self.calcMeshNormals(upMesh)
                #upMesh.GetPointData().SetNormals(normals)

                boolUpMapper = vtk.vtkPolyDataMapper()
                boolUpMapper.SetInputData(upMesh)
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

                bottomMesh = triangle2.GetOutput()
                bottomMesh = self.calcMeshNormals(bottomMesh)
                #bottomMesh.GetPointData().SetNormals(normals2)

                boolDownMapper = vtk.vtkPolyDataMapper()
                boolDownMapper.SetInputData(bottomMesh)
                boolDownMapper.Update()
                boolDownActor.SetMapper(boolDownMapper)

        return boolUpActor,boolDownActor

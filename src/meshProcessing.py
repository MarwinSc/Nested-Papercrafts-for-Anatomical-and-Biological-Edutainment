import projector as projector
import vtkmodules.all as vtk
import numpy as np
import os
import meshInteraction
import meshio
import util
from boolean import boolean_interface
from mu3d.mu3d import Graph
from PyQt5 import QtWidgets

#Class responsible for mesh related processing steps and I/O.
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
        hull.Update()

        triangleFilter = vtk.vtkTriangleFilter()
        triangleFilter.SetInputData(hull.GetOutput())
        triangleFilter.Update()

        mesh = util.subdivideMesh(triangleFilter.GetOutput())

        mesh = util.shrinkWrap(mesh,polysAppended)
        paperMesh = util.offsetMesh(mesh)

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
            mesh = self.shiftUVsToOrigin(mesh)
            mesh = self.normalizeUV(mesh)
            mesh = self.calcMeshNormals(mesh)

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputData(mesh)
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)

            #actor.GetProperty().SetColor([1.0, 1.0, 1.0])
            actor.GetProperty().BackfaceCullingOn()
            actor.GetProperty().SetOpacity(0.5)

            self.tempPaperActor = actor

            #just to write the model with normalized uvs
            filename = os.path.join(self.dirname, "../out/3D/unfolded/model.obj")
            objWriter = vtk.vtkOBJWriter()
            objWriter.SetInputData(actor.GetMapper().GetInput())
            objWriter.SetFileName(filename)
            objWriter.Write()

            return True

    def appendAllMeshes(self, actorList):

        append = vtk.vtkAppendPolyData()
        clean = vtk.vtkCleanPolyData()

        for a in actorList:
            append.AddInputData(a.GetMapper().GetInput())

        clean.SetInputConnection(append.GetOutputPort())
        clean.Update()

        return clean.GetOutput()

    ##For each loaded structure shrinkwrap the global paper mesh onto them
    def createDedicatedMeshes(self,inflateStrucList,actorList):

        for a in range(len(actorList)):

            if inflateStrucList[a]:

                mesh = self.tempPaperActor.GetMapper().GetInput()

                poly = util.shrinkWrap(mesh, actorList[a].GetMapper().GetInput(),False)
                #poly = util.smoothMesh(mesh,actorList[a].GetMapper().GetInput(),20,0.1)

            else:
                poly = self.tempPaperActor.GetMapper().GetInput()

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputData(poly)

            actor = vtk.vtkActor()
            actor.SetMapper(mapper)

            self.dedicatedPaperMeshes.append(actor)

    ##Imports a previous unfolded mesh, with hardcodet filename, could be changed todo.
    def importUnfoldedMesh(self, name, boolReturn = False):

        importer = vtk.vtkOBJReader()

        #hardcodet file path todo
        filename = os.path.join(self.dirname, "../out/3D/unfolded/"+name+".obj")
        importer.SetFileName(filename)

        importer.Update()
        mesh = importer.GetOutput()

        mesh = self.shiftUVsToOrigin(mesh)
        mesh = self.normalizeUV(mesh)
        mesh = self.calcMeshNormals(mesh)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(mesh)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        # actor.GetProperty().SetColor([1.0,1.0,1.0])
        actor.GetProperty().BackfaceCullingOn()
        actor.GetProperty().SetOpacity(0.5)

        if not boolReturn:
            self.tempPaperActor = actor
        return actor

    def normalizeUV(self,mesh):
        textureCoordinates = mesh.GetPointData().GetTCoords()
        newTCoords = vtk.vtkFloatArray()
        newTCoords.SetNumberOfComponents(2)

        gmin, gmax = self.getMinMaxUV(mesh)

        # Normalize the uvs.
        for i in range(int(textureCoordinates.GetNumberOfTuples() / 3)):
            for j in range(3):
                u = ((textureCoordinates.GetTuple2((mesh.GetCell(i).GetPointId(j)))[0])-gmin) / (gmax-gmin)
                v = ((textureCoordinates.GetTuple2((mesh.GetCell(i).GetPointId(j)))[1])-gmin) / (gmax-gmin)
                newTCoords.InsertNextTuple2(u, v)

        mesh.GetPointData().SetTCoords(newTCoords)
        return mesh

    def getMinMaxUV(self,mesh, axis = False):
        textureCoordinates = mesh.GetPointData().GetTCoords()

        if(axis):
            uMax = 0
            uMin = 0
            vMax = 0
            vMin = 0

            for i in range(int(textureCoordinates.GetNumberOfTuples() / 3)):
                uvs = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
                textureCoordinates.GetTuple(mesh.GetCell(i).GetPointId(0), uvs[0])
                textureCoordinates.GetTuple(mesh.GetCell(i).GetPointId(1), uvs[1])
                textureCoordinates.GetTuple(mesh.GetCell(i).GetPointId(2), uvs[2])

                #normtemp = np.amax((np.array(uvs)))
                umax = np.amax(np.array(uvs)[:,0])
                umin = np.amin(np.array(uvs)[:,0])
                vmax = np.amax(np.array(uvs)[:,1])
                vmin = np.amin(np.array(uvs)[:,1])

                #abstemp = np.amax(abs(np.array(uvs)))
                if umax > uMax: uMax = umax
                if umin < uMin: uMin = umin
                if vmax > vMax: vMax = vmax
                if vmin < vMin: vMin = vmin

            return uMin,uMax,vMin,vMax

        else:
            gmax = 0.0
            gmin = 0.0
            # Find maxima to normalize uvs.
            for i in range(int(textureCoordinates.GetNumberOfTuples() / 3)):
                uvs = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
                textureCoordinates.GetTuple(mesh.GetCell(i).GetPointId(0), uvs[0])
                textureCoordinates.GetTuple(mesh.GetCell(i).GetPointId(1), uvs[1])
                textureCoordinates.GetTuple(mesh.GetCell(i).GetPointId(2), uvs[2])

                #normtemp = np.amax((np.array(uvs)))

                max = np.amax(np.array(uvs))
                min = np.amin(np.array(uvs))
                #abstemp = np.amax(abs(np.array(uvs)))
                if max > gmax: gmax = max
                if min < gmin: gmin = min

            return gmin, gmax

    def shiftUVsToOrigin(self,mesh):
        uMin, uMax, vMin, vMax = self.getMinMaxUV(mesh, axis=True)
        textureCoordinates = mesh.GetPointData().GetTCoords()
        newTCoords = vtk.vtkFloatArray()
        newTCoords.SetNumberOfComponents(2)

        for i in range(int(textureCoordinates.GetNumberOfTuples() / 3)):
            for j in range(3):
                if uMin > 0:
                    u = ((textureCoordinates.GetTuple2((mesh.GetCell(i).GetPointId(j)))[0]) - uMin)
                else:
                    u = ((textureCoordinates.GetTuple2((mesh.GetCell(i).GetPointId(j)))[0]) + abs(uMin))

                if vMin > 0:
                    v = ((textureCoordinates.GetTuple2((mesh.GetCell(i).GetPointId(j)))[1]) - vMin)
                else:
                    v = ((textureCoordinates.GetTuple2((mesh.GetCell(i).GetPointId(j)))[1]) + abs(vMin))
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
            try:
                mesh = (self.projector.projectPerTriangle(a,inflateStruc,actorList,idx, resolution))
                meshes.append(mesh)
                self.projector.createUnfoldedPaperMesh(mesh, self.tempPaperActor, idx)
            except:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setText("One or more meshes could not be projected.")
                msgBox.exec()
            idx+=1
        return meshes

    def cutMeshWithPlanes(self,mesh,cutPlanes):

        results = []

        for plane in cutPlanes:

            clip = vtk.vtkClipPolyData()
            clip.GenerateClipScalarsOn()
            clip.GenerateClippedOutputOn()
            clip.SetInputData(mesh)

            clip.SetClipFunction(plane)
            clip.SetValue(0)
            clip.Update()

            fillHoles = vtk.vtkFillHolesFilter()
            fillHoles.SetInputData(clip.GetOutput())
            fillHoles.SetHoleSize(10000.0)
            fillHoles.Update()

            filledMesh = fillHoles.GetOutput()
            filledMesh = self.calcMeshNormals(filledMesh)
            # filledMesh.GetPointData().SetNormals(normals2)
            filledMesh.GetPoints().GetData().Modified()

            results.append(filledMesh)

        return results

    def booleanCGAL(self,nestedlist):

        mesh = self.tempPaperActor.GetMapper().GetInput()
        #todo change to list or hierachical logic
        cutout = self.createPapermesh(nestedlist, True)
        cutPlanes = []

        #for the "upper" part of the mesh
        plane = vtk.vtkPlane()
        plane.SetOrigin(cutout.GetMapper().GetInput().GetCenter())
        plane.SetNormal(0.0, 0.0, 1.0)
        cutPlanes.append(plane)

        #for the "lower" part of the mesh
        plane = vtk.vtkPlane()
        plane.SetOrigin(cutout.GetMapper().GetInput().GetCenter())
        plane.SetNormal(0.0, 0.0, -1.0)
        cutPlanes.append(plane)

        meshPieces = self.cutMeshWithPlanes(mesh, cutPlanes)

        util.writeStl(cutout.GetMapper().GetInput(), "cutout")

        for i in range(len(meshPieces)):
            # write stl and convert to off
            util.writeStl(meshPieces[i],"boolMesh{}".format(i))
            inPath = os.path.join(self.dirname, "../out/3D/boolMesh{}.stl".format(i))
            outPath = os.path.join(self.dirname, "../out/3D/boolMesh{}.off".format(i))
            util.meshioIO(inPath,outPath)

        #right now only for testing todo write and bool all meshes
        cgal = boolean_interface.Boolean_Interface()
        mesh = os.path.join(self.dirname, "../out/3D/boolMesh1.off")
        cutout = os.path.join(self.dirname, "../out/3D/cutout.off")

        #todo split boolean and unfolding so first all meshes are calculated and afterwards the meshes are unfolded on user input?
        cgal.boolean(mesh, cutout)
        graph = Graph()
        filename = os.path.join(self.dirname, "difference.off")
        graph.load(filename)
        if not graph.unfold(50000, 0):
            print("failed to unfold :(")
        else:
            print("succesfully unfolded :)")
            filename = os.path.join(self.dirname, "../out/3D/unfolded/difference.obj")
            gluetabs_filename = os.path.join(self.dirname, "../out/3D/unfolded/gluetabs_difference.obj")
            graph.save(filename, gluetabs_filename)

        #set as dedicated meshes for unfolding
        #self.dedicatedPaperMeshes = [upper,lower]

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

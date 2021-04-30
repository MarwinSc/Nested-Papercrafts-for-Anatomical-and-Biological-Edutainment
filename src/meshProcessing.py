import vtkmodules.all as vtk
import numpy as np
import os
from mesh import Mesh
import util
from boolean import boolean_interface
from mu3d.mu3dpy.mu3d import Graph
from PyQt5 import QtWidgets

class MeshProcessing():
    '''
    Class responsible for mesh related processing steps and I/O.
    '''

    dirname = os.path.dirname(__file__)

    tempPaperActor = vtk.vtkActor()

    #meshInteractor = meshInteraction.MeshInteraction(dedicatedPaperMeshes)

    def mu3dUnfoldPaperMesh(self, actor, graph, iterations):
        '''
        Forwards the mesh to the mu3d wrapper to unfold it.
        :param actor: The vtk actor containing the mesh to unfold.
        :param graph: The wrapped mu3d graph object.
        :param iterations: The iterations for the unfolding.
        :return: If the unfolding is successful the vtk actor containing the unfolded mesh is returned.
        '''
        inpath = os.path.join(self.dirname, "../out/3D/papermesh.stl")
        outpath = os.path.join(self.dirname, "../out/3D/papermesh.off")
        util.writeStl(actor.GetMapper().GetInput(),"papermesh")

        util.meshioIO(inpath,outpath)

        graph.load(outpath)
        if not graph.unfold(iterations, 0):
            msgBox = QtWidgets.QMessageBox()
            msgBox.setText("failed to unfold :( in {} iterations".format(iterations))
            msgBox.exec()
            return None
        else:
            print("succesfully unfolded :) in {} iterations".format(iterations))

            filename = os.path.join(self.dirname, "../out/3D/unfolded/model.obj")
            gluetabs_filename = os.path.join(self.dirname, "../out/3D/unfolded/gluetabs.obj")

            graph.save(filename, gluetabs_filename)

            filename = os.path.join(self.dirname, "../out/3D/unfolded/model.obj")
            mesh = util.readObj(filename)
            mesh = self.normalizeUV(mesh)
            mesh = self.calcMeshNormals(mesh)

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputData(mesh)
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)

            #actor.GetProperty().SetColor([1.0, 1.0, 1.0])
            actor.GetProperty().BackfaceCullingOn()
            actor.GetProperty().SetOpacity(0.5)

            #just to write the model with normalized uvs
            util.writeObj(actor.GetMapper().GetInput(), "unfolded/model")
            return actor

    def createDedicatedMeshes(self, hierarchy):
        '''
        For each loaded structure in the given hierarchical mesh create a projection mesh depending on the chosen projection method for this structure.
        :param hierarchy: a hierarchical mesh object
        :return:
        '''
        meshes = hierarchy.meshes
        for a in range(len(meshes)):

            if meshes[a].projectionMethod == Mesh.ProjectionMethod.Inflate:
                mesh = hierarchy.unfoldedActor.GetMapper().GetInput()
                mesh = util.smoothMesh(mesh,meshes[a].mesh,iterations=15,relaxation=0.1)
                poly = mesh

            elif meshes[a].projectionMethod == Mesh.ProjectionMethod.Cube:
                poly = util.projectMeshToBounds(hierarchy.unfoldedActor.GetMapper().GetInput())

            elif meshes[a].projectionMethod == Mesh.ProjectionMethod.Clipping:
                poly = hierarchy.unfoldedActor.GetMapper().GetInput()

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputData(poly)

            actor = vtk.vtkActor()
            actor.SetMapper(mapper)

            meshes[a].projectionActor = actor

    def importUnfoldedMesh(self, name):
        '''
        Imports a previous unfolded .obj mesh.
        :param name: The name of the mesh without file extension.
        :return:
        '''
        filename = os.path.join(self.dirname, "../out/3D/unfolded/"+name+".obj")
        mesh = util.readObj(filename)

        mesh = self.normalizeUV(mesh)
        mesh = self.calcMeshNormals(mesh)

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(mesh)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        # actor.GetProperty().SetColor([1.0,1.0,1.0])
        actor.GetProperty().BackfaceCullingOn()
        actor.GetProperty().SetOpacity(0.5)

        return actor

    def normalizeUV(self,mesh):
        '''
        Normalizes the UVs of the given mesh while keeping the ratio between U and V coordinates.
        :param mesh:
        :return: mesh with normalized uvs.
        '''

        mesh = self.shiftUVsToOrigin(mesh)

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
        '''
        Helper method for normalizeUVs and shiftUVsToOrigin.
        :param mesh:
        :param axis:
        :return: Depending on axis, either the min and max of the UVs for each axis,
        or the global min and max.
        '''
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

                umax = np.amax(np.array(uvs)[:,0])
                umin = np.amin(np.array(uvs)[:,0])
                vmax = np.amax(np.array(uvs)[:,1])
                vmin = np.amin(np.array(uvs)[:,1])

                if umax > uMax: uMax = umax
                if umin < uMin: uMin = umin
                if vmax > vMax: vMax = vmax
                if vmin < vMin: vMin = vmin

            return uMin,uMax,vMin,vMax

        else:
            gmax = 0.0
            gmin = 0.0
            for i in range(int(textureCoordinates.GetNumberOfTuples() / 3)):
                uvs = [[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]
                textureCoordinates.GetTuple(mesh.GetCell(i).GetPointId(0), uvs[0])
                textureCoordinates.GetTuple(mesh.GetCell(i).GetPointId(1), uvs[1])
                textureCoordinates.GetTuple(mesh.GetCell(i).GetPointId(2), uvs[2])

                max = np.amax(np.array(uvs))
                min = np.amin(np.array(uvs))
                if max > gmax: gmax = max
                if min < gmin: gmin = min

            return gmin, gmax

    def shiftUVsToOrigin(self,mesh):
        '''
        Shifts the UVs, so all coordinates are positive.
        Necessary before normalize UVs to get constant results.
        :param mesh:
        :return:
        '''
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
        '''
        Filles in the normals of a vtk polydata.
        :param polydata: 
        :return: 
        '''''
        normals = vtk.vtkPolyDataNormals()
        normals.ComputePointNormalsOff()
        normals.ComputeCellNormalsOn()
        normals.ConsistencyOn()
        normals.SetInputData(polydata)
        normals.Update()
        return normals.GetOutput()

    def cutMeshWithPlanes(self,mesh,cutPlanes):
        '''
        Cuts a Mesh with the given planes.
        :param mesh:
        :param cutPlanes:
        :return: A list containing the new mesh segments.
        '''
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
        #todo change to hierachy
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
        '''
        Debug method to unfold a mesh int het out/3D/ folder given by its name.
        :param name:
        :return:
        '''
        inpath = os.path.join(self.dirname, "../out/3D/" + name + ".stl")
        outpath = os.path.join(self.dirname, "../out/3D/papermesh.off")
        util.meshioIO(inpath,outpath)

        graph = Graph()
        graph.load(outpath)
        if not graph.unfold(50000, 0):
            print("failed to unfold :(")
        else:
            print("succesfully unfolded :)")
            filename = os.path.join(self.dirname, "../out/3D/unfolded/" + name + ".obj")
            gluetabs_filename = os.path.join(self.dirname, "../out/3D/unfolded/gluetabs_" + name + ".obj")

            graph.save(filename, gluetabs_filename)

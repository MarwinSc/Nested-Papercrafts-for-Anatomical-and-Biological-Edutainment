import vtkmodules.all as vtk
import numpy as np
import os
from projectionStructure import ProjectionStructure
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

    def mu3dUnfoldPaperMesh(self, mesh, graph, iterations):
        '''
        Forwards the mesh to the mu3d wrapper to unfold it.
        :param actor: The vtk actor containing the mesh to unfold.
        :param graph: The wrapped mu3d graph object.
        :param iterations: The iterations for the unfolding.
        :return: If the unfolding is successful the vtk actor containing the unfolded mesh is returned.
        '''
        inpath = os.path.join(self.dirname, "../out/3D/papermesh.stl")
        outpath = os.path.join(self.dirname, "../out/3D/papermesh.off")
        util.writeStl(mesh,"papermesh")

        util.meshioIO(inpath,outpath)

        graph.load(outpath)
        while not graph.unfold(iterations, 0):
            #msgBox = QtWidgets.QMessageBox()
            #msgBox.setText("failed to unfold :( in {} iterations".format(iterations))
            #msgBox.exec()
            print("failed to unfold :( in {} iterations".format(iterations))
        else:
            print("succesfully unfolded :) in {} iterations".format(iterations))

            filename = os.path.join(self.dirname, "../out/3D/unfolded/model.obj")
            gluetabs_filename = os.path.join(self.dirname, "../out/3D/unfolded/gluetabs.obj")

            graph.save(filename, gluetabs_filename)

            filename = os.path.join(self.dirname, "../out/3D/unfolded/model2D.obj")
            gluetabs_filename = os.path.join(self.dirname, "../out/3D/unfolded/gluetabs.obj")
            gluetabs_mirrored_filename = os.path.join(self.dirname, "../out/3D/unfolded/gluetabs_mirrored.obj")

            graph.save_all(filename, gluetabs_filename,gluetabs_mirrored_filename)

            filename = os.path.join(self.dirname, "../out/3D/unfolded/model.obj")
            mesh = util.readObj(filename)
            tCoords = mesh.GetPointData().GetTCoords()

            #mesh = self.normalizeUV(mesh)
            mesh = util.calcMeshNormals(mesh)

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputData(mesh)
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)

            #actor.GetProperty().SetColor([1.0, 1.0, 1.0])
            actor.GetProperty().BackfaceCullingOn()
            actor.GetProperty().SetOpacity(0.5)


            return actor


    def createDedicatedMesh(self,mesh,unfolded):
        '''
        Creates a dedicated paper mesh to project a given mesh onto it
        :param mesh: the anatomical structure
        :param unfolded: the unfolded papermesh for that structure
        :return:
        '''
        if mesh.projectionMethod == ProjectionStructure.ProjectionMethod.Inflate:
            mesh = util.smoothMesh(unfolded.GetMapper().GetInput(), mesh.mesh , iterations=15, relaxation=0.1)
            poly = mesh

        elif mesh.projectionMethod == ProjectionStructure.ProjectionMethod.Cube:
            poly = util.projectMeshToBounds(unfolded.GetMapper().GetInput())
            util.writeStl(poly,"cubeDedicatedMesh")

        elif mesh.projectionMethod == ProjectionStructure.ProjectionMethod.Clipping:
            mesh = util.cleanMesh(unfolded.GetMapper().GetInput())
            poly = mesh

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(poly)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        return actor


    def importUnfoldedMesh(self, name):
        '''
        Imports a previous unfolded .obj mesh.
        :param name: The name of the mesh without file extension.
        :return:
        '''
        filename = os.path.join(self.dirname, "../out/3D/unfolded/"+name+".obj")
        mesh = util.readObj(filename)

        mesh = self.normalizeUV(mesh)
        mesh = util.calcMeshNormals(mesh)

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


    def cutMeshWithPlanes(self,mesh,cutPlanes, centerPoint):
        '''
        Cuts a Mesh with the given planes.
        If nor planes are given default horizontal plane through the centerPoint is used.
        :param mesh:
        :param cutPlanes:
        :return: A list containing the new mesh segments.
        '''

        if not cutPlanes:
            cutPlanes = []

            # for the "upper" part of the mesh
            plane = vtk.vtkPlane()
            plane.SetOrigin(centerPoint)
            plane.SetNormal(0.0, 0.0, 1.0)
            cutPlanes.append(plane)

            # for the "lower" part of the mesh
            plane = vtk.vtkPlane()
            plane.SetOrigin(centerPoint)
            plane.SetNormal(0.0, 0.0, -1.0)
            cutPlanes.append(plane)

        results = []

        for plane in cutPlanes:

            clip = vtk.vtkClipPolyData()
            clip.GenerateClipScalarsOn()
            clip.GenerateClippedOutputOn()
            clip.SetInputData(mesh)

            clip.SetClipFunction(plane)
            clip.SetValue(0)
            clip.Update()

            #fillHoles = vtk.vtkFillHolesFilter()
            #fillHoles.SetInputData(clip.GetOutput())
            #fillHoles.SetHoleSize(10000.0)
            #fillHoles.Update()

            #filledMesh = fillHoles.GetOutput()
            #filledMesh = self.calcMeshNormals(filledMesh)
            #filledMesh.GetPoints().GetData().Modified()

            results.append(clip.GetOutput())

        return results

    def booleanCGAL(self,mesh,cutout):
        '''
        Subtracts cutout from mesh with cgal boolean and remeshes the cutsurface.
        :param mesh: filepath of the .off file, from the mesh that should be modified
        :param cutout: filepath of the .off file, from the mesh that should be subtracted
        :return:
        '''

        bool = boolean_interface.Boolean_Interface()
        bool.boolean(mesh, cutout)
        '''
        graph = Graph()
        filename = os.path.join(self.dirname, "difference.off")

        graph.load(filename)
        #todo use iterations from ui
        if not graph.unfold(10000, 0):
            print("failed to unfold :(")
        else:
            print("succesfully unfolded :)")
            filename = os.path.join(self.dirname, "../out/3D/unfolded/difference.obj")
            gluetabs_filename = os.path.join(self.dirname, "../out/3D/unfolded/gluetabs_difference.obj")
            graph.save(filename, gluetabs_filename)
        '''
        #set as dedicated meshes for unfolding
        #self.dedicatedPaperMeshes = [upper,lower]

    def unfoldTest(self, name = "upper"):
        '''
        Debug method to unfold a mesh int het out/3D/ folder given by its name.
        :param name:
        :return:
        '''
        outpath = os.path.join(self.dirname, "../out/3D/papermesh.off")

        graph = Graph()
        graph.load(outpath)
        if not graph.unfold(50000, 0):
            print("failed to unfold :(")
        else:
            print("succesfully unfolded :)")
            filename = os.path.join(self.dirname, "../out/3D/unfolded/" + name + ".obj")
            gluetabs_filename = os.path.join(self.dirname, "../out/3D/unfolded/gluetabs_" + name + ".obj")

            graph.save(filename, gluetabs_filename)


    def writeObjMeshes(self,hm):
        '''
        bools and cuts input meshes
        writes the meshes pre and post boolean to the disc
        :return:
        '''
        for child in hm.children:

            trans = vtk.vtkTransform()
            trans.RotateX(-90.)
            transform = vtk.vtkTransformPolyDataFilter()
            transform.SetTransform(trans)

            name = "objTestCuts/originalPapermesh"
            util.writeObj(util.cleanMesh(child.papermesh),name)
            #-again rotated
            name = "objTestCuts/originalPapermesh_rotated"
            transform.SetInputData(util.cleanMesh(child.papermesh))
            transform.Update()
            util.writeObj(transform.GetOutput(),name)

            name = "objTestCuts/originalPapermeshChild"
            util.writeObj(util.cleanMesh(child.children[0].papermesh), name)
            #-shrink mesh
            name = "objTestCuts/originalPapermeshChild_shrunk"
            cutout = util.offsetMesh(child.children[0].papermesh,-2.5)
            util.writeObj(cutout, name)

            #-again rotated
            name = "objTestCuts/originalPapermeshChild_rotated"
            transform.SetInputData(util.cleanMesh(child.children[0].papermesh))
            transform.Update()
            util.writeObj(transform.GetOutput(),name)
            name = "objTestCuts/originalPapermeshChild_shrunk_rotated"
            cutout = util.offsetMesh(transform.GetOutput() ,-2.5)
            util.writeObj(cutout, name)

            centerPoint = child.children[0].papermesh.GetCenter()

            cutPlanes = []

            # for the "upper" part of the mesh
            plane = vtk.vtkPlane()
            plane.SetOrigin(centerPoint)
            plane.SetNormal(0.0, 0.0, 1.0)
            cutPlanes.append(plane)

            # for the "lower" part of the mesh
            plane = vtk.vtkPlane()
            plane.SetOrigin(centerPoint)
            plane.SetNormal(-0.0, -0.0, -1.0)
            cutPlanes.append(plane)

            pieces = self.cutMeshWithPlanes(child.papermesh,cutPlanes,centerPoint)

            for i,piece in enumerate(pieces):
                name = "papermeshPiece{}".format(i)
                file = "objTestCuts/PapermeshCuts/"+name
                piece = util.cleanMesh(piece)
                util.writeObj(piece, file)

                filepath = os.path.join(self.dirname, "../out/3D/" + file + ".obj")
                meshPath = os.path.join(self.dirname, "../out/3D/" + name + ".off")
                util.meshioIO(filepath,meshPath)

                util.writeStl(child.children[0].papermesh, "tempCutout")
                inPath = os.path.join(self.dirname, "../out/3D/tempCutout.stl")
                outPath = os.path.join(self.dirname, "../out/3D/tempCutout.off")
                util.meshioIO(inPath, outPath)

                print(meshPath,outPath)
                self.booleanCGAL(meshPath,outPath)
                inPath = os.path.join(self.dirname,"difference.off")
                outPath = os.path.join(self.dirname,"../out/3D/objTestCuts/BooledPapermeshCuts/booledPapermeshPiece{}".format(i)+".obj")
                util.meshioIO(inPath,outPath)

                #-------- again for rotated meshes

                file = file + "_rotated"
                transform.SetInputData(piece)
                transform.Update()
                util.writeObj(transform.GetOutput(), file)

                filepath = os.path.join(self.dirname, "../out/3D/" + file + ".obj")
                meshPath = os.path.join(self.dirname, "../out/3D/" + name + "_rotated.off")
                util.meshioIO(filepath,meshPath)

                transform.SetInputData(child.children[0].papermesh)
                transform.Update()
                util.writeStl(transform.GetOutput(), "tempCutout_rotated")
                inPath = os.path.join(self.dirname, "../out/3D/tempCutout_rotated.stl")
                outPath = os.path.join(self.dirname, "../out/3D/tempCutout_rotated.off")
                util.meshioIO(inPath, outPath)

                self.booleanCGAL(meshPath,outPath)
                inPath = os.path.join(self.dirname,"difference.off")
                outPath = os.path.join(self.dirname,"../out/3D/objTestCuts/BooledPapermeshCuts/booledPapermeshPiece{}".format(i)+"_rotated.obj")
                util.meshioIO(inPath,outPath)





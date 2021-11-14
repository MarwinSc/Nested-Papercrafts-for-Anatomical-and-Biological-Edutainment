import os
import trimesh
from mim.mim import meshB_inside_meshA
import vtkmodules.all as vtk
import util
import projector
from imageProcessing import ImageProcessor
from mu3d.mu3dpy.mu3d import Graph
from boolean import boolean_interface
from PyQt5 import QtWidgets


class HierarchicalMesh(object):
    """
    Represents a hierarchy of meshes. Tree structure.
    Children linked to parents and parents to children.
    Each HierarchicalMesh is associated with one or multiple Mesh instances
    and a Papermesh.
    """
    dirname = os.path.dirname(__file__)

    def __init__(self, parent, meshes, meshProcessor, faces=50):
        """
        Initialises a hierarchical mesh.
        :param self: this
        :param parent: parent hierarchical mesh.
        :param mesh: mesh or list of meshes associated with this level
        :return: None
        """

        self.children = []
        self.parent = parent
        self.meshProcessor = meshProcessor
        self.imageProcessor = ImageProcessor()

        # new sub node of the tree
        if meshes:
            if isinstance(meshes, list):
                self.meshes = meshes
            else:
                self.meshes = [meshes]

            # generates a papermesh for the loaded structures
            # and writes it to the disk as papermeshLevelTemp.stl and .off
            self.mesh = self.generatePaperMesh(x=faces)
            self.setName(self.writePapermeshStlAndOff("Temp"))

        # no structures thus anchor of the tree
        else:
            self.meshes = []
            self.mesh = None
            self.name = 'Anchor'
            self.file = None
            self.offname = None

        # list holding the unfolded pieces of the papermesh once cut and unfold
        self.unfoldedActors = []
        self.gluetabs = []
        self.modelFiles = []
        self.gluetabFiles = []
        self.mirrorgtFiles = []

    def setName(self, filename):
        '''
        Helper to set all name variables.
        :param filename:
        :return:
        '''
        self.name = filename[filename.rfind('/') + 1:]
        self.offname = filename[:filename.rfind('.')] + ".off"
        self.file = filename

    def render(self, level, renderer):
        """
        Adds actors for the papermeshes to the given renderer based on the level.
        :param level:
        :param renderer:
        :return:
        """
        if level == 0:
            if self.mesh:
                # renderer.AddActor(self.mesh)
                self.mesh.SetVisibility(True)

            for child in self.children:
                child.render(level, renderer)
        else:
            if self.mesh:
                self.mesh.SetVisibility(False)

            for child in self.children:
                child.render(level - 1, renderer)

    def inside(self, mesh):
        """
        Checks if the given mesh is inside this mesh.
        :param mesh: Hierarchical mesh that is checked
        :return: True if the given mesh is inside this mesh
        """

        print("checking if ", mesh.name, " is inside of ", self.name)
        return meshB_inside_meshA(bytes(self.offname, 'utf-8'), bytes(mesh.offname, 'utf-8'))


    def add(self, mesh):
        """
        Adds the given mesh to the hierarchy
        :param mesh: Hierarchical Mesh that should be added to the hierarchy
        :return: None
        """
        # is the mesh inside any of the children?
        # if so add it to the first child we encounter
        for child in self.children:
            if child.add(mesh):
                return True

        # if this object does not have a mesh
        # it's top level therefore we can add any mesh here
        # if it was not added to any children
        # or we are nto top level and it's inside the current mesh
        if self.mesh is None or self.inside(mesh):
            # this means any of this children is potentially a child of the new mesh
            temp = self.children.copy()
            self.children.clear()
            self.appendChild(mesh)
            for tmp_child in temp:
                # write the temp child as .stl and .off and change the name temporarily so it wont be equal to mesh.name
                tmp_child.setName(tmp_child.writePapermeshStlAndOff('Temp'))
                # either add as child to the new hierarchical mesh, or add to the same level as the new hm.
                if not mesh.add(tmp_child):
                    self.appendChild(tmp_child)
            return True

        return False

    def reName(self):
        '''
        Name the mesh according to its current position in the tree
        and also write .stl and .off file with the new name
        also calls reName() of all children.
        :return:
        '''
        for child in self.children:
            idx = "{}_{}".format(child.getLevel(), child.getChildIdx())
            child.setName(child.writePapermeshStlAndOff(idx))
            child.reName()

    def appendChild(self, newChild):
        '''
        Helper method to trigger side effects necessary when appending a child.
        :param newChild:
        :return:
        '''
        self.children.append(newChild)
        newChild.parent = self
        self.reName()

    def recursive_difference(self, convexHull_cutout = False):
        """
        "Cuts" out children of this mesh from this mesh.
        Recursively "cuts" out children of children of children of children ...
        Calls the cutsurface-triangulation for each mesh piece,
        and updates the amount of Buttons in the GUI.
        :return: None
        """
        if self.mesh is not None and self.children:

            try:
                self.boolThenCut(convexHull_cutout)
                #self.cutThenBool()
            except Exception:
                raise Exception("Failure on Cut")

            self.clearUiButtons()
            for piece in self.meshPieces:
                self.updateUiButtons(piece)

        # if all children are cut out save it and call boolean for children
        for child in self.children:
            child.recursive_difference(convexHull_cutout)

    def boolThenCut(self, convexCutout = False):
        '''
        First calculates the difference with the inner mesh, then cuts the given papermesh open and uses boolean.cpp methods to mesh the resulting two boundary loops.
        '''

        final_mesh = trimesh.load(self.file)

        if convexCutout:
            #hull = vtk.vtkHull()
            #hull.SetInputData(self.children[0].papermesh)
            #hull.AddCubeFacePlanes()
            #hull.AddRecursiveSpherePlanes(2)
            #hull.Update()
            delauney = vtk.vtkDelaunay3D()
            delauney.SetInputData(self.children[0].papermesh)
            delauney.Update()
            surfaceFilter = vtk.vtkDataSetSurfaceFilter()
            surfaceFilter.SetInputData(delauney.GetOutput())
            surfaceFilter.Update()

            util.writeStl(surfaceFilter.GetOutput(),"hull")

            tri_child = trimesh.load(os.path.join(self.dirname, "../out/3D/hull.stl"))
        else:
            tri_child = trimesh.load(self.children[0].file)
        # final_mesh = final_mesh.difference(tri_child, engine='blender')
        #final_mesh = trimesh.boolean.difference([final_mesh, tri_child], engine="blender")
        final_mesh = trimesh.boolean.difference([final_mesh, tri_child], engine="scad")

        filename = os.path.join(self.dirname, "../out/3D/differenced_" + self.name)
        final_mesh.export(filename)

        final_mesh_vtk = util.readStl(filename)

        # ----cut planes
        # centerPoint = self.children[0].papermesh.GetCenter()
        cutNormal = self.cutPlane.GetNormal()
        invertedNormal = (-1.0 * float(cutNormal[0]), -1.0 * float(cutNormal[1]), -1.0 * float(cutNormal[2]))
        invPlane = vtk.vtkPlane()
        invPlane.SetOrigin(self.cutPlane.GetOrigin())
        invPlane.SetNormal(invertedNormal[0], invertedNormal[1], invertedNormal[2])
        planes = [self.cutPlane, invPlane]
        print(self.cutPlane.GetOrigin())
        # ----------

        meshPieces = self.meshProcessor.cutMeshWithPlanes(final_mesh_vtk, planes, None)
        paths = []
        paths2 = []
        minDistance = []
        avgDistance = []

        for i in range(len(meshPieces)):
            # write stl and convert to off
            print("name " + self.name)
            util.writeStl(meshPieces[i], "piece{}_".format(i) + self.name[:self.name.rfind('.')])
            inPath = os.path.join(self.dirname, "../out/3D/piece{}_".format(i) + self.name)
            outPath = os.path.join(self.dirname,
                                   "../out/3D/piece{}_".format(i) + self.name[:self.name.rfind('.')] + ".off")
            util.meshioIO(inPath, outPath)
            paths.append(outPath)
            paths2.append(inPath)

            min, max, avg = self.meshProcessor.getAverageEdgeLength(meshPieces[i])
            minDistance.append(min)
            avgDistance.append(avg)

        self.meshPieces = []

        for i, path in enumerate(paths):
            normal = planes[i].GetNormal()
            origin = planes[i].GetOrigin()
            #TODO should be the same value for the whole Hierarchy, but maybe relativ to the first node
            threshold = 30.0 #avgDistance[i]/100.0
            boolean_interface.Boolean_Interface().triangulateCut(paths[i], threshold, -1.0 * float(normal[0]),
                                                                 -1.0 * float(normal[1]), -1.0 * float(normal[2]),
                                                                 float(origin[0]), float(origin[1]), float(origin[2]))
            inPath = os.path.join(self.dirname, "../out/3D/connected.off")
            outPath = paths2[i]
            util.meshioIO(inPath, outPath)

            #bool other children out of each piece
            #TODO this has to be changed if we want to support multiple cuts per parent
            for child in self.children[1:]:
                final_mesh = trimesh.load(outPath)
                tri_child = trimesh.load(child.file)
                # final_mesh = final_mesh.difference(tri_child, engine='blender')
                # final_mesh = trimesh.boolean.difference([final_mesh, tri_child], engine="blender")
                final_mesh = trimesh.boolean.difference([final_mesh, tri_child], engine="scad")
                final_mesh.export(outPath)

            mesh = util.readStl(outPath)

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputData(mesh)
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetOpacity(0.15)
            self.meshPieces.append(actor)

    def cutThenBool(self):
        '''
        First cuts the papermesh open, then uses vtk methods to mesh the resulting single boundary loop, then calculates the difference between each of the pieces and the inner mesh.
        '''

        # ----cut planes
        # centerPoint = self.children[0].papermesh.GetCenter()
        cutNormal = self.cutPlane.GetNormal()
        invertedNormal = (-1.0 * float(cutNormal[0]), -1.0 * float(cutNormal[1]), -1.0 * float(cutNormal[2]))
        invPlane = vtk.vtkPlane()
        invPlane.SetOrigin(self.cutPlane.GetOrigin())
        invPlane.SetNormal(invertedNormal[0], invertedNormal[1], invertedNormal[2])
        planes = [self.cutPlane, invPlane]
        print(self.cutPlane.GetOrigin())
        # ----------

        meshPieces = self.meshProcessor.cutMeshWithPlanes(self.papermesh, planes, None)
        self.meshPieces = []

        for i in range(len(meshPieces)):

            fillHole = vtk.vtkFillHolesFilter()
            fillHole.SetInputData(meshPieces[i])
            fillHole.SetHoleSize(1000.0)
            fillHole.Update()
            normals = vtk.vtkPolyDataNormals()
            normals.SetInputData(fillHole.GetOutput())
            normals.ComputeCellNormalsOn()
            normals.Update()
            mesh = normals.GetOutput()
            mesh = util.cleanMesh(mesh)

            util.writeStl(mesh, "piece{}_".format(i) + self.name[:self.name.rfind('.')])
            path = os.path.join(self.dirname, "../out/3D/piece{}_".format(i) + self.name)

            # --------boolean difference
            final_mesh = trimesh.load(path)
            for child in self.children:
                tri_child = trimesh.load(child.file)
                # final_mesh = final_mesh.difference(tri_child, engine='blender')
                final_mesh = trimesh.boolean.difference([final_mesh, tri_child], engine="scad")

            filename = os.path.join(self.dirname, "../out/3D/differenced_{}_".format(i) + self.name)
            final_mesh.export(filename)

            normal = planes[i].GetNormal()
            origin = planes[i].GetOrigin()
            outPath = os.path.join(self.dirname,
                                   "../out/3D/differenced_{}_".format(i) + self.name[:self.name.rfind('.')] + ".off")
            util.meshioIO(filename, outPath)
            boolean_interface.Boolean_Interface().merge_adjacent_vertices_by_distance(outPath, 1.0,
                                                                                      -1.0 * float(normal[0]),
                                                                                      -1.0 * float(normal[1]),
                                                                                      -1.0 * float(normal[2]),
                                                                                      float(origin[0]),
                                                                                      float(origin[1]),
                                                                                      float(origin[2]))
            inPath = os.path.join(self.dirname, "../out/3D/merged.off")
            outPath = os.path.join(self.dirname, "../out/3D/merged.stl")
            util.meshioIO(inPath, outPath)

            final_mesh_vtk = util.cleanMesh(util.readStl(outPath))
            # util.writeStl(final_mesh_vtk, "debug_differenced_{}_".format(i) + self.name[:self.name.rfind('.')])

            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputData(final_mesh_vtk)
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetOpacity(0.15)
            self.meshPieces.append(actor)

    def renderStructures(self, renderer):
        '''
        Renders all loaded structures in the hierarchy to the given renderer.
        :param renderer:
        :return:
        '''
        for m in self.meshes:
            renderer.AddActor(m.getActor())
        for child in self.children:
            child.renderStructures(renderer)

    def renderPaperMeshes(self, renderer):
        '''
        Renders all papermeshes in the hierarchy to the given renderer.
        :param renderer:
        :return:
        '''
        if self.unfoldedActors:
            for actor in self.unfoldedActors:
                renderer.AddActor(actor)

        elif hasattr(self, "meshPieces"):
            for actor in self.meshPieces:
                renderer.AddActor(actor)

        elif hasattr(self, "unfoldedActor"):
            renderer.AddActor(self.unfoldedActor)

        else:
            renderer.AddActor(self.mesh)

        for child in self.children:
            child.renderPaperMeshes(renderer)


    def generatePaperMesh(self,x=50):

        meshes = [m.mesh for m in self.meshes]
        polysAppended = util.appendMeshes(meshes)

        hull = vtk.vtkHull()
        hull.SetInputData(polysAppended)
        hull.AddCubeFacePlanes()
        hull.Update()
        triangleFilter = vtk.vtkTriangleFilter()
        triangleFilter.SetInputData(hull.GetOutput())
        triangleFilter.Update()

        mesh = util.subdivideMesh(triangleFilter.GetOutput(),4)
        mesh = util.shrinkWrap(mesh, polysAppended)

        util.writeStl(mesh,"pre_approximate")
        inPath = os.path.join(self.dirname, "../out/3D/pre_approximate.stl")
        outPath = os.path.join(self.dirname, "../out/3D/pre_approximate.off")
        util.meshioIO(inPath,outPath)

        threshold = 0.00032*x-0.00039
        # 0.017 -> about 52, 0.012 -> ca 40 faces in the resulting mesh
        boolean_interface.Boolean_Interface().simplify(outPath,threshold)

        outPath = os.path.join(self.dirname, "../out/3D/simplified.stl")
        inPath = os.path.join(self.dirname, "../out/3D/simplified.off")
        util.meshioIO(inPath,outPath)

        self.papermesh = util.cleanMesh(util.readStl(outPath))

        '''
        com = vtk.vtkCenterOfMass()
        com.SetInputData(mesh)
        com.Update()
        com = com.GetCenter()

        transform = vtk.vtkTransform()
        transform.Translate(com[0],com[1],com[2])
        transform.Scale(1.1, 1.1, 1.1)
        transform.Translate(-com[0],-com[1],-com[2])
        transform.Update()

        scale = vtk.vtkTransformPolyDataFilter()
        scale.SetTransform(transform)
        scale.SetInputData(mesh)
        scale.Update()
        self.papermesh = scale.GetOutput()

        '''

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(self.papermesh)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetOpacity(0.15)
        return actor


    '''
    def generatePaperMesh(self, convex = False, adaptive = False):
        
        #Generates a papermesh for the loaded structures in self.meshes.
        #As a side effect the papermesh itself is saved to self.papermesh.
        #:return: A vtk actor of the generated papermesh.
        
        meshes = [m.mesh for m in self.meshes]
        polysAppended = util.appendMeshes(meshes)
        hull = vtk.vtkHull()
        hull.SetInputData(polysAppended)
        hull.AddCubeFacePlanes()
        hull.Update()
        triangleFilter = vtk.vtkTriangleFilter()
        triangleFilter.SetInputData(hull.GetOutput())
        triangleFilter.Update()

        mesh = util.subdivideMesh(triangleFilter.GetOutput())
        mesh = util.shrinkWrap(mesh, polysAppended)

        if adaptive:
            #TODO if adaptive is used the edge length has to be based on the average edge length of the respective mesh
            mesh = util.adaptiveSubdivideMesh(mesh, maxEdgeLength=30.0, numberOfTriangles=100, numberOfPasses=10)
            mesh = util.shrinkWrap(mesh, polysAppended)
            mesh = util.adaptiveSubdivideMesh(mesh, maxEdgeLength=15.0, numberOfTriangles=100, numberOfPasses=10)
            mesh = util.shrinkWrap(mesh, polysAppended)


        mesh = util.cleanMesh(mesh)

        if convex:
            delauney = vtk.vtkDelaunay3D()
            delauney.SetInputData(mesh)
            delauney.Update()
            surfaceFilter = vtk.vtkDataSetSurfaceFilter()
            surfaceFilter.SetInputData(delauney.GetOutput())
            surfaceFilter.Update()
            mesh = surfaceFilter.GetOutput()

        self.papermesh = mesh

        com = vtk.vtkCenterOfMass()
        com.SetInputData(mesh)
        com.Update()
        com = com.GetCenter()

        transform = vtk.vtkTransform()
        transform.Translate(com[0],com[1],com[2])
        transform.Scale(1.1, 1.1, 1.1)
        transform.Translate(-com[0],-com[1],-com[2])
        transform.Update()

        scale = vtk.vtkTransformPolyDataFilter()
        scale.SetTransform(transform)
        scale.SetInputData(mesh)
        scale.Update()
        self.papermesh = scale.GetOutput()

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(self.papermesh)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetOpacity(0.15)
        return actor
    '''

    def unfoldWholeHierarchy(self, iterations):
        '''
        Calls the respective unfolding methodology, depending on if the papermeshes were cut.
        '''
        if hasattr(self, 'meshPieces'):
            self.unfoldPaperMeshPieces(iterations)
        elif hasattr(self, 'papermesh'):
            self.unfoldPaperMesh(iterations)
        for child in self.children:
            child.unfoldWholeHierarchy(iterations)

    def unfoldPaperMesh(self, iterations):
        '''
        Calls the mu3dUnfoldPaperMesh() method of the given meshProcessor and afterward the createDedicatedPaperMesh()
        to unfold self.papermesh.
        :param meshProcessor:
        :param iterations:
        :return:
        '''
        # idx = "{}_{}".format(self.getLevel(), self.getChildIdx())
        # self.writePapermeshStlAndOff(idx)
        self.graph = Graph()

        unfoldedActor = self.meshProcessor.mu3dUnfoldPaperMesh(self.getChildIdx(), self.papermesh, self.graph, iterations)
        if unfoldedActor:
            self.label.setText("Unfolded")
            self.unfoldedActor = unfoldedActor

            # import labels
            # for some reason crashes if imported as .obj
            filename = os.path.join(self.dirname, "../out/3D/unfolded/gluetabs.obj")
            outpath = os.path.join(self.dirname, "../out/3D/unfolded/gluetabs.stl")
            util.meshioIO(filename, outpath)
            self.gluetab = util.readStl(outpath)
            self.modelFile = os.path.join(self.dirname, "../out/3D/unfolded/model{}.obj".format(self.getChildIdx()))
            self.gluetabFile = os.path.join(self.dirname, "../out/3D/unfolded/gluetabs{}.obj".format(self.getChildIdx()))
            self.mirrorgtFile = os.path.join(self.dirname, "../out/3D/unfolded/gluetabs_mirrored{}.obj".format(self.getChildIdx()))

    def unfoldPaperMeshPieces(self, iterations):
        '''
        Calls the mu3dUnfoldPaperMesh() method of the given meshProcessor and afterward the createDedicatedPaperMesh()
        to unfold each mesh in self.meshPieces.
        :param meshProcessor:
        :param iterations:
        :return:
        '''
        unfoldedString = ""
        for i, piece in enumerate(self.meshPieces):
            # idx = "{}_{}_Piece{}".format(self.getLevel(), self.getChildIdx(), i)
            # self.writePapermeshStlAndOff(idx)
            self.graph = Graph()
            unfoldedActor = self.meshProcessor.mu3dUnfoldPaperMesh(self.getChildIdx() + "_" + str(i), piece.GetMapper().GetInput(), self.graph, iterations)
            if unfoldedActor:
                unfoldedString += "Unfolded, "
                self.unfoldedActors.append(unfoldedActor)

                # import labels
                # for some reason crashes if imported as .obj
                filename = os.path.join(self.dirname, "../out/3D/unfolded/gluetabs.obj")
                outpath = os.path.join(self.dirname, "../out/3D/unfolded/gluetabs.stl")
                util.meshioIO(filename, outpath)
                self.gluetabs.append(util.readStl(outpath))
                self.modelFiles.append(os.path.join(self.dirname, "../out/3D/unfolded/model{}.obj".format(self.getChildIdx() + "_" + str(i))))
                self.gluetabFiles.append(os.path.join(self.dirname,
                                                "../out/3D/unfolded/gluetabs{}.obj".format(self.getChildIdx() + "_" + str(i))))
                self.mirrorgtFiles.append(os.path.join(self.dirname, "../out/3D/unfolded/gluetabs_mirrored{}.obj".format(
                    self.getChildIdx() + "_" + str(i))))

            else:
                unfoldedString += "Not Unfolded, "
        self.label.setText(unfoldedString)

    def project(self, resolution):
        '''
        Projects all anatomical structures of the hierarchy onto their respective unfolded papermesh.
        Creating the unfolded paper template for each papermesh.
        '''
        if self.parent:
            self.clearUiButtons()
            # creating unfoldings per mesh piece
            if hasattr(self, "meshPieces"):
                for j, unfolded_piece in enumerate(self.unfoldedActors):
                    unfolded_piece = self.project_helper(unfolded_piece, resolution, self.getChildIdx()+"_"+str(j), j)
                    self.updateUiButtons(unfolded_piece)
                    print("Projected a hierarchy level.")

            # create single unfolding if mesh isn't cut
            else:
                self.unfoldedActor = self.project_helper(self.unfoldedActor, resolution, self.getChildIdx())
                self.updateUiButtons(self.unfoldedActor)
                print("Projected a hierarchy level.")

        self.getChildIdx()
        for child in self.children:
            child.project(resolution)

    def project_helper(self, actor, resolution, textureIdx, j=-1):
        '''
        Helper for projectMethod for either pieces of a cut mesh or an uncut papermesh.
        Handels calling the necessary methods for projection with a given actor and the meshes of the current hierarchical level.
        '''
        for i, mesh in enumerate(self.meshes):
            # if there are multiple pieces
            if j < 0:
                dedicatedMesh = self.meshProcessor.createDedicatedMesh(mesh, actor)
                projection = projector.projectPerTriangle(dedicatedMesh, mesh.getActor(), textureIdx + "_" + str(i), resolution)
                # createUnfoldedPaperMesh is deprecated but still used for the previewTexture
                img, previewTexture = projector.createUnfoldedPaperMesh(projection, actor, self.gluetab, textureIdx  + "_" + str(i))
            # if the papermesh is complete
            else:
                map_innerMesh_to_cutPlane = self.create_projection_dictionary_for_cutout()
                dedicatedMesh = self.meshProcessor.createDedicatedMesh(mesh, actor, hull = self.cut_boundingBox(actor.GetMapper().GetInput(),j), map = map_innerMesh_to_cutPlane)
                #dedicatedMesh = self.meshProcessor.createDedicatedMesh(mesh, actor)
                projection = projector.projectPerTriangle(dedicatedMesh, mesh.getActor(), textureIdx + "_" + str(i), resolution)
                # createUnfoldedPaperMesh is deprecated but still used for the previewTexture
                img, previewTexture = projector.createUnfoldedPaperMesh(projection, actor, self.gluetabs[j], textureIdx + "_" + str(i))

            if j < 0:
                front_output, back_output, labelIDs_output = projector.renderFinalOutput(self.modelFile, self.gluetabFile, self.mirrorgtFile, textureIdx + "_" + str(i),
                                            projection.GetMapper().GetInput().GetPointData().GetTCoords())
            else:
                front_output, back_output, labelIDs_output = projector.renderFinalOutput(self.modelFiles[j], self.gluetabFiles[j], self.mirrorgtFiles[j], textureIdx + "_" + str(i),
                                            projection.GetMapper().GetInput().GetPointData().GetTCoords())

            previewTexture = projector.mask(util.VtkToNp(previewTexture))
            front_output = util.VtkToNp(front_output)
            # to multiply textures of successive meshes
            if i > 0:
                previousPreviewTexture = self.multiplyProjections(previewTexture, previousPreviewTexture)
                previousFrontOutput = self.multiplyProjections(front_output,previousFrontOutput)
            # for the first mesh or if we have only one
            else:
                previousPreviewTexture = self.multiplyProjections(previewTexture, None)
                previousFrontOutput = self.multiplyProjections(front_output, None)

            print("Projected a mesh")

        # subtract the labelIDs rendering from the multiplied texture, masking front and back and save to disk
        previousFrontOutput = previousFrontOutput - util.VtkToNp(labelIDs_output)

        previousFrontOutput, back_output = projector.mask(previousFrontOutput, util.VtkToNp(back_output))

        filename = os.path.join(self.dirname, "../out/2D/unfolding_back{}.png".format(textureIdx))
        dy, dx, dz = back_output.shape
        util.writeImage(util.NpToVtk(back_output, dx, dy, dz), filename)
        filename = os.path.join(self.dirname, "../out/2D/unfolding{}.png".format(textureIdx))
        dy, dx, dz = previousFrontOutput.shape
        util.writeImage(util.NpToVtk(previousFrontOutput, dx, dy, dz), filename)

        # Set Texture for preview
        texture = vtk.vtkTexture()
        dy, dx, dz = previousPreviewTexture.shape
        previousPreviewTexture = util.NpToVtk(previousPreviewTexture, dx, dy, dz)
        castFilter = vtk.vtkImageCast()
        castFilter.SetInputData(previousPreviewTexture)
        castFilter.SetOutputScalarTypeToUnsignedChar()
        castFilter.Update()
        texture.SetInputData(castFilter.GetOutput())
        actor.GetMapper().SetInputData(self.meshProcessor.normalizeUV(actor.GetMapper().GetInput()))
        actor.SetTexture(texture)
        actor.Modified()
        actor.GetProperty().SetOpacity(1.0)

        return actor

    def cut_boundingBox(self, mesh, i):
        hull = vtk.vtkHull()
        hull.SetInputData(mesh)
        hull.AddCubeFacePlanes()
        hull.Update()
        cutNormal = self.cutPlane.GetNormal()
        invertedNormal = (-1.0 * float(cutNormal[0]), -1.0 * float(cutNormal[1]), -1.0 * float(cutNormal[2]))
        invPlane = vtk.vtkPlane()
        invPlane.SetOrigin(self.cutPlane.GetOrigin())
        invPlane.SetNormal(invertedNormal[0], invertedNormal[1], invertedNormal[2])
        planes = [self.cutPlane, invPlane]
        hullPieces = self.meshProcessor.cutMeshWithPlanes(hull.GetOutput(), planes, None)
        fillHoles = vtk.vtkFillHolesFilter()
        fillHoles.SetInputData(hullPieces[i])
        fillHoles.SetHoleSize(10000.0)
        fillHoles.Update()
        filledMesh = fillHoles.GetOutput()
        filledMesh = util.calcMeshNormals(filledMesh)
        filledMesh.GetPoints().GetData().Modified()
        return filledMesh

    def create_projection_dictionary_for_cutout(self):
        '''
        creates a dictionary that has the vertex positions of the inner mesh as key, and the projected position onto the cut plane as value.
        '''
        dictionary = {}
        newGeometry = vtk.vtkPolyData()
        newPoints = vtk.vtkPoints()

        for child in self.children:
            if hasattr(child, "meshPieces"):
                for j, unfolded_piece in enumerate(child.unfoldedActors):
                    poly = unfolded_piece.GetMapper().GetInput()
                    points = poly.GetPoints()
                    for i in range(points.GetNumberOfPoints()):
                        projected = [0.0, 0.0, 0.0]
                        self.cutPlane.GeneralizedProjectPoint(points.GetPoint(i), projected)
                        dictionary[points.GetPoint(i)] = projected
                        newPoints.InsertNextPoint(projected)
            else:
                # support also cut inner meshes or does unfolded actor remain unchanged after cutting
                poly = child.unfoldedActor.GetMapper().GetInput()
                points = poly.GetPoints()
                for i in range(points.GetNumberOfPoints()):
                    projected = [0.0,0.0,0.0]
                    self.cutPlane.GeneralizedProjectPoint(points.GetPoint(i),projected)
                    dictionary[points.GetPoint(i)] = projected
                    newPoints.InsertNextPoint(projected)

        newGeometry.SetPoints(newPoints)
        util.writeStl(newGeometry,"debu_projected_points")
        return dictionary

    def multiplyProjections(self,image,prev_image):
        '''
        Multiplies a previously written image and a given image and writes the result to disk.
        :param image: a numpy array, the previously written file
        :param prev_image: the image to multiply with
        :param idx: number to name the created image file.
        '''
        # multiply created unfoldings

        result = self.imageProcessor.optimizedBrighten(image)
        if prev_image is not None:
            result = self.imageProcessor.normalizeMultiplication(prev_image, result)
        return result

    def getAllMeshes(self, asActor=True):
        '''
        :param asActor: if true the meshes are returned as vtk actors, if false as ProjectionStructure-objects.
        :return: All meshes of this level and below as vtk actors or Mesh-objects.
        '''
        if asActor:
            actors = [m.getActor() for m in self.meshes]
        else:
            actors = [m for m in self.meshes]

        for child in self.children:
            actors.extend(child.getAllMeshes())
        return actors

    def writePapermeshStlAndOff(self, levelIdx):
        '''
        Writes the papermesh in .stl format and .off format to the disk.
        named papermesh{level}.
        :return: The full path to the stl file.
        '''
        name = "papermeshLevel{}".format(levelIdx)
        util.writeStl(util.cleanMesh(self.papermesh), name)
        inpath = os.path.join(self.dirname, "../out/3D/papermeshLevel{}.stl".format(levelIdx))
        outpath = os.path.join(self.dirname, "../out/3D/papermeshLevel{}.off".format(levelIdx))
        util.meshioIO(inpath, outpath)
        return inpath

    def writeAllUnfoldedMeshesInHierarchy(self):
        '''
        Writes the unfoldedMeshes of the whole hierarchy to the disk,
        to be able to import them later on.
        '''
        if self.parent:
            if hasattr(self, "meshPieces"):
                for i, unfolded_piece in enumerate(self.unfoldedActors):
                    util.writeObj(unfolded_piece.GetMapper().GetInput(),
                                  "unfolded_{}_piece_{}".format(self.getChildIdx(), i))
                    util.writeObj(self.gluetabs[i], "gluetabs_{}_piece_{}".format(self.getChildIdx(), i))
            else:
                util.writeObj(self.unfoldedActor.GetMapper().GetInput(), "unfolded_{}".format(self.getChildIdx()))
                util.writeObj(self.gluetab, "gluetabs_{}".format(self.getChildIdx()))

        for child in self.children:
            child.writeAllUnfoldedMeshesInHierarchy()

    def importAllUnfoldedMeshesInHierarchy(self):
        '''
        Imports the unfoldedMeshes written in the previous session,
        for debugging purposes.
        '''
        if self.parent:
            directory = os.path.join(self.dirname, r"../out/3D/unfolded")
            unfoldedString = ""
            for filename in os.listdir(directory):
                if hasattr(self, "meshPieces"):
                    for i in range(len(self.meshPieces)):
                        if filename == "model{}.obj".format(self.getChildIdx() + "_" + str(i)):
                            mesh = util.readObj(os.path.join(directory, filename))
                            mapper = vtk.vtkPolyDataMapper()
                            mapper.SetInputData(mesh)
                            actor = vtk.vtkActor()
                            actor.SetMapper(mapper)
                            actor.GetProperty().SetOpacity(0.5)
                            self.unfoldedActors.append(actor)
                            unfoldedString += "Unfolded, "
                            self.modelFiles.append(filename)
                        elif filename == "gluetabs{}.obj".format(self.getChildIdx() + "_" + str(i)):
                            mesh = util.readObj(os.path.join(directory, filename))
                            self.gluetabs.append(mesh)
                            self.gluetabFiles.append(filename)
                        elif filename == "gluetabs_mirrored{}.obj".format(self.getChildIdx() + "_" + str(i)):
                            self.mirrorgtFiles.append(filename)

                else:
                    if filename == "model{}.obj".format(self.getChildIdx()):
                        mesh = util.readObj(os.path.join(directory, filename))
                        mapper = vtk.vtkPolyDataMapper()
                        mapper.SetInputData(mesh)
                        actor = vtk.vtkActor()
                        actor.SetMapper(mapper)
                        actor.GetProperty().SetOpacity(0.5)
                        self.unfoldedActor = actor
                        self.modelFile = filename
                        unfoldedString = "Unfolded"
                    elif filename == "gluetabs{}.obj".format(self.getChildIdx()):
                        mesh = util.readObj(os.path.join(directory, filename))
                        self.gluetab = mesh
                        self.gluetabFile = filename
                    elif filename == "gluetabs_mirrored{}.obj".format(self.getChildIdx()):
                        self.mirrorgtFile = filename

            self.label.setText(unfoldedString)

        for child in self.children:
            child.importAllUnfoldedMeshesInHierarchy()

    def toString(self):
        '''
        Prints the Hierarchy Tree in the console.
        :return:
        '''
        print("level: ", self.getLevel(), "childIdx: ", self.getChildIdx())
        names = [m.filename.split("/")[-1] for m in self.meshes]
        for name in names:
            print(name, end=", ")
        print("")
        if hasattr(self, "cutPlane"):
            print("cutPlane, Normal: {} Origin: {}".format(self.cutPlane.GetNormal(), self.cutPlane.GetOrigin()))
        for child in self.children:
            child.toString()

    def getChildIdx(self):
        '''
        the concatenated childIdx of all parents to this mesh, to be able to assign a unique name encoding the position in the hierarchy to the papermesh.
        e.g.: 100 ... mesh is second child of parent, which is first child of its parent, which is the first child of the anchor mesh.
        :return:
        '''
        if not self.parent:
            return 0
        childIdx = "{}".format(self.parent.children.index(self))
        parent = self.parent
        # if there is another level above the parent
        # sum the childIdx
        while parent is not None:
            if parent.parent:
                childIdx += "{}".format(parent.parent.children.index(parent))
            parent = parent.parent
        return childIdx

    def getLevel(self):
        """
        :return: the level of this hierarchical mesh.
        """
        if not self.parent:
            return 0
        parent = self.parent
        level = 0
        while parent is not None:
            parent = parent.parent
            level += 1
        return level

    def toLevelDict(self, hierarchical_dict):
        """
        :return:
        """
        current_level = self.getLevel()
        print(current_level, self.file)
        if current_level not in hierarchical_dict:
            hierarchical_dict[current_level] = []

        hierarchical_dict[current_level].append(self.file)

        for child in self.children:
            child.toLevelDict(hierarchical_dict)

        return hierarchical_dict

    def toList(self):
        level_dict = self.toLevelDict({})
        del level_dict[0]
        return [v for k, v in level_dict.items()]

    def toActorDict(self, hierarchical_dict):
        """
        :return:
        """
        current_level = self.getLevel()
        if current_level not in hierarchical_dict:
            hierarchical_dict[current_level] = []

        for mesh in self.meshes:
            hierarchical_dict[current_level].append(mesh.getActor())

        for child in self.children:
            child.toActorDict(hierarchical_dict)

        return hierarchical_dict

    def toActorList(self):
        actorDict = self.toActorDict({})
        del actorDict[0]
        list = [v for k, v in actorDict.items()]
        return list

    def addCutPlanes(self, viewpoints, bds):
        '''
        Adds the cuttingplanes to the respective level in the hierarchy.
        '''
        if self.parent:
            plane = vtk.vtkPlane()
            #x0, y0, z0 = (bds[self.getLevel() - 1][1] - abs(bds[self.getLevel() - 1][0])) / 2., (
            #            bds[self.getLevel() - 1][3] - abs(bds[self.getLevel() - 1][2])) / 2., (
            #                         bds[self.getLevel() - 1][5] - abs(bds[self.getLevel() - 1][4])) / 2.
            meshes = []
            for child in self.children:
                meshes.append(child.papermesh)
            if meshes:
                mesh = util.appendMeshes(meshes)
                com = vtk.vtkCenterOfMass()
                com.SetInputData(mesh)
                com.Update()
                x0,y0,z0 = com.GetCenter()

                plane.SetOrigin(x0, y0, z0)
                plane.SetNormal(viewpoints[self.getLevel() - 1][0], viewpoints[self.getLevel() - 1][1],
                                viewpoints[self.getLevel() - 1][2])
                self.cutPlane = plane

        for child in self.children:
            child.addCutPlanes(viewpoints, bds)

    def setUpAdditionalUiElements(self):
        '''
        Adds Buttons to the GUI to switch the visibility of the papermeshes.
        '''
        self.ui_elements_box = QtWidgets.QGroupBox()
        self.ui_elements_layout = QtWidgets.QHBoxLayout()
        self.ui_elements_box.setLayout(self.ui_elements_layout)
        self.updateUiButtons(self.mesh)

    def clearUiButtons(self):
        while self.ui_elements_layout.count():
            child = self.ui_elements_layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()

    def updateUiButtons(self, actor):
        visibilityButton = QtWidgets.QPushButton("Render")
        visibilityButton.setCheckable(True)
        visibilityButton.setChecked(True)

        def onVisibilityChange():
            if (visibilityButton.isChecked()):
                actor.SetVisibility(True)
            else:
                actor.SetVisibility(False)

        visibilityButton.clicked.connect(onVisibilityChange)
        self.ui_elements_layout.addWidget(visibilityButton)

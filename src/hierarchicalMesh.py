import os
import trimesh
from mim.mim import meshB_inside_meshA
import vtkmodules.all as vtk
import util
from mu3d.mu3dpy.mu3d import Graph
from boolean import boolean_interface

class HierarchicalMesh(object):
    """
    Represents a hierarchy of meshes. Tree structure.
    Children linked to parents and parents to children.
    Each HierarchicalMesh is associated with one or multiple Mesh instances
    and a Papermesh.
    """
    dirname = os.path.dirname(__file__)

    def __init__(self, parent, meshes, meshProcessor):
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

        # new sub node of the tree
        if meshes:
            if isinstance(meshes,list):
                self.meshes = meshes
            else:
                self.meshes = [meshes]

            # generates a papermesh for the loaded structures
            # and writes it to the disk as papermeshLevelTemp.stl and .off
            self.mesh = self.generatePaperMesh()
            self.setName(self.writePapermeshStlAndOff("Temp"))

        # no structures thus anchor of the tree
        else:
            self.meshes = []
            self.mesh = None
            self.name = 'Anchor'
            self.file = None
            self.offname = None

    def setName(self,filename):
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
                child.render(level-1, renderer)

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
                #write the temp child as .stl and .off and change the name temporarily so it wont be equal to mesh.name
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

    def appendChild(self,newChild):
        '''
        Helper method to trigger side effects necessary when appending a child.
        :param newChild:
        :return:
        '''
        self.children.append(newChild)
        newChild.parent = self
        self.reName()

    def recursive_difference(self):
        """
        "Cuts" out children of this mesh from this mesh.
        Recursively "cuts" out children of children of children of children ...
        :return: None
        """
        if self.mesh is not None and self.children:
            final_mesh = trimesh.load(self.file)
            for child in self.children:
                tri_child = trimesh.load(child.file)
                #final_mesh = final_mesh.difference(tri_child, engine='blender')
                final_mesh = trimesh.boolean.difference([final_mesh, tri_child], engine="scad")

            filename = os.path.join(self.dirname, "../out/3D/differenced_" + self.name)
            final_mesh.export(filename)

            final_mesh_vtk = util.readStl(filename)

            #----cut planes
            #centerPoint = self.children[0].papermesh.GetCenter()
            cutNormal = self.cutPlane.GetNormal()
            invertedNormal = (-1.0 * float(cutNormal[0]), -1.0 * float(cutNormal[1]), -1.0 * float(cutNormal[2]))
            invPlane = vtk.vtkPlane()
            invPlane.SetOrigin(self.cutPlane.GetOrigin())
            invPlane.SetNormal(invertedNormal[0],invertedNormal[1],invertedNormal[2])
            planes = [self.cutPlane, invPlane]
            print(self.cutPlane.GetOrigin())
            #----------
            '''
            centerPoint = self.children[0].papermesh.GetCenter()
            planes = []

            # for the "upper" part of the mesh
            plane = vtk.vtkPlane()
            plane.SetOrigin(centerPoint)
            plane.SetNormal(0.0, 0.0, 1.0)
            planes.append(plane)

            # for the "lower" part of the mesh
            plane = vtk.vtkPlane()
            plane.SetOrigin(centerPoint)
            plane.SetNormal(0.0, 0.0, -1.0)
            planes.append(plane)
            '''

            meshPieces = self.meshProcessor.cutMeshWithPlanes(final_mesh_vtk, planes, None)
            paths = []
            paths2 = []
            for i in range(len(meshPieces)):
                # write stl and convert to off
                print("name " + self.name)
                util.writeStl(meshPieces[i], "piece{}_".format(i) + self.name[:self.name.rfind('.')])
                inPath = os.path.join(self.dirname, "../out/3D/piece{}_".format(i) + self.name)
                outPath = os.path.join(self.dirname, "../out/3D/piece{}_".format(i) + self.name[:self.name.rfind('.')]+".off")
                util.meshioIO(inPath, outPath)
                paths.append(outPath)
                paths2.append(inPath)

            self.meshPieces = []

            for i,path in enumerate(paths):
                normal = planes[i].GetNormal()
                origin = planes[i].GetOrigin()
                boolean_interface.Boolean_Interface().triangulateCut(paths[i],-1.0 * float(normal[0]),-1.0 * float(normal[1]),-1.0 * float(normal[2]), float(origin[0]), float(origin[1]), float(origin[2]))
                inPath = os.path.join(self.dirname, "../out/3D/connected.off")
                outPath = paths2[i]
                util.meshioIO(inPath, outPath)

                mesh = util.readStl(outPath)

                mapper = vtk.vtkPolyDataMapper()
                mapper.SetInputData(mesh)
                actor = vtk.vtkActor()
                actor.SetMapper(mapper)
                actor.GetProperty().SetOpacity(0.15)
                self.meshPieces.append(actor)

        # if all children are cut out save it and call boolean for children
        for child in self.children:
            child.recursive_difference()

    def renderStructures(self,renderer):
        '''
        Renders all loaded structures in the hierarchy to the given renderer.
        :param renderer:
        :return:
        '''
        for m in self.meshes:
            renderer.AddActor(m.getActor())
        for child in self.children:
            child.renderStructures(renderer)

    def renderPaperMeshes(self,renderer):
        '''
        Renders all papermeshes in the hierarchy to the given renderer.
        :param renderer:
        :return:
        '''
        if hasattr(self,"meshPieces"):
            for actor in self.meshPieces:
                renderer.AddActor(actor)
        else:
            renderer.AddActor(self.mesh)
        for child in self.children:
            child.renderPaperMeshes(renderer)

    def generatePaperMesh(self):
        '''
        Generates a papermesh for the loaded structures in self.meshes.
        As a side effect the papermesh itself is saved to self.papermesh.
        :return: A vtk actor of the generated papermesh.
        '''
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
        mesh = util.cleanMesh(mesh)
        mesh = util.shrinkWrap(mesh,polysAppended)
        self.papermesh = util.offsetMesh(mesh)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(self.papermesh)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetOpacity(0.15)
        return actor

    def unfoldWholeHierarchy(self, iterations):
        if hasattr(self,'meshPieces'):
            self.unfoldPaperMeshPieces(iterations)
        elif hasattr(self,'papermesh'):
            self.unfoldPaperMesh(iterations)
        for child in self.children:
            child.unfoldWholeHierarchy(iterations)

    def unfoldPaperMesh(self, iterations):
        '''
        Calls the mu3dUnfoldPaperMesh() method of the given meshProcessor and afterward the createDedicatedPaperMesh()
        to create a projectionMesh for each mesh in self.meshes.
        :param meshProcessor:
        :param iterations:
        :return:
        '''
        idx = "{}_{}".format(self.getLevel(), self.getChildIdx())
        self.writePapermeshStlAndOff(idx)
        self.graph = Graph()
        unfoldedActor = self.meshProcessor.mu3dUnfoldPaperMesh(self.papermesh, self.graph, iterations)
        if unfoldedActor:
            self.label.setText("Unfolded")
            self.unfoldedActor = unfoldedActor
            #for mesh in self.meshes:
            #    mesh.dedicatedMeshes = [self.meshProcessor.createDedicatedMesh(mesh,unfoldedActor)]

    def unfoldPaperMeshPieces(self, iterations):
        self.unfoldedActors = []
        unfoldedString = ""
        for i,piece in enumerate(self.meshPieces):
            idx = "{}_{}_Piece{}".format(self.getLevel(), self.getChildIdx(), i)
            self.writePapermeshStlAndOff(idx)
            self.graph = Graph()
            unfoldedActor = self.meshProcessor.mu3dUnfoldPaperMesh(piece.GetMapper().GetInput(), self.graph, iterations)
            if unfoldedActor:
                unfoldedString += "Unfolded, "
                self.unfoldedActors.append(unfoldedActor)
            else:
                unfoldedString += "Not Unfolded, "
        self.label.setText(unfoldedString)

        #for mesh in self.meshes:
        #    mesh.dedicatedMeshes = []
        #    for actor in self.unfoldedActors:
        #        mesh.dedicatedMeshes.append(self.meshProcessor.createDedicatedMesh(mesh, actor))

    def getAllMeshes(self, asActor = True):
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
        util.meshioIO(inpath,outpath)
        return inpath


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
        if hasattr(self,"cutPlane"):
            print("cutPlane, Normal: {} Origin: {}".format(self.cutPlane.GetNormal(),self.cutPlane.GetOrigin()))
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

        if self.parent:
            plane = vtk.vtkPlane()
            x0, y0, z0 = (bds[self.getLevel()-1][1] - abs(bds[self.getLevel()-1][0])) / 2., (bds[self.getLevel()-1][3] - abs(bds[self.getLevel()-1][2])) / 2., (bds[self.getLevel()-1][5] - abs(bds[self.getLevel()-1][4])) / 2.
            plane.SetOrigin(x0,y0,z0)
            plane.SetNormal(viewpoints[self.getLevel()-1][0], viewpoints[self.getLevel()-1][1], viewpoints[self.getLevel()-1][2])
            self.cutPlane = plane

        for child in self.children:
            child.addCutPlanes(viewpoints, bds)



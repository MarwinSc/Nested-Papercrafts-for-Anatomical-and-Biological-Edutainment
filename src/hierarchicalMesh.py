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
        if self.mesh is not None:
            final_mesh = trimesh.load(self.file)
            for child in self.children:
                tri_child = trimesh.load(child.file)
                #final_mesh = final_mesh.difference(tri_child, engine='blender')
                final_mesh = trimesh.boolean.difference([final_mesh, tri_child], engine="blender")

            filename = os.path.join(self.dirname, "../out/3D/differenced_" + self.name)
            final_mesh.export(filename)

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
        if hasattr(self,'papermesh'):
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
            self.unfoldedActor = unfoldedActor
            self.meshProcessor.createDedicatedMeshes(self)

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
        util.writeStl(self.papermesh, name)
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
        for child in self.children:
            child.toString()

    def getChildIdx(self):
        '''
        the concatenated childIdx of all parents to this mesh, to be able to assign a unique name  encoding the position in the hierarchy to the papermesh.
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

    def difference(self):

        #currently cuting first child
        centerPoint = self.children[0].papermesh.GetCenter()

        meshPieces = self.meshProcessor.cutMeshWithPlanes(self.papermesh,None,centerPoint)

        util.writeStl(self.children[0].papermesh, "tempCutout")
        inPath = os.path.join(self.dirname, "../out/3D/tempCutout.stl")
        outPath = os.path.join(self.dirname, "../out/3D/tempCutout.off")
        util.meshioIO(inPath, outPath)

        for i in range(len(meshPieces)):
            # write stl and convert to off
            util.writeStl(meshPieces[i],"tempMeshPiece{}".format(i))
            inPath = os.path.join(self.dirname, "../out/3D/tempMeshPiece{}.stl".format(i))
            outPath = os.path.join(self.dirname, "../out/3D/tempMeshPiece{}.off".format(i))
            util.meshioIO(inPath,outPath)

        mesh = os.path.join(self.dirname, "../out/3D/tempMeshPiece1.off")
        cutout = os.path.join(self.dirname, "../out/3D/tempCutout.off")

        self.meshProcessor.booleanCGAL(mesh,cutout)

    def union(self):

        cube = vtk.vtkCubeSource()
        cube.SetXLength(1000.0)
        cube.SetYLength(1000.0)
        cube.SetZLength(0.001)
        cube.SetCenter(self.children[0].papermesh.GetCenter())
        cube.Update()

        util.writeStl(cube.GetOutput(),"testCube")
        inPath = os.path.join(self.dirname, "../out/3D/testCube.stl")
        outPath = os.path.join(self.dirname, "../out/3D/testCube.off")
        util.meshioIO(inPath, outPath)

        meshPath = os.path.join(self.dirname, "../out/3D/papermeshLevel1.off")

        bool = boolean_interface.Boolean_Interface()
        bool.union(meshPath, outPath)

        self.intersection()

    def intersection(self):

        meshPath = os.path.join(self.dirname, "../out/3D/papermeshLevel0.off")
        cutoutPath = os.path.join(self.dirname, "../out/3D/union.off")

        bool = boolean_interface.Boolean_Interface()
        bool.boolean(meshPath,cutoutPath)


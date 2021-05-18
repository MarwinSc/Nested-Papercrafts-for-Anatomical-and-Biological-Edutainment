import os
import trimesh
from mim.mim import meshB_inside_meshA
import vtkmodules.all as vtk
import util
from mu3d.mu3dpy.mu3d import Graph


class HierarchicalMesh(object):
    """
    Represents a hierarchy of meshes. Tree structure.
    Children linked to parents and parents to children.
    A hierarchical mesh is associated with a paper mesh.
    Additionally it can have n-meshes that are projected on a paper mesh.
    """
    dirname = os.path.dirname(__file__)

    def __init__(self, parent, mesh, filename):
        """
        Initialises a hierarchical mesh.
        :param self: this
        :param mesh: mesh associated with this
        :param name: name of the mesh
        :return: None
        """
        self.mesh = mesh  # paper mesh
        self.children = []  # chain of paper meshes (inside of this paper mesh)
        self.projections = []  # add structures that are projected on this paper mesh here
        self.generated_by_plane_cut = ["../viewpointselection/data/urdf/bunny-front.urdf",
                                       "../viewpointselection/data/urdf/bunny-back.urdf",
                                       "../viewpointselection/data/urdf/bunny-sphere.urdf"]
        self.parent = parent  # parent paper mesh
        if parent is None:
            self.name = filename
        else:
            self.name = filename[filename.rfind('/') + 1:]
            self.offname = filename[:filename.rfind('.')] + ".off"
            print(self.name, "added to ", parent.name)

        self.file = filename
        self.dirname = os.path.dirname(__file__)

    '''
    def __init__(self, parent, meshes):
        """
        Initialises a hierarchical mesh.
        :param self: this
        :param parent: parent hierarchical mesh.
        :param meshes: mesh or list of meshes associated with this level
        :return: None
        """

        if isinstance(meshes,list):
            self.meshes = meshes
        else:
            self.meshes = [meshes]

        self.mesh = self.generatePaperMesh()
        self.children = []
        self.parent = parent

        filename = self.writePapermeshStlAndOff()

        if parent is None:
            self.name = filename[filename.rfind('/')+1:]
        else:
            self.name = filename[filename.rfind('/')+1:]
            self.offname = filename[:filename.rfind('.')] + ".off"
            print(self.name, "added to ", parent.name)

        self.file = filename
    '''

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

    def inside(self, file):
        """
        Checks if the given mesh is inside this mesh.
        :param file: file of mesh that is checked
        :return: True if the given mesh is inside this mesh
        """
        meshB = file[:file.rfind('.')] + ".off"
        return meshB_inside_meshA(bytes(self.offname, 'utf-8'), bytes(meshB, 'utf-8'))

    # TODO: at best here off file needs to exist!!
    def add(self, mesh, filename):
        """
        Adds the given mesh to the hierarchy
        :param mesh: Mesh that should be added to the hierarchy
        :param filename: Filename of the mesh
        :return: None
        """
        # is the mesh inside any of the children?
        # if so add it to the first child we encounter
        for child in self.children:
            if child.add(mesh, filename):
                return True

        # if this object does not have a mesh
        # it's top level therefore we can add any mesh here
        # if it was not added to any children
        # or we are not top level and it's inside the current mesh
        if self.mesh is None or self.inside(filename):
            # this means any of this children is potentially a child of the new mesh
            temp = self.children.copy()
            self.children.clear()
            self.children.append(HierarchicalMesh(self, mesh, filename))
            for tmp_child in temp:
                filename_off = filename[:filename.rfind('.')] + ".off" # here we need the off file to check the hierarchy
                if meshB_inside_meshA(bytes(filename_off, 'utf-8'), bytes(tmp_child.offname, 'utf-8')):
                    self.children[0].children.append(tmp_child)
                    tmp_child.parent = self.children[0]
                else:
                    self.children.append(tmp_child)

            return True

        return False

    # TODO change this to work on the simplified convex shapes
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
                # final_mesh = final_mesh.difference(tri_child, engine='blender')
                final_mesh = trimesh.boolean.difference([final_mesh, tri_child], engine="blender")

            filename = os.path.join(self.dirname, "../out/3D/differenced_" + self.name)
            final_mesh.export(filename)

        # if all children are cut out save it and call boolean for children
        for child in self.children:
            child.recursive_difference()

    def appendMesh(self, mesh):
        '''
        Append a mesh-object to the self.meshes list and update the papermesh.
        :param mesh:
        :return:
        '''
        self.meshes.append(mesh)
        self.mesh = self.generatePaperMesh()

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
        mesh = util.shrinkWrap(mesh, polysAppended)
        self.papermesh = util.offsetMesh(mesh)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputData(self.papermesh)
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetOpacity(0.15)
        return actor

    def unfoldWholeHierarchy(self, meshProcessor, iterations):
        self.unfoldPaperMesh(meshProcessor, iterations)
        for cild in self.children:
            cild.unfoldWholeHierarchy(meshProcessor, iterations)

    def unfoldPaperMesh(self, meshProcessor, iterations):
        '''
        Calls the mu3dUnfoldPaperMesh() method of the given meshProcessor and afterward the createDedicatedPaperMesh()
        to create a projectionMesh for each mesh in self.meshes.
        :param meshProcessor:
        :param iterations:
        :return:
        '''
        self.writePapermeshStlAndOff()
        self.graph = Graph()
        unfoldedActor = meshProcessor.mu3dUnfoldPaperMesh(self.mesh, self.graph, iterations)
        if unfoldedActor:
            self.unfoldedActor = unfoldedActor
            meshProcessor.createDedicatedMeshes(self)

    def getAllMeshes(self, asActor=True):
        '''
        :param asActor: if true the meshes are returned as vtk actors, if false as Mesh-objects.
        :return: All meshes of this level and below as vtk actors or Mesh-objects.
        '''
        if asActor:
            actors = [m.getActor() for m in self.meshes]
        else:
            actors = [m for m in self.meshes]

        for child in self.children:
            actors.extend(child.getAllMeshes())
        return actors

    def writePapermeshStlAndOff(self):
        '''
        Writes the papermesh in .stl format and .off format to the disk.
        named papermesh{level}.
        :return: The full path to the stl file.
        '''
        name = "papermeshLevel{}".format(self.getLevel())
        util.writeStl(self.papermesh, name)
        inpath = os.path.join(self.dirname, "../out/3D/papermeshLevel{}.stl".format(self.getLevel()))
        outpath = os.path.join(self.dirname, "../out/3D/papermeshLevel{}.off".format(self.getLevel()))
        util.meshioIO(inpath, outpath)
        return inpath

    def toString(self):
        '''
        Prints the Hierarchy Tree in the console.
        :return:
        '''
        print("level: ", self.getLevel())
        names = [m.filename.split("/")[-1] for m in self.meshes]
        for name in names:
            print(name, end=", ")
        print("")
        for child in self.children:
            child.toString()

    def getLevel(self):
        """
        :return: the level of this hierarchical mesh.
        """
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

        hierarchical_dict[current_level].append(self.mesh)

        for child in self.children:
            child.toLevelDict(hierarchical_dict)

        return hierarchical_dict

    def toList(self):
        level_dict = self.toLevelDict({})
        del level_dict[0]
        return [v for k, v in level_dict.items()]

    def convexify(self):
        """
        Saves convex versions of the papermesh to self.convex_paper_mesh
        :return: None
        """
        pass

    def simplify(self):
        """
        Saves simplified versions of self.convex_paper_mesh to self.simplified_paper_mesh
        :return: None
        """
        pass

    def cut_with_plane(self, plane):
        """
        Saves all objects generated by the cut to a list self.generated_by_plane_cut
        :param plane: plane defined as Ax + By + Cz + D
        :return:
        """
        pass

    def unfold_hierarchical_model(self, path):
        """
        Unfolds all objects from self.generated_by_plane_cut.
        :param path: path where unfolded models are saved
        :return:
        """
        pass

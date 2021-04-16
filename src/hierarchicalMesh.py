import os
import trimesh
from mim.mim import meshB_inside_meshA


class HierarchicalMesh(object):
    """
    Represents a hierarchy of meshes. Tree structure.
    Children linked to parents and parents to children.
    Each HierarchicalMesh is associated with a vtkActor.
    """

    def __init__(self, parent, mesh, filename):
        """
        Initialises a hierarchical mesh.
        :param self: this
        :param mesh: mesh associated with this
        :return: None
        """
        self.mesh = mesh
        self.children = []
        self.parent = parent

        if parent is None:
            self.name = filename
        else:
            self.name = filename[filename.rfind('/')+1:]
            self.offname = filename[:filename.rfind('.')] + ".off"
            print(self.name, "added to ", parent.name)

        self.file = filename
        self.dirname = os.path.dirname(__file__)

    def render(self, level, renderer):
        """
        Adds actors to the given renderer based on the level.
        :param level:
        :param renderer:
        :return:
        """
        if level == 0:
            if self.mesh is not None:
                # renderer.AddActor(self.mesh)
                self.mesh.SetVisibility(True)

            for child in self.children:
                child.render(level, renderer)
        else:
            if self.mesh is not None:
                self.mesh.SetVisibility(False)

            for child in self.children:
                child.render(level-1, renderer)

    def inside(self, file):
        """
        Checks if the given mesh is inside this mesh.
        :param file: file of mesh that is checked
        :return: True if the given mesh is inside this mesh
        """
        meshB = file[:file.rfind('.')] + ".off"
        print(self.offname, meshB)
        return meshB_inside_meshA(bytes(self.offname, 'utf-8'), bytes(meshB, 'utf-8'))

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
        # or we are nto top level and it's inside the current mesh
        if self.mesh is None or self.inside(filename):
            # this means any of this children is potentially a child of the new mesh
            temp = self.children.copy()
            self.children.clear()
            self.children.append(HierarchicalMesh(self, mesh, filename))
            for tmp_child in temp:
                if not self.children[0].add(tmp_child.mesh, tmp_child.file):
                    self.children.append(tmp_child)

            return True

        return False


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

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
        :param name: name of the mesh
        :return: None
        """
        self.mesh = mesh
        self.children = []
        self.parent = parent
        if parent is None:
            self.name = filename
        else:
            self.name = filename[filename.rfind('/')+1:]
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

    def inside(self, mesh):
        """
        Checks if the given mesh is inside this mesh.
        :param mesh: mesh that is checked
        :return: True if the given mesh is inside this mesh
        """

        result = meshB_inside_meshA("", "")

        return True

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
            if child.inside(mesh):
                child.add(mesh, filename)
                return

        # we only get here if it is not inside any of the children
        # check if the mesh is inside this mesh?
        # if not raise an error for now, otherwise the hierarchy needs
        # to be reorganized, ignored for now
        if not self.inside(mesh):
            raise RuntimeError('Mesh is outside of this mesh..')

        self.children.append(HierarchicalMesh(self, mesh, filename))

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

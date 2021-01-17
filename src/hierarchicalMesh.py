import vtkmodules.all as vtk


class HierarchicalMesh(object):
    """
    """

    def __init__(self, parent, mesh):
        """
        Initialises a hierarchical mesh.
        :param self: this
        :param mesh: mesh associated with this
        :return: None
        """
        self.mesh = mesh
        self.children = []
        self.parent = parent

    def inside(self, mesh):
        """
        Checks if the given mesh is inside this mesh.
        :param mesh: mesh that is checked
        :return: True if the given mesh is inside this mesh
        """

        return True

    def add(self, mesh):
        """
        Adds the given mesh to the hierarchy
        :param mesh: Mesh that should be added to the hierarchy
        :return: None
        """
        # is the mesh inside any of the children?
        # if so add it to the first child we encounter
        for child in self.children:
            #if child.self.inside(mesh):
            child.add(mesh)
            return

        # we only get here if it is not inside any of the children
        # check if the mesh is inside this mesh?
        # if not raise an error for now, otherwise the hierarchy needs
        # to be reorganized, ignored for now
        #if not self.inside(mesh):
        #    raise RuntimeError('Mesh is outside of this mesh..')

        self.children.append(HierarchicalMesh(self, mesh))

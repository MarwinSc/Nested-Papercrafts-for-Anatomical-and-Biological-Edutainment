class HierarchicalMesh(object):
    """
    Represents a hierarchy of meshes. Tree structure.
    Children linked to parents and parents to children.
    Each HierarchicalMesh is associated with a vtkActor.
    """

    def __init__(self, parent, mesh, name):
        """
        Initialises a hierarchical mesh.
        :param self: this
        :param mesh: mesh associated with this
        :return: None
        """
        self.mesh = mesh
        self.children = []
        self.parent = parent
        self.name = name

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

        return True

    def add(self, mesh, name):
        """
        Adds the given mesh to the hierarchy
        :param mesh: Mesh that should be added to the hierarchy
        :return: None
        """
        # is the mesh inside any of the children?
        # if so add it to the first child we encounter
        for child in self.children:
            if child.inside(mesh):
                child.add(mesh, name)
                return

        # we only get here if it is not inside any of the children
        # check if the mesh is inside this mesh?
        # if not raise an error for now, otherwise the hierarchy needs
        # to be reorganized, ignored for now
        if not self.inside(mesh):
            raise RuntimeError('Mesh is outside of this mesh..')

        print(name, "added to ", self.name)
        self.children.append(HierarchicalMesh(self, mesh, name))

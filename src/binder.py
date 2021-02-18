import ctypes
import meshio
import vtkmodules.all as vtk
from vtkmodules.numpy_interface.dataset_adapter import numpy_support
import os
from sys import platform

class Binder():

    dirname = os.path.dirname(os.path.abspath(__file__))

    def __init__(self):
        if platform == "darwin":
            filepath = os.path.join(self.dirname, "libMU3D.so")
            self.mu3d = ctypes.CDLL(filepath)
        elif platform == "win32":
            filepath = os.path.join(self.dirname, "mu3d.dll")
            self.mu3d = ctypes.CDLL(filepath)
            #self.mu3d.python_interface.argtypes = [ctypes.c_int,ctypes.c_int]
        else:
            print("Sorry, No Linux Support")


    def onUnfold(self):

        filename = os.path.join(self.dirname, "../out/3D/papermesh.stl")

        mesh = meshio.read(
            filename,  # string, os.PathLike, or a buffer/open file
            file_format="stl",  # optional if filename is a path; inferred from extension
        )

        filename = os.path.join(self.dirname, "../out/3D/papermesh.off")

        meshio.write(
            filename,  # str, os.PathLike, or buffer/ open file
            mesh,
            # file_format="vtk",  # optional if first argument is a path; inferred from extension
        )

        print(filename)
        print(self.mu3d.python_interface(30000,1000))


    def unfoldTest(self):

        filename = os.path.join(self.dirname, "../out/3D/upper.stl")
#        filename = os.path.join(self.dirname, "../out/3D/papermesh_cleaned.off")


        mesh = meshio.read(
            filename,  # string, os.PathLike, or a buffer/open file
            file_format="stl",  # optional if filename is a path; inferred from extension
        )

        filename = os.path.join(self.dirname, "../out/3D/papermesh.off")

        meshio.write(
            filename,  # str, os.PathLike, or buffer/ open file
            mesh,
            # file_format="vtk",  # optional if first argument is a path; inferred from extension
        )

#       print(filename)
        print(self.mu3d.python_interface(200000,1000))

#       importer = vtk.vtkOBJReader()
#       filename = os.path.join(self.dirname, "../out/3D/unfolded/model.obj")
#       importer.SetFileName(filename)
#       importer.Update()
#       mesh = importer.GetOutput()

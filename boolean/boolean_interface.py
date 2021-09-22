import ctypes
import os

class Boolean_Interface(object):
    """
    Wrapper for CGAL to use boolean and remesh functionality
    """

    def __init__(self):
        """

        """
        dir_path = os.path.dirname(os.path.realpath(__file__))
        self.boolean_interface = ctypes.CDLL(dir_path + "/boolean_interface.dll")
        self.boolean_interface._boolean_interface.restype = ctypes.c_void_p
        self.boolean_interface._boolean.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_char_p]
        self.boolean_interface._boolUnion.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_char_p]
        self.boolean_interface._triangulateCut.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_float,ctypes.c_float,ctypes.c_float, ctypes.c_float,ctypes.c_float,ctypes.c_float]
        self.boolean_interface._merge.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_float, ctypes.c_float,ctypes.c_float,ctypes.c_float, ctypes.c_float,ctypes.c_float,ctypes.c_float]

        self.obj = self.boolean_interface._boolean_interface()

    def merge_adjacent_vertices_by_distance(self,path, threshold, normal_x,normal_y,normal_z, origin_x, origin_y, origin_z):

        self.boolean_interface._merge(self.obj, path.encode(), threshold, normal_x, normal_y, normal_z, origin_x, origin_y, origin_z)

    def triangulateCut(self,path,normal_x,normal_y,normal_z, origin_x, origin_y, origin_z):
        """

        :param path: Path to mesh piece that should be triangulated
        :param normal_x: x-coordinate of the cutting planes normal
        :param normal_y: y-coordinate of the cutting planes normal
        :param normal_z: z-coordinate of the cutting planes normal
        :return:
        """
        self.boolean_interface._triangulateCut(self.obj, path.encode(), normal_x, normal_y, normal_z, origin_x, origin_y, origin_z)

    def boolean(self,firstpath,secondpath):
        """

        """
        self.boolean_interface._boolean(self.obj, firstpath.encode(), secondpath.encode())

    def union(self,firstpath,secondpath):
        self.boolean_interface._boolUnion(self.obj, firstpath.encode(), secondpath.encode())
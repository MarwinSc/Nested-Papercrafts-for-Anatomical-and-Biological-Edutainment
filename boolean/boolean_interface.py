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


        self.obj = self.boolean_interface._boolean_interface()

    def boolean(self,firstpath,secondpath):
        """

        """
        self.boolean_interface._boolean(self.obj, firstpath.encode(), secondpath.encode())

    def union(self,firstpath,secondpath):
        self.boolean_interface._boolUnion(self.obj, firstpath.encode(), secondpath.encode())
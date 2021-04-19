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
        #todo rename cgal to boolean
        self.boolean_interface = ctypes.CDLL(dir_path + "/cgal_interface.dll")
        self.boolean_interface._cgal_interface.restype = ctypes.c_void_p
        self.boolean_interface._boolean.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_char_p]

        self.obj = self.boolean_interface._cgal_interface()

    def boolean(self,firstpath,secondpath):
        """

        """
        self.boolean_interface._boolean(self.obj, firstpath.encode(), secondpath.encode())
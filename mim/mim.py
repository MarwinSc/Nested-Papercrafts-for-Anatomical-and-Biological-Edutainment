import ctypes
import os
from sys import platform

def loadLib():
    dir_path = os.path.dirname(os.path.realpath(__file__))
    if platform == "linux" or platform == "linux2" or platform == "darwin":
        lib = ctypes.CDLL(dir_path + "/libmim.so")
    elif platform == "win32":
        lib = ctypes.WinDLL(dir_path + "/mim.dll")
    else:
        raise Exception("unsupported platform")
    return lib


def meshB_inside_meshA(fileA, fileB):
    tmp = loadLib()
    tmp.meshB_inside_of_meshA.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
    tmp.meshB_inside_of_meshA.restype = ctypes.c_bool
    return tmp.meshB_inside_of_meshA(fileA, fileB)


def hull_of_mesh(file_in, file_out):
    tmp = loadLib()
    tmp.convex_hull_of_mesh.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
    tmp.convex_hull_of_mesh.restype = ctypes.c_bool
    return tmp.convex_hull_of_mesh(file_in, file_out)


def simplify_mesh(file_in, file_out, simplification_rate):
    tmp = loadLib()
    tmp.simplify_mesh.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_float]
    tmp.simplify_mesh.restype = ctypes.c_bool
    return tmp.simplify_mesh(file_in, file_out, simplification_rate)


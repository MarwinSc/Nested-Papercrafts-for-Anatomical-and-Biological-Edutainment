import ctypes
import os


def meshB_inside_meshA(fileA, fileB):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    tmp = ctypes.WinDLL(dir_path + "/mim.dll")
    tmp.meshB_inside_of_meshA.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
    tmp.meshB_inside_of_meshA.restype = ctypes.c_bool
    return tmp.meshB_inside_of_meshA(fileA, fileB)


def hull_of_mesh(file_in, file_out):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    tmp = ctypes.WinDLL(dir_path + "/mim.dll")
    tmp.convex_hull_of_mesh.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
    tmp.convex_hull_of_mesh.restype = ctypes.c_bool
    return tmp.convex_hull_of_mesh(file_in, file_out)

import ctypes
import os


def meshB_inside_meshA(fileA, fileB):
    dir_path = os.path.dirname(os.path.realpath(__file__))
    tmp = ctypes.WinDLL(dir_path + "/mim.dll")
    tmp.meshB_inside_of_meshA.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
    tmp.meshB_inside_of_meshA.restype = ctypes.c_bool
    tmp.meshB_inside_of_meshA(fileA, fileB)

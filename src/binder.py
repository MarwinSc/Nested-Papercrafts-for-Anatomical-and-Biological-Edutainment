import ctypes

class Binder():

    def init(self):
        lib = ctypes.CDLL('lib.so')
        counter = lib.hello(4)
        print(counter)
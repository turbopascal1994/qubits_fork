cimport cython
from libcpp.vector cimport vector

cdef extern from "./hello.h":
    vector[int] hello()

def call_hello():
    return hello()

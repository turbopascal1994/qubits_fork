from libcpp.vector cimport vector
from libcpp cimport complex

cdef extern from "kernel.h":
    cdef cppclass Kernel:
        Kernel(double, double, double, double, double, double, int) except +
        vector[int] CreateStartSCALLOP(int, double)

cdef extern from "test.h":
    cdef cppclass TestCPP:
        TestCPP(double, double, double) except +
        vector[double complex] get()
        void eig(int)

cdef class PyTestCPP:
    cdef TestCPP * c_ptr

    def __cinit__(self, double h, double w01, double w12):
        self.c_ptr = new TestCPP(h, w01, w12)

    def __dealloc__(self):
        del self.c_ptr

    def get(self):
        return self.c_ptr.get()
    
    def eig(self, int n):
        return self.c_ptr.eig(n)


cdef class PyKernel:

    cdef Kernel * c_ptr

    def __cinit__(
        self,
        double tstep,
		double w01,
		double w12,
		double wt,
		double w,
		double T,
		int Type
    ):
        self.c_ptr = new Kernel(tstep, w01, w12, wt, w, T, Type)

    def __dealloc__(self):
        del self.c_ptr

    def CreateStartSCALLOP(self, int M, double AmpThreshold):
        return self.c_ptr.CreateStartSCALLOP(M, AmpThreshold)

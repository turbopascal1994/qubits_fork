import os

path_to_libs = "D:\\c++\\NNSU\\Linev\\matlab\\code\\python_experiments\\OHFSS"

os.add_dll_directory(path_to_libs)
os.add_dll_directory("C:/Program Files (x86)/Intel/oneAPI/mkl/2022.0.2/lib/intel64")
os.add_dll_directory("C:/Program Files (x86)/Intel/oneAPI/mkl/latest/redist/intel64")

os.environ['PATH'] = path_to_libs + os.pathsep + os.environ['PATH']

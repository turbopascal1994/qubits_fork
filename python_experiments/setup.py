from xml.etree.ElementInclude import include
from setuptools import setup, Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import glob
import os
from os.path import join as jp

ext_1 = Extension(
    "OHFSS.Hello",
    sources=["src/hello/hello.cpp", "src/hello/test.pyx"],
    include_dirs=["src/hello"],
    language="c++"
)

ext_2 = Extension(
    "OHFSS.Kernel",
    sources=[
        "src/kernel/linalg.cpp",
        "src/kernel/kernel.cpp",
        "src/kernel/PyKernel.pyx",
    ],
    depends=glob.glob(jp(os.path.abspath('src'), '*.h')),
    include_dirs=[
        "D:\\c++\\NNSU\\Linev\\matlab\\code\\python_experiments\\src\\hello",
        "D:\\c++\\NNSU\\Linev\\matlab\\code\\python_experiments\\src\\kernel",
        "C:/Program Files (x86)/Intel/oneAPI/mkl/2022.0.2/include",
        "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/include/oneapi",
        "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/include/oneapi/mkl"
    ],
    extra_compile_args=[
        '-wd4267',
        '-wd4244',
        '-wd4101',
        '-wd4996',
        '/std:c++17',
        '-GS'
    ],
    extra_link_args=[
        '-NXCompat',
        '-DynamicBase',
        '-IGNORE:4197'
    ],
    library_dirs=[
        "C:/Program Files (x86)/Intel/oneAPI/mkl/2022.0.2/lib/intel64",
        "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/redist/intel64",
    ],
    libraries=[
        "mkl_blacs_ilp64_dll",
        "mkl_blacs_lp64_dll",
        "mkl_blas95_ilp64",
        "mkl_blas95_lp64",
        "mkl_cdft_core",
        "mkl_cdft_core_dll",
        "mkl_core",
        "mkl_core_dll",
        "mkl_intel_ilp64",
        "mkl_intel_ilp64_dll",
        "mkl_intel_lp64",
        "mkl_intel_lp64_dll",
        "mkl_intel_thread_dll",
        "mkl_lapack95_ilp64",
        "mkl_lapack95_lp64",
        "mkl_pgi_thread",
        "mkl_pgi_thread_dll",
        "mkl_rt",
        "mkl_scalapack_ilp64",
        "mkl_scalapack_ilp64_dll",
        "mkl_scalapack_lp64",
        "mkl_scalapack_lp64_dll",
        "mkl_sequential",
        "mkl_sequential_dll",
        "mkl_sycl",
        "mkl_sycl_dll",
        "mkl_sycld",
        "mkl_sycld_dll",
        "mkl_tbb_thread",
        "mkl_tbb_thread_dll",
        "mkl_tbb_threadd",
        "mkl_tbb_threadd_dll",
    ],
    language="c++"
)

EXTENSIONS = [ext_1, ext_2]

setup(
    name="OHFSS",
    cmdclass={"build_ext": build_ext},
    ext_modules=EXTENSIONS,
)

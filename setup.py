from distutils.core import setup, Extension

from Cython.Distutils import build_ext
import pkg_resources

data_dir = pkg_resources.resource_filename("autowrap", "data_files")

ext = Extension("pllpy",
                sources = ['src/pllpy.pyx', 'src/pllml.cpp'],
                language="c++",
                include_dirs = [data_dir],
                libraries=['pll-sse3-pthreads'],
                extra_compile_args=['-std=c++11'],
               )

setup(cmdclass={'build_ext':build_ext},
      name="pllpy",
      version="0.0.1",
      ext_modules = [ext],
      install_requires=[
        'autowrap',
        'cython',
    ],
     )
try:
    from setuptools import setup, Extension
except ImportError:
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
      version="0.1.13",
      author='Kevin Gori',
      author_email='kgori@ebi.ac.uk',
      description='Wrapper for Phylogenetic Likelihood Library',
      url='https://github.com/kgori/pllpy.git',
      ext_modules = [ext],
      install_requires=[
        'autowrap',
        'cython',
      ],
      scripts=['bin/pll']
     )

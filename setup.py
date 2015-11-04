try:
    from setuptools import setup, Extension
except ImportError:
    from distutils.core import setup, Extension

from Cython.Distutils import build_ext
import pkg_resources
import platform, re, subprocess

data_dir = pkg_resources.resource_filename("autowrap", "data_files")

def is_clang(bin):
    proc = subprocess.Popen([bin, '-v'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    output = '\n'.join([stdout, stderr])
    return not re.search(r'clang', output) is None

class my_build_ext(build_ext):
    def build_extensions(self):
        binary = self.compiler.compiler[0]
        if is_clang(binary):
            for e in self.extensions:
                e.extra_compile_args.append('-stdlib=libc++')
                if platform.system() == 'Darwin':
                    e.extra_compile_args.append('-mmacosx-version-min=10.7')
                    e.extra_link_args.append('-mmacosx-version-min=10.7')

        build_ext.build_extensions(self)

ext = Extension("pllpy",
                sources = ['src/pllpy.pyx', 'src/pllml.cpp'],
                language="c++",
                include_dirs = [data_dir],
                libraries=['pll-sse3'],
                extra_compile_args=['-std=c++11'],
               )

setup(cmdclass={'build_ext':my_build_ext},
      name="pllpy",
      version="0.1.16",
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

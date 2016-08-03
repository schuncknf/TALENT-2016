#!/usr/bin/env python

from distutils.core import setup, Extension
from distutils.cmd import Command
shf_dir = "../.."
shf_build_dir = shf_dir + ""
shf_include_dir = shf_dir + "/cpp/include"
shf = Extension('_shf', ['shf.i'],
                include_dirs = [shf_include_dir ,'/usr/include/python2','/usr/include','./armanpy-0.1.4/include'],
                libraries = ['m', 'z', 'armadillo', 'gsl', 'gslcblas'],
                extra_compile_args = ['-O3', '-march=native', '-DARMA_NO_DEBUG', '-std=gnu++11'],
                #extra_compile_args = ['-g', '-fno-inline-functions', '-O0', '-D_GLIBCXX_DEBUG', '-std=gnu++11'],
                extra_objects = [shf_build_dir + '/libshf.a'],
                depends = [shf_build_dir + '/libshf.a']
                )

setup(
    name = 'shf',
    version = '0.1a',
    author = 'DamnDucksTowels',
    description = """Hartree-Fock lib for python""",
    ext_modules = [ shf ],
    py_modules = [ 'shf' ],
#    cmdclass={
#        'test': TestClass,
#    },
)

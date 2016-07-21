#!/usr/bin/env python

from distutils.core import setup, Extension
from distutils.cmd import Command
hffs_dir = "../../cpp"
hffs_build_dir = hffs_dir + ""
hffs_include_dir = hffs_dir + "/include"
hffs = Extension('_hffs', ['hffs.i'],
                include_dirs = [hffs_include_dir ,'/usr/include/python2','/usr/include','./armanpy-0.1.4/include'],
                libraries = ['m', 'z', 'armadillo'],
                #extra_compile_args = ['-O3', '-march=native', '-DARMA_NO_DEBUG'],
                extra_compile_args = ['-g', '-fno-inline-functions', '-O0', '-D_GLIBCXX_DEBUG'],
                extra_objects = [hffs_build_dir + '/libhffs.a'],
                depends = [hffs_build_dir + '/libhffs.a']
                )

setup(
    name = 'hffs',
    version = '0.1a',
    author = 'DamnDucksTowels',
    description = """Hartree-Fock lib for python""",
    ext_modules = [ hffs ],
    py_modules = [ 'hffs' ],
#    cmdclass={
#        'test': TestClass,
#    },
)

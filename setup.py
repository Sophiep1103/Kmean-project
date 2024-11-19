from setuptools import setup, Extension

module = Extension("mykmeanssp", sources=["spkmeansmodule.c", "spkmeans.c"])

setup(
    name = 'mykmeanssp',
    version = '1.0',
    description = 'spkmeans',
    ext_modules = [module]
)

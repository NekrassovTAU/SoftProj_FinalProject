from setuptools import setup, find_packages, Extension

setup(
    name='spkmeansmodule',
    version='1.0.0',
    author='Daniel Nekrassov and Lior Grinberg',
    description='spkmeans module created in C',
    install_requires=['invoke'],
    packages=find_packages(),
    ext_modules=[Extension('spkmeansmodule', ['spkmeansmodule.c', 'spkmeans.c'],
                           include_dirs=['include'])]
)
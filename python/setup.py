from distutils.core import setup
from distutils.extension import Extension
import os
import sys
import platform
import numpy

openmm_dir = '@OPENMM_DIR@'
imageplugin_header_dir = '@IMAGEPLUGIN_HEADER_DIR@'
imageplugin_library_dir = '@IMAGEPLUGIN_LIBRARY_DIR@'

# setup extra compile and link arguments on Mac
extra_compile_args = ['-std=c++11']
extra_link_args = []

if platform.system() == 'Darwin':
    extra_compile_args += ['-stdlib=libc++']
    extra_link_args += ['-stdlib=libc++', '-Wl', '-rpath', openmm_dir+'/lib']
    if 'MACOSX_DEPLOYMENT_TARGET' not in os.environ and platform.processor() != 'arm':
                extra_compile_args += ['-mmacosx-version-min=10.7']
                extra_link_args += ['-mmacosx-version-min=10.7']

extension = Extension(name='_imageplugin',
                      sources=['ImagePluginWrapper.cpp'],
                      libraries=['OpenMM', 'OpenMMImage'],
                      include_dirs=[os.path.join(openmm_dir, 'include'), imageplugin_header_dir, numpy.get_include()],
                      library_dirs=[os.path.join(openmm_dir, 'lib'), imageplugin_library_dir],
                      extra_compile_args=extra_compile_args,
                      extra_link_args=extra_link_args
                     )

setup(name='imageplugin',
      version='1.1',
      py_modules=['imageplugin'],
      ext_modules=[extension],
     )

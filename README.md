OpenMM Image Plugin
=====================

This plugin is based on the OpenMM 7.7.0, the minimum version of CMake is 3.17. 

Building The Plugin
===================

This project uses [CMake](http://www.cmake.org) for its build system.  To build it, follow these
steps:

1. Create a directory in which to build the plugin.

2. Set environmental variables such as CXXFLAGS='-std=c++11', CC=gcc, CXX=g++, 
OPENMM_CUDA_COMPILER=$(which nvcc).

3. Run the CMake GUI or ccmake, specifying your new directory as the build directory and the top
level directory of this project as the source directory. (This CMakeLists.txt only supports building 
dynamic libraries.)

4. Press "Configure".

5. Set OPENMM_DIR to point to the directory where OpenMM is installed.  This is needed to locate
the OpenMM header files and libraries.

6. Set CMAKE_INSTALL_PREFIX to the directory where the plugin should be installed.  Usually,
this will be the same as OPENMM_DIR, so the plugin will be added to your OpenMM installation.

7. OpenCL Platform is not supported.

8. If you plan to build the CUDA platform, make sure that CUDA_TOOLKIT_ROOT_DIR is set correctly
and that IMAGE_BUILD_CUDA_LIB is selected.

9. Press "Configure" again if necessary, then press "Generate".

10. Use the build system you selected to build and install the plugin. 
Performing three sets of commands 'make / make install / make PythonInstall' will install the plugin
as 'imageplugin' package in your python. The swig version should be higher than 3.0.

Note: If you want to use the ImageCustomIntegrator in OpenMM 8.0, you need to 
copy the folder asmjit of OpenMM to the library folder in the plugin for substitution. 
The CMakeList.txt also needs to be modified accordingly. 
If you update the imageplugin version in setup.py, you should use pip uninstall the old version, 
then install the new version.

Usage
==========

The image charge method is an effienct method to deal with the surface polarization. Here we reconstructed 
the image charge method in OpenMM 7.7.0 or higher versions with Langevin and Custom Integrators. Note that 
the plugin can only used in the system with two parallel conductor planes. Refer to 
[this articles](https://doi.org/10.1073/pnas.2020615118), 
[this articles](https://doi.org/10.1063/5.0040172), 
[this articles](https://pubs.acs.org/doi/10.1021/acs.jpcc.9b06635) for more details of this method.

### ImageLangevinIntegrator
```python
    from imageplugin import ImageLangevinIntegrator

    integrator = ImageLangevinIntegrator(temperature, freq, timestep)
    integrator.setCellSize(zbox)  # zbox is the distance between two surfaces
    nRealAtoms = system.getNumParticles()
    nbforce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == NonbondedForce][0]

    for i in range(nRealAtoms):
        (q, sig, eps) = nbforce.getParticleParameters(i)
        newAtom = topology.addAtom('IM', None, newResidue)
        idxat = system.addParticle(0*dalton)
        idxat2 = nbforce.addParticle(-q, imsig, imeps)
        # add image pairs
        integrator.setImagePair(idxat, i)
```

### ImageCustomIntegrator
This integrator with same function of CustomIntegrator of OpenMM offers additional functionality for moving 
image particles. It can be used to construct any type of integrator, such as Velocity-Verlet integrator or 
Langevin integrator, with the added capability to move image particles.
```
from imageplugin import CustomIntegrator
```

### ImageMTSLangevinIntegrator

This integraotr which is provied in the script imagemtsintegrator.py is inherited from 
the ImageCustomIntegrator. The MTSLangevinIntegrator implements the RESPA multiple time 
step algorithm, which is usually used in conjunction with the AMOEBA force field.

```python
from pytools/imagemtsintegrator import *

for f in system.getForces():
    if isinstance(f, AmoebaVdwForce):
        f.setForceGroup(1)
    elif isinstance(f, AmoebaMultipoleForce):
        f.setForceGroup(2)
imageInteg = ImageMTSLangevinIntegrator(temperature, freq, timestep, [(0,8), (1,1), (2,1)])
integrator.setCellSize(zbox)
mtpForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == AmoebaMultipoleForce][0]
for i in range(nRealAtoms):
    (q, dip, quad, axisType, atomZ, atomX, atomY, thole, dampFactor, polarity) = mtpForce.getMultipoleParameters(i)
    newAtom = topology.addAtom('IM', next(atoms).element, newResidue)
    idxat = system.addParticle(0*dalton)
    dip1 = Quantity((-dip[0], -dip[1], -dip[2])).in_units_of(nanometer*elementary_charge)
    quad1 = Quantity(tuple(-d for d in quad)).in_units_of(nanometer**2*elementary_charge)
    idxat2 = mtpForce.addMultipole(-q, dip1, quad1, axisType, atomZ+nRealAtoms, atomX+nRealAtoms, atomY+nRealAtoms,
                                   thole, dampFactor, polarity)
    # add image pairs
    imageInteg.setImagePair(idxat, i)
```

### Monte Carlo Barostat in Z-Direction

The MCZbarastae allows for controlled pressure only in the z-direction for a system with two 
parallel plates. Only one of plates is moved proportionally with the box.
Refer to 
[this article](https://pubs.acs.org/doi/10.1021/acs.jpcc.0c00299) for more details.
```python
from imageplugin import MCZBarostat

system.addForce(MCZBarostat(pressure, temperature, barostatInterval))
```

### SlabCorrection of the nonpolarization surface of Amoeba force field

The SlabCorrection fixes the system polariztion when using vaccum boundary condition and Amoeba force field. The 
SlabeCorretion takes into account not only the particle's charge dipole moments (q*z), but also the fixed and 
induced dipole moments of the Amoeba force field. There is an option (bool useAmoebaDip) of choosing whether to use 
dipole moments generated by Amoeba force field.
Refer to [this article](https://doi.org/10.1063/1.479595) for more details. 

```python
from imageplugin import SlabCorrection
# The first parameter is whether to use dipoles of all partcles to compute the slab correction. 
# The second parameter is whether to use dipole moments generated by Amoeba force field.
system.addForce(SlabCorrection(True, True))
```
Note: When the SlabCorrection is not in the same force group with AmoebaMultipoleFroce, the tuple of force group of 
SlabCorrection should be placed before the tuple of AmoebaMultipoleForce in integrators of RESPA algorithm class 
((Image)MSTIntegrator and (Image)MSTLangevinIntegrator).

### The scirpt for converting the Tinker prm file to OpenMM xml file.

The tinker2openmm.py is modified based on one scirpt file of [the repo](https://github.com/Inniag/openmm-scripts-amoeba). 
Here the modified script supports the Z-Bisector local axis type and some other
parameters, and adds some options which can be seen in detail with "--help".
Typically, use the script as follow.
```
python3 tinker2openmm.py -resname chcl3 -input_xyz=chloroform -input_prm=chcl3 -xml_bond=number
```

For more detailed usage, please refer to the scripts.


Examples and citation
=====================

Examples from the following work are provided in `examples` to demonstrate the usage of this plugin.
Please cite this article if you find this plugin useful.


License
=======

This is part of the OpenMM molecular simulation toolkit originating from
Simbios, the NIH National Center for Physics-Based Simulation of
Biological Structures at Stanford, funded under the NIH Roadmap for
Medical Research, grant U54 GM072970. See https://simtk.org.

Portions copyright (c) 2023 the Authors.

Authors: Ruochao Wang

Contributors: 
Part of the code comes from [scychon's openmm_constV](https://github.com/scychon/openmm_constV).

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the
Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
USE OR OTHER DEALINGS IN THE SOFTWARE.


import itertools
import numpy as np
from openmm.app import *
from openmm import *
from openmm.unit import *
import time
import sys
# from mdtools import vvintegrator
from imageplugin import ImageCustomIntegrator
# from constvplugin import *
from imagemtsintegrator import *

time_start = time.time()
start_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time_start))

# print(sys.getrecursionlimit())
# sys.setrecursionlimit(1500)

# integrator parameters
pressure = 1*atmospheres
meltT = 400*kelvin
annealT1 = 360*kelvin
annealT2 = 330*kelvin
temperature = 298 * kelvin
freq = 1.0 / picosecond
timestep = 1.0 * femtosecond
barostatInterval = 25
cutoff = 1.0 * nanometers
vdw_cutoff = 1.2 * nanometers
ewald_error_tolerance = 0.0005
polarization = 'mutual'
mutualInducedTargetEpsilon = 0.0005
maxIteration = 100
meltSteps = 2000000
annealSteps = 1000000
equilibrationSteps = 10000000
samplingSteps = 100000000
zmin = 0
zmax = 1.8623


# force field
ffField = 'pstfsi10_pAlt.xml'
forceField = ForceField(ffField)

# input files and output path
input_name = 'pstfsi10-h2o_cw51'
pdb = PDBFile('eq_pstfsi10Li-h2o_npt51.pdb')
output = 'output/'

#####################
# Specific setup for box dimension and potential difference of the electrodes
#####################
print('Building system...')
platform = Platform.getPlatformByName('CUDA')
platformProperties = {'Precision': 'mixed', 'DeviceIndex': '0'}


# set the electric filed
constV = 0*volt  # potential difference between two parallel electrodes
# electric field calculated for the gap between the two electrode
constEfield = (constV*elementary_charge*AVOGADRO_CONSTANT_NA/(zmax*nanometer)).in_units_of(kilojoule_per_mole/nanometer)

# ##########################
# Setup the simulation system topology and initial positions
##########################
topology = pdb.topology
position = pdb.positions

# Get the dimension values of the crystallographic unit cell
boxvec = pdb.topology.getUnitCellDimensions()._value
print(boxvec)
if (boxvec[2]-zmax*2) != 0:
    pdb.topology.setUnitCellDimensions(Vec3(boxvec[0], boxvec[1], zmax*2)*nanometer)

# create system with real particles only
system = forceField.createSystem(topology, nonbondedMethod=PME, nonbondedCutoff=cutoff, vdwCutoff=vdw_cutoff,
                                 ewaldErrorTolerance=ewald_error_tolerance, polarization=polarization,
                                 mutualInducedTargetEpsilon=mutualInducedTargetEpsilon,
                                 mutualInducedMaxIterations=maxIteration,
                                 constraints=HBonds, rigidWater=True)

nRealAtoms = system.getNumParticles()
vdwForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == AmoebaVdwForce][0]
mtpForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == AmoebaMultipoleForce][0]
hbForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == HarmonicBondForce][0]
abForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == CustomBondForce][0]
'''
add exclusion for electrode atoms
In this example only the wall atoms have zero mass.
If other real particles have zero mass, you need to specify other routine to set the exclusion between the wall atoms
This step is may be skipped since the interactions between the wall atoms would be a constant throughout
the simulation and does not affect the force on real atoms.
But excluding the wall interaction ensures legitimate value for the system potential energy
'''
# first find all the wall atoms
grpBottom = []
grpTop = []
for i in range(nRealAtoms):
    if system.getParticleMass(i) == 0*dalton:
        if position[i][2]==0*nanometer:
            grpBottom.append(i)
        else:
            grpTop.append(i)

# add exclusion between atoms in right electrode (at z=zmax)
for igrpTop in grpTop:
    exvdwParticleIndex = []
    for jgrpTop in grpTop:
        dpos = position[igrpTop] - position[jgrpTop]
        # periodic boundary condition
        dr = dpos - (np.round(np.asarray(dpos._value)/boxvec)*boxvec)*nanometer
        # ignore the vdw force between atoms with bonds
        if np.linalg.norm(dr._value)*nanometer < .5*nanometer and igrpTop != jgrpTop:
            exvdwParticleIndex.append(jgrpTop)
    vdwForce.setParticleExclusions(igrpTop, exvdwParticleIndex)

for igrpBot in grpBottom:
    exvdwParticleIndex = []
    for jgrpBot in grpBottom:
        dpos = position[igrpBot] - position[jgrpBot]
        # periodic boundary condition
        dr = dpos - (np.round(np.asarray(dpos._value)/boxvec)*boxvec)*nanometer
        # ignore the vdw force between atoms with bonds
        if np.linalg.norm(dr._value)*nanometer < .5*nanometer and igrpBot != jgrpBot:
            exvdwParticleIndex.append(jgrpBot)
    vdwForce.setParticleExclusions(igrpBot, exvdwParticleIndex)

#####################
#  set force group
#####################
for f in system.getForces():
    if (isinstance(f, AmoebaVdwForce)
            or isinstance(f, AmoebaGeneralizedKirkwoodForce)
            or isinstance(f, AmoebaWcaDispersionForce)):
        f.setForceGroup(1)
    elif isinstance(f, AmoebaMultipoleForce):
        f.setForceGroup(2)

# get force group
for i in range(system.getNumForces()):
    f = system.getForce(i)
    # f.setForceGroup(i)
    print(str(f.getForceGroup())+' '+str(type(f)))

#####################
#  Warming up
#####################
print("Melting...")
integrator = MTSLangevinIntegrator(meltT, freq, timestep, [(0,4), (1,1), (2,1)])
simulation = Simulation(topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(position)
state_init = simulation.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True)
# minimize and equilibrate
# print('Performing the energy minimization...')
# simulation.minimizeEnergy()
simulation.step(meltSteps)
stat_melt = simulation.context.getState(getPositions=True)
position_melt = stat_melt.getPositions()

#####################
#  Annealing
#####################
integrator = MTSLangevinIntegrator(annealT1, freq, timestep, [(0,4), (1,1), (2,1)])
simulation = Simulation(topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(position_melt)
simulation.step(annealSteps)
stat_anneal1 = simulation.context.getState(getPositions=True)
position_anneal1 = stat_anneal1.getPositions()
print("Annealing to 330K...")
integrator = MTSLangevinIntegrator(annealT2, freq, timestep, [(0,4), (1,1), (2,1)])
simulation = Simulation(topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(position_anneal1)
simulation.step(annealSteps)
stat_anneal2 = simulation.context.getState(getPositions=True)
position_anneal2 = stat_anneal2.getPositions()
print("Annealing to 298K...")
integrator = MTSLangevinIntegrator(temperature, freq, timestep, [(0,4), (1,1), (2,1)])
simulation = Simulation(topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(position_anneal2)
simulation.step(annealSteps)
stat_anneal3 = simulation.context.getState(getPositions=True)
position_anneal3 = stat_anneal3.getPositions()
PDBFile.writeFile(simulation.topology, position_anneal3, open(output+'anneal_'+input_name+'.pdb', 'w'))

#####################
#  Image Integrator
#####################
imageInteg = ImageMTSLangevinIntegrator(temperature, freq, timestep, [(0,4), (1,1), (2,1)])
# add image charge particles to the topology and position
print("Building image system...")
atoms = topology.atoms()
newChain = topology.addChain()
newResidue = topology.addResidue('IM', newChain)
print('The number of chains of new topology is: '+str(topology.getNumChains()))
imsig = 1*nanometer
imeps = 0*kilojoule_per_mole
print('number of real atoms is ' + str(system.getNumParticles()))
# create external force object if the external electric field by the walls is not zero
# (either dur to potential difference or charged wall)
# not for the amoeba force field because not consider dipoles and quadrupoles!
if constEfield._value > 0:
    constVForce = CustomExternalForce('efield*q*z')
    constVForce.addGlobalParameter('efield', constEfield._value)
    constVForce.addPerParticleParameter('q')

# add image particles for all real particles regardless of the charge
for i in range(nRealAtoms):
    (q, dip, quad, axisType, atomZ, atomX, atomY, thole, dampFactor, polarity) = mtpForce.getMultipoleParameters(i)
    newAtom = topology.addAtom('IM', next(atoms).element, newResidue)
    pos = position_anneal3[i].value_in_unit(nanometer)
    if q != -q:
        # position the image charge particles at the mirror location with respect to the left wall at z=0
        position_anneal3.append((pos[0], pos[1], -pos[2])*nanometer)
    else:
        # shift the position of image particle for non-charged to avoid zero distance divergence
        # between the real an image particles
        # (This is required for the wall atoms)
        position_anneal3.append((pos[0], pos[1], -pos[2]+0.001)*nanometer)
    idxat = system.addParticle(0*dalton)
    dip1 = Quantity((-dip[0], -dip[1], dip[2])).in_units_of(nanometer*elementary_charge)
    quad1 = Quantity(tuple(-d for d in quad)).in_units_of(nanometer**2*elementary_charge)
    idxat2 = mtpForce.addMultipole(-q, dip1, quad1, axisType, atomZ+nRealAtoms, atomX+nRealAtoms, atomY+nRealAtoms,
                                   thole, dampFactor, polarity)
    # print("idxat: "+ str(idxat) + ", idxat2: "+ str(idxat2))
    # add image pairs
    imageInteg.setImagePair(idxat, i)
    # add fake bond between image and parent so that they are always in the same periodic cell
    hbForce.addBond(idxat, i, 0, 0)
    # add charged particles into constVForce to apply the E_field by the walls
    # (either potential difference or charged wall)
    if constEfield._value > 0:
        idxres = constVForce.addParticle(i, [q])

# print the image pairs
# print(imageInteg.getImagePairs())

# add the E-field external force to the system
if constEfield._value > 0:
    system.addForce(constVForce)

##########################
# Add exclusion for image particles
##########################
# keep the electrostatic part (q) but ignore the vdW interaction
print('number of all atoms is '+str(system.getNumParticles()))
# print('number of nonbonded force particles is '+str(nbforce.getNumParticles()))

# There are no electrostatic interactions between real atoms within 1-2 and 1-3.
# There are no vdw interactions between image particles.
# Need to set covalent maps for image particles to construct the exception between 1-2 mpoles.
for i in range(nRealAtoms):
    covalent = mtpForce.getCovalentMaps(i)
    # print(covalent)
    for j in range(len(covalent)):
        mtpForce.setCovalentMap(nRealAtoms+i, j, covalent[j])
    # print(mtpForce.getCovalentMaps(nRealAtoms+i))

imgType = vdwForce.addParticleType(imsig, imeps)
# print(vdwForce.getNumParticleTypes())
for i in range(nRealAtoms):
    # exclusions = vdwForce.getParticleExclusions(i)
    # para = vdwForce.getParticleParameters(i)
    # print(para[5])
    # print(vdwForce.getParticleTypeParameters(para[5]))  # bug in OpenMM7.7
    img = vdwForce.addParticle(i+nRealAtoms, imgType, 1.0)

#########################
# Create simulation object
#########################
simulation = Simulation(topology, system, imageInteg, platform, platformProperties)
simulation.context.setPositions(position_anneal3)
# write the pdb file before beginning
print('The kinetic energy after annealing is '+str(stat_anneal3.getKineticEnergy()))
print('The potential energy after annealing is '+str(stat_anneal3.getPotentialEnergy()))
# minimize and equilibrate
# print('Performing the energy minimization...')
# simulation.minimizeEnergy()

#############################
# Actual simulation routine
#############################
print('Equilibrating...')
simulation.step(equilibrationSteps)
state_eq = simulation.context.getState(getEnergy=True, getForces=True, getPositions=True, enforcePeriodicBox=True)
position_eq = state_eq.getPositions()
PDBFile.writeFile(simulation.topology, position_eq, open(output+'eq_'+input_name+'.pdb', 'w'))

print('Simulating and sampling...')
simulation.reporters.append(DCDReporter(output+'samp_' +input_name+'.dcd', 2000))
simulation.reporters.append(StateDataReporter(output+'samp_'+input_name+'.log', 2000, totalSteps=samplingSteps,
                                              step=True, potentialEnergy=True, kineticEnergy=True,
                                              totalEnergy=True, temperature=True, speed=True, remainingTime=True,
                                              density=True, separator='\t'))
simulation.reporters.append(CheckpointReporter(output+'samp_'+input_name+'.chk', 50000))


# simulation.reporters[1] is the context of StateDataReporter(). Add the initial state before the equilibration.
simulation.reporters[1].report(simulation, state_eq)

#simulation.step(samplingSteps)

# log file to report the energy components
with open(output+'samp_'+input_name+'_ener.log', 'w') as enerlog:
    # write the header for the energy log file
    enerlog.write('# Energy log file\n')
    enerlog.write('# x1 : time (ps)\n')
    for j in range(system.getNumForces()):
        f = system.getForce(j)
        enerlog.write('# x'+str(j+2) + ' : ' + str(type(f)) + ' (kJ/mol) '+'forceGroup: '+str(f.getForceGroup())+'\n')

    with open(output+'samp'+input_name+'_dipole.log', 'w') as dipolog:
        dipolog.write("Induced dipoles of each particles.\n")

        for i in range(1, int(samplingSteps/2000)+1):
            simulation.step(2000)
            enerlog.write(str(i*2000))
            for j in range(system.getNumForces()):
                f = system.getForce(j)
                enerlog.write('  ' + str(simulation.context.getState(getEnergy=True, groups=2**j).getPotentialEnergy().value_in_unit(kilojoule_per_mole)))
            enerlog.write('\n')
            enerlog.flush()  # 刷新缓冲区，即将缓冲区中的数据立刻写入文件，同时清空缓冲区

            dipolog.write(str(i*2000)+'\n')
            for k in mtpForce.getInducedDipoles(simulation.context):
                dipolog.write(str(k) + '\n')
            dipolog.write("Total dipoles: "+str(mtpForce.getTotalDipoles(simulation.context))+'\n')
            dipolog.flush()


state_final = simulation.context.getState(getEnergy=True, getForces=True, getPositions=True,
                                    enforcePeriodicBox=True)
position_final = state_final.getPositions()
# boxVectors = state.getPeriodicBoxVectors()
PDBFile.writeFile(simulation.topology, position_final, open(output+'samp_'+input_name+'.pdb', 'w'))


time_end = time.time()
end_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time_end))
print('Done!')
print('Running date: ', start_time, '-', end_time)
print('Running Time (min): ', (time_end-time_start)/60)

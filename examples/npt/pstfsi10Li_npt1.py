from openmm.app import *
from openmm import *
from openmm.unit import *
import time
import sys
from mcbarostate import *
import numpy as np
import itertools

time_start = time.time()
start_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time_start))

# print(sys.getrecursionlimit())
sys.setrecursionlimit(2500)

# integration parameters
# pressure = Vec3(1.0, 1.0, 1.0) * atmospheres
pressure = 1*atmospheres
temperature = 298 * kelvin
freq = 1.0 / picosecond
timestep = 1 * femtosecond
barostatInterval = 25
cutoff = 1.0 * nanometers
vdw_cutoff = 1.2 * nanometers
ewald_error_tolerance = 0.0005
polarization = 'mutual'
mutualInducedTargetEpsilon = 0.0005
inducedMaxIteration = 100
equilibrationSteps = 6000000
zmin = 0
zmax = 5.0

# force field
ffField = 'pstfsi10_pAlt.xml'
forceField = ForceField(ffField)
imsig = 1*nanometer
imeps = 0*kilojoule/mole

# input files and output path
input_name = 'pstfsi10Li-h2o_npt1'
pdb = PDBFile('2pstfsi10Li-H2O_39_42.pdb')
output = 'output/'
# simulation options
print('Building system...')
platform = Platform.getPlatformByName('CUDA')
platformProperties = {'Precision': 'mixed', 'DeviceIndex': '0'}
# if there is no connect info of graphene residue, need to creat a file to store the bond info of graphene
# pdb.topology.loadBondDefinitions('GrapheneConnectivity.xml')
# pdb.topology.createStandardBonds()

topology = pdb.topology
position = pdb.positions
boxvec = pdb.topology.getUnitCellDimensions()._value
if (boxvec[2]-zmax*2) != 0:  #_value remove the unit
    topology.setUnitCellDimensions(Vec3(boxvec[0], boxvec[1], zmax*2)*nanometer)

system = forceField.createSystem(topology, nonbondedMethod=PME, nonbondedCutoff=cutoff, vdwCutoff=vdw_cutoff,
                                 ewaldErrorTolerance=ewald_error_tolerance, polarization=polarization,
                                 mutualInducedTargetEpsilon=mutualInducedTargetEpsilon, 
                                 mutualInducedMaxIterations=inducedMaxIteration,
                                 constraints=HBonds, rigidWater=True)
#system.addForce(MonteCarloBarostat(pressure, temperature, barostatInterval))

for f in system.getForces():
        if (isinstance(f, AmoebaMultipoleForce)
            or isinstance(f, AmoebaVdwForce)
            or isinstance(f, AmoebaGeneralizedKirkwoodForce)
            or isinstance(f, AmoebaWcaDispersionForce)):
            f.setForceGroup(1)

integrator = MTSLangevinIntegrator(temperature, freq, timestep, [(0,4), (1,1)])
nRealAtoms = system.getNumParticles()
#for f in [system.getForce(i) for i in range(system.getNumForces())]:
#    print(type(f))
vdwForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == AmoebaVdwForce][0]
numvdwPair = vdwForce.getNumTypePairs()
numvdwParticle = vdwForce.getNumParticles()
vdwMethod = vdwForce.getNonbondedMethod()  # 1 means CutoffPeriodic, 0 means NoCutoff.
vdwFuc = vdwForce.getPotentialFunction()
# nbforce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == NonbondedForce][0]
multipoleForce = [f for f in [system.getForce(i) for i in range(system.getNumForces())] if type(f) == AmoebaMultipoleForce][0]
print(multipoleForce.getMutualInducedMaxIterations())

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

# set force group to report each energy components
for i in range(system.getNumForces()):
    f = system.getForce(i)
#    f.setForceGroup(i)
    print(str(f.getForceGroup()) + str(type(f)))


simulation = Simulation(topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(position)
barostat = Barostat(simulation, pressure, temperature, barostatInterval)

# minimize and equilibrate
print('Performing the energy minimization...')
simulation.minimizeEnergy()
state_mini = simulation.context.getState(getEnergy=True,getForces=True,getVelocities=True,getPositions=True)
position_mini = state_mini.getPositions()
# PDBFile.writeFile(topology, position_mini, open(output+miniPDB+'.pdb', 'w+'))

simulation.reporters.append(PDBReporter(output+'trj_'+input_name+'.pdb', 10000))
simulation.reporters.append(StateDataReporter(output+'eq_'+input_name+'.log', 10000, totalSteps=equilibrationSteps, step=True,
                                              speed=True, progress=True, potentialEnergy=True, temperature=True,
                                              volume=True, remainingTime=True, density=True, separator='\t'))
print('Equilibrating...')
barostat.step_poly(equilibrationSteps)

state = simulation.context.getState(getEnergy=True, getForces=True, getPositions=True,
                                    enforcePeriodicBox=True)
position = state.getPositions()
boxVectors = state.getPeriodicBoxVectors()
PDBFile.writeFile(simulation.topology, position, open(output+'eq_'+input_name+'.pdb', 'w'))
#with open(output+'eq_'+input_name+'.log', 'a') as f:
#    f.write("The periodic box vectors are " + str(boxVectors)+'\n')
#    f.write("The distance between two graphene is "+str(barostat.aveDist))

time_end = time.time()
end_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(time_end))
print('Done!')
print('Running date: ', start_time, '-', end_time)
print('Running Time (min): ', (time_end-time_start)/60)

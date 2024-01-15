import openmm.unit as unit
import openmm as mm
import random
import math
import numpy as np
from cvforce import *


class Barostat(object):

    def __init__(self, simeq, pressure = 1.0*unit.bar, temperature = 298.0 * unit.kelvin, barofreq = 25):
        self.simeq = simeq
        self.temperature = temperature
        self.barofreq = barofreq
        self.RT = unit.BOLTZMANN_CONSTANT_kB*temperature*unit.AVOGADRO_CONSTANT_NA
        self.naccept = 0
        self.ntrials = 0
        self.celldim = self.simeq.topology.getUnitCellDimensions()
        self.lenscale = self.celldim[2]*0.01
        self.pressure = pressure*self.celldim[0]*self.celldim[1]*unit.AVOGADRO_CONSTANT_NA
        self.numres = simeq.topology.getNumResidues()
        self.numRealRes = simeq.topology.getNumResidues()
        self.firstRealRes = 0
        self.resmass = np.zeros(self.numres)*unit.dalton
        self.resNumAtoms = np.zeros(self.numres,int)
        self.resFirstAtomIdxs = np.zeros(self.numres,int)
        self.firstidx = 0
        self.lastidx = simeq.topology.getNumAtoms()
        self.aveDist = 0.0
        for at in simeq.topology.atoms():
            self.resmass[at.residue.index] += simeq.system.getParticleMass(at.index)
            self.resNumAtoms[at.residue.index] += 1
        idxAtom = 0
        for i in range(self.numres):
            self.resFirstAtomIdxs[i] = idxAtom
            # Exclude fixed residues from barostat acceptance criteria
            if self.resmass[i] == 0*unit.dalton:
                self.numRealRes -= 1
                if self.firstidx == idxAtom:
                    self.firstidx += self.resNumAtoms[i]
                    self.firstRealRes += 1
                elif self.lastidx == simeq.topology.getNumAtoms():
                    self.lastidx = idxAtom
            idxAtom += self.resNumAtoms[i]
        self.numVirtRes = self.numres - self.numRealRes

    def getAcceptRatio(self):
        return self.naccept/self.ntrials

    def metropolis(self,pecomp):
        if pecomp < 0*self.RT:
            return True
        elif (random.uniform(0.0,1.0) < np.exp(-(pecomp)/self.RT)):
            return True
        else:
            return False

    def metropolis_modify(self,pecomp):
        if pecomp < 0*self.RT:
            return True
        elif (random.uniform(0.4,1.0) < np.exp(-(pecomp)/self.RT)):
            return True
        else:
            return False

    def step_poly(self,nstep):
        niter = int(nstep/self.barofreq)
        numSample = 0
        for istep in range(niter):
            #cvforce.SlabCorrection(self.simeq.context.getSystem(),
            #                       self.simeq.context.getState(getPositions=True).getPeriodicBoxVolume().value_in_unit((unit.nanometer)**3))
            self.simeq.step(self.barofreq)
            self.ntrials += 1
            statebak = self.simeq.context.getState(getEnergy=True,getPositions=True)
            oldE = statebak.getPotentialEnergy()
            oldpos = statebak.getPositions()
            newpos0 = np.asarray(oldpos.value_in_unit(unit.nanometer))
            newpos = np.asarray(oldpos.value_in_unit(unit.nanometer))
            boxvec = statebak.getPeriodicBoxVectors()
            oldboxlen = boxvec[2][2]
            deltalen = self.lenscale*(random.uniform(0,1)*2.-1)
            newboxlen = oldboxlen+deltalen
            #idvirtres = 0
            for i in range(self.numres):
                fidx = self.resFirstAtomIdxs[i]
                lidx = self.resFirstAtomIdxs[i]+self.resNumAtoms[i]
                if self.resmass[i] == 0*unit.dalton:
                    if newpos0[fidx,2] > oldboxlen/unit.nanometer/4:
                        newpos[fidx:lidx,2] += 0.5*deltalen/unit.nanometer
                else:
                    centerZ = 0.0
                    for iatom in range(self.resNumAtoms[i]):
                        centerZ += newpos0[fidx+iatom,2]
                    centerZ = centerZ/self.resNumAtoms[i]
                    centerZ_new = centerZ*(newboxlen/oldboxlen)
                    posdel = centerZ_new - centerZ
                    newpos[fidx:lidx,2] += posdel
            self.simeq.context.setPositions(newpos)
            self.simeq.context.setPeriodicBoxVectors(boxvec[0],boxvec[1],mm.Vec3(0,0,newboxlen/unit.nanometer)*unit.nanometer)
            statenew = self.simeq.context.getState(getEnergy=True,getPositions=True)
            newE = statenew.getPotentialEnergy()
            energy = newE.value_in_unit(unit.kilojoules_per_mole)
            if math.isnan(energy):
                print("Energy is Nan!")
            elif math.isinf(energy):
                print("Energy is infinite!")
            w = newE-oldE + self.pressure*0.5*deltalen - self.numRealRes * self.RT*np.log(newpos[-1,2]/newpos0[-1,2])
            #print(str(newpos[-1,2]-newpos0[-1,2])+' '+str(0.5*deltalen))
            if self.metropolis(w):
                self.naccept += 1
                if istep > 160000:
                    self.aveDist += newpos[-1,2]
                    numSample += 1
                print('w is '+str(w)+', newE: '+str(newE)+', oldE: '+str(oldE)+", dE: "+str(newE-oldE))
                print(str(self.getAcceptRatio())+', '+str(newboxlen))
                print("the z coord of moved C is "+str(newpos[-1,2])+'\n')
            else:
                self.simeq.context.setPositions(oldpos)
                self.simeq.context.setPeriodicBoxVectors(boxvec[0],boxvec[1],mm.Vec3(0,0,oldboxlen/unit.nanometer)*unit.nanometer)
            if self.ntrials >= 10:
                if (self.naccept < 0.25*self.ntrials) :
                    self.lenscale /= 1.1
                    self.ntrials = 0
                    self.naccept = 0
                    print("Lenscale is too large, accept ratio is smaller than 0.25, reset lenscale "+str(self.lenscale))
                elif self.naccept > 0.75*self.ntrials :
                    self.lenscale = min(self.lenscale*1.1, newboxlen*0.3)
                    self.ntrials = 0
                    self.naccept = 0
                    print("Lenscale is too small, accept ratio is smaller than 0.75, reset lenscale "+str(self.lenscale))
        # self.aveDist = self.aveDist/numSample
        # print("The average Distance between two graphene is "+str(self.aveDist))

    def step_nogra(self, nstep):
        niter = int(nstep / self.barofreq)
        numSample = 0
        for istep in range(niter):
            # cvforce.SlabCorrection(self.simeq.context.getSystem(),
            #                       self.simeq.context.getState(getPositions=True).getPeriodicBoxVolume().value_in_unit((unit.nanometer)**3))
            self.simeq.step(self.barofreq)
            self.ntrials += 1
            statebak = self.simeq.context.getState(getEnergy=True, getPositions=True)
            oldE = statebak.getPotentialEnergy()
            oldpos = statebak.getPositions()
            newpos0 = np.asarray(oldpos.value_in_unit(unit.nanometer))
            newpos = np.asarray(oldpos.value_in_unit(unit.nanometer))
            boxvec = statebak.getPeriodicBoxVectors()
            oldboxlen = boxvec[2][2]
            deltalen = self.lenscale * (random.uniform(0, 1) * 2. - 1)
            newboxlen = oldboxlen + deltalen
            # idvirtres = 0
            for i in range(self.numres):
                fidx = self.resFirstAtomIdxs[i]
                lidx = self.resFirstAtomIdxs[i] + self.resNumAtoms[i]
                if self.resmass[i] == 0 * unit.dalton:
                    if newpos0[fidx, 2] > oldboxlen / unit.nanometer / 4:
                        newpos[fidx:lidx, 2] += 0.5 * deltalen / unit.nanometer
                else:
                    centerZ = 0.0
                    for iatom in range(self.resNumAtoms[i]):
                        centerZ += newpos0[fidx + iatom, 2]
                    centerZ = centerZ / self.resNumAtoms[i]
                    centerZ_new = centerZ * (newboxlen / oldboxlen)
                    posdel = centerZ_new - centerZ
                    newpos[fidx:lidx, 2] += posdel
            self.simeq.context.setPositions(newpos)
            self.simeq.context.setPeriodicBoxVectors(boxvec[0], boxvec[1],
                                                     mm.Vec3(0, 0, newboxlen / unit.nanometer) * unit.nanometer)
            statenew = self.simeq.context.getState(getEnergy=True, getPositions=True)
            newE = statenew.getPotentialEnergy()
            w = newE - oldE + self.pressure * 0.5 * deltalen - self.numRealRes * self.RT * np.log(newboxlen/oldboxlen)
            # print(str(newpos[-1,2]-newpos0[-1,2])+' '+str(0.5*deltalen))
            if self.metropolis(w):
                self.naccept += 1
                if istep > 160000:
                    self.aveDist += newpos[-1, 2]
                    numSample += 1
                print('w is ' + str(w) + ', newE: ' + str(newE) + ', oldE: ' + str(oldE) + ", dE: " + str(
                    newE - oldE))
                print(str(self.getAcceptRatio()) + ', ' + str(newboxlen))
                # print("the z coord of moved C is " + str(newpos[-1, 2]) + '\n')
            else:
                self.simeq.context.setPositions(oldpos)
                self.simeq.context.setPeriodicBoxVectors(boxvec[0], boxvec[1],
                                                         mm.Vec3(0, 0, oldboxlen / unit.nanometer) * unit.nanometer)
            if self.ntrials >= 10:
                if (self.naccept < 0.25 * self.ntrials):
                    self.lenscale /= 1.1
                    self.ntrials = 0
                    self.naccept = 0
                    print("Lenscale is too large, accept ratio is smaller than 0.25, reset lenscale " + str(
                        self.lenscale))
                elif self.naccept > 0.75 * self.ntrials:
                    self.lenscale = min(self.lenscale * 1.1, newboxlen * 0.3)
                    self.ntrials = 0
                    self.naccept = 0
                    print("Lenscale is too small, accept ratio is smaller than 0.75, reset lenscale " + str(
                        self.lenscale))
        # self.aveDist = self.aveDist/numSample
        # print("The average Distance between two graphene is "+str(self.aveDist))
    def step_1(self,nstep, istep=1):
        if nstep%self.barofreq == 0:
            # cvforce.SlabCorrection(self.simeq.context.getSystem(),
            #                       self.simeq.context.getState(getPositions=True).getPeriodicBoxVolume().value_in_unit((unit.nanometer)**3))
            # self.simeq.step(self.barofreq)
            self.ntrials += 1
            statebak = self.simeq.context.getState(getEnergy=True, getPositions=True)
            oldE = statebak.getPotentialEnergy()
            oldpos = statebak.getPositions()
            newpos0 = np.asarray(oldpos.value_in_unit(unit.nanometer))
            newpos = np.asarray(oldpos.value_in_unit(unit.nanometer))
            boxvec = statebak.getPeriodicBoxVectors()
            oldboxlen = boxvec[2][2]
            deltalen = self.lenscale * (random.uniform(0, 1) * 2. - 1)
            newboxlen = oldboxlen + deltalen
            # idvirtres = 0
            for i in range(self.numres):
                fidx = self.resFirstAtomIdxs[i]
                lidx = self.resFirstAtomIdxs[i] + self.resNumAtoms[i]
                if self.resmass[i] == 0 * unit.dalton:
                    if newpos0[fidx, 2] > oldboxlen / unit.nanometer / 4:
                        newpos[fidx:lidx, 2] += 0.5 * deltalen / unit.nanometer
                else:
                    centerZ = 0.0
                    for iatom in range(self.resNumAtoms[i]):
                        centerZ += newpos0[fidx + iatom, 2]
                    centerZ = centerZ / self.resNumAtoms[i]
                    centerZ_new = centerZ * (newboxlen / oldboxlen)
                    posdel = centerZ_new - centerZ
                    newpos[fidx:lidx, 2] += posdel
            self.simeq.context.setPositions(newpos)
            self.simeq.context.setPeriodicBoxVectors(boxvec[0], boxvec[1],
                                                     mm.Vec3(0, 0, newboxlen / unit.nanometer) * unit.nanometer)
            statenew = self.simeq.context.getState(getEnergy=True, getPositions=True)
            newE = statenew.getPotentialEnergy()
            w = newE - oldE + self.pressure * 0.5 * deltalen - self.numRealRes * self.RT * np.log(
                newpos[-1, 2] / newpos0[-1, 2])
            # print(str(newpos[-1,2]-newpos0[-1,2])+' '+str(0.5*deltalen))
            if self.metropolis(w):
                self.naccept += 1
                # if istep > 160000:
                #     self.aveDist += newpos[-1, 2]
                    # numSample += 1
                print('w is ' + str(w) + ', newE: ' + str(newE) + ', oldE: ' + str(oldE) + ", dE: " + str(newE - oldE))
                print(str(self.getAcceptRatio()) + ', ' + str(newboxlen))
                print("the z coord of moved C is " + str(newpos[-1, 2]) + '\n')
            else:
                self.simeq.context.setPositions(oldpos)
                self.simeq.context.setPeriodicBoxVectors(boxvec[0], boxvec[1],
                                                         mm.Vec3(0, 0, oldboxlen / unit.nanometer) * unit.nanometer)
            if self.ntrials >= 10:
                if (self.naccept < 0.25 * self.ntrials):
                    self.lenscale /= 1.1
                    self.ntrials = 0
                    self.naccept = 0
                    print("Lenscale is too large, accept ratio is smaller than 0.25, reset lenscale " + str(
                        self.lenscale))
                elif self.naccept > 0.75 * self.ntrials:
                    self.lenscale = min(self.lenscale * 1.1, newboxlen * 0.3)
                    self.ntrials = 0
                    self.naccept = 0
                    print("Lenscale is too small, accept ratio is smaller than 0.75, reset lenscale " + str(
                        self.lenscale))
        # self.aveDist = self.aveDist/numSample
        # print("The average Distance between two graphene is "+str(self.aveDist))
        else:
            self.simeq.step(istep)

    # modified metropolis
    def step_poly_mod(self,nstep):
        niter = int(nstep/self.barofreq)
        numSample = 0
        for istep in range(niter):
            #cvforce.SlabCorrection(self.simeq.context.getSystem(),
            #                       self.simeq.context.getState(getPositions=True).getPeriodicBoxVolume().value_in_unit((unit.nanometer)**3))
            self.simeq.step(self.barofreq)
            self.ntrials += 1
            statebak = self.simeq.context.getState(getEnergy=True,getPositions=True)
            oldE = statebak.getPotentialEnergy()
            oldpos = statebak.getPositions()
            newpos0 = np.asarray(oldpos.value_in_unit(unit.nanometer))
            newpos = np.asarray(oldpos.value_in_unit(unit.nanometer))
            boxvec = statebak.getPeriodicBoxVectors()
            oldboxlen = boxvec[2][2]
            deltalen = self.lenscale*(random.uniform(0,1)*2.-1)
            newboxlen = oldboxlen+deltalen
            #idvirtres = 0
            for i in range(self.numres):
                fidx = self.resFirstAtomIdxs[i]
                lidx = self.resFirstAtomIdxs[i]+self.resNumAtoms[i]
                if self.resmass[i] == 0*unit.dalton:
                    if newpos0[fidx,2] > oldboxlen/unit.nanometer/4:
                        newpos[fidx:lidx,2] += 0.5*deltalen/unit.nanometer
                else:
                    centerZ = 0.0
                    for iatom in range(self.resNumAtoms[i]):
                        centerZ += newpos0[fidx+iatom,2]
                    centerZ = centerZ/self.resNumAtoms[i]
                    centerZ_new = centerZ*(newboxlen/oldboxlen)
                    posdel = centerZ_new - centerZ
                    newpos[fidx:lidx,2] += posdel
            self.simeq.context.setPositions(newpos)
            self.simeq.context.setPeriodicBoxVectors(boxvec[0],boxvec[1],mm.Vec3(0,0,newboxlen/unit.nanometer)*unit.nanometer)
            statenew = self.simeq.context.getState(getEnergy=True,getPositions=True)
            newE = statenew.getPotentialEnergy()
            energy = newE.value_in_unit(unit.kilojoules_per_mole)
            if math.isnan(energy):
                print("Energy is Nan!")
            elif math.isinf(energy):
                print("Energy is infinite!")
            w = newE-oldE + self.pressure*0.5*deltalen - self.numRealRes * self.RT*np.log(newpos[-1,2]/newpos0[-1,2])
            #print(str(newpos[-1,2]-newpos0[-1,2])+' '+str(0.5*deltalen))
            if self.metropolis_modify(w):
                self.naccept += 1
                if istep > 160000:
                    self.aveDist += newpos[-1,2]
                    numSample += 1
                print('w is '+str(w)+', newE: '+str(newE)+', oldE: '+str(oldE)+", dE: "+str(newE-oldE))
                print(str(self.getAcceptRatio())+', '+str(newboxlen))
                print("the z coord of moved C is "+str(newpos[-1,2])+'\n')
            else:
                self.simeq.context.setPositions(oldpos)
                self.simeq.context.setPeriodicBoxVectors(boxvec[0],boxvec[1],mm.Vec3(0,0,oldboxlen/unit.nanometer)*unit.nanometer)
            if self.ntrials >= 10:
                if (self.naccept < 0.25*self.ntrials) :
                    self.lenscale /= 1.1
                    self.ntrials = 0
                    self.naccept = 0
                    print("Lenscale is too large, accept ratio is smaller than 0.25, reset lenscale "+str(self.lenscale))
                elif self.naccept > 0.75*self.ntrials :
                    self.lenscale = min(self.lenscale*1.1, newboxlen*0.3)
                    self.ntrials = 0
                    self.naccept = 0
                    print("Lenscale is too small, accept ratio is smaller than 0.75, reset lenscale "+str(self.lenscale))
        # self.aveDist = self.aveDist/numSample
        # print("The average Distance between two graphene is "+str(self.aveDist))

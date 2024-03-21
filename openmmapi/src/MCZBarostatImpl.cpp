/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010-2020 Stanford University and the Authors.      *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "openmm/internal/MCZBarostatImpl.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/internal/OSRngSeed.h"
#include "openmm/Context.h"
#include "openmm/ImageKernels.h"
#include "openmm/OpenMMException.h"
#include "openmm/reference/SimTKOpenMMUtilities.h"
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>

using namespace OpenMM;
using namespace std;

MCZBarostatImpl::MCZBarostatImpl(const MCZBarostat& owner) : owner(owner), step(0) {
}

void MCZBarostatImpl::initialize(ContextImpl& context) {
    if (!context.getSystem().usesPeriodicBoundaryConditions())
        throw OpenMMException("A barostat cannot be used with a non-periodic system");
    kernel = context.getPlatform().createKernel(ApplyMCZBarostatKernel::Name(), context);
    kernel.getAs<ApplyMCZBarostatKernel>().initialize(context.getSystem(), owner);
    Vec3 box[3];
    context.getPeriodicBoxVectors(box[0], box[1], box[2]);
    double volume = box[0][0]*box[1][1]*box[2][2];
    lenScale = 0.01*box[2][2];
    numAttempted = 0;
    numAccepted = 0;
    SimTKOpenMMUtilities::setRandomNumberSeed(owner.getRandomNumberSeed());

    numRealRes = context.getMolecules().size();
    numVirtAtoms = 0;
    for(int i=0; i<context.getSystem().getNumParticles(); i++){
        if(context.getSystem().getParticleMass(i) == 0.0){
            numRealRes--;
            numVirtAtoms++;
        }
    }
    cout<<"The number of real residues is "<<numRealRes<<endl;
}

void MCZBarostatImpl::updateContextState(ContextImpl& context, bool& forcesInvalid) {
    if (++step < owner.getFrequency() || owner.getFrequency() == 0)
        return;
    step = 0;

    // Compute the current potential energy.

    int groups = context.getIntegrator().getIntegrationForceGroups();
    double initialEnergy = context.getOwner().getState(State::Energy, false, groups).getPotentialEnergy();

    // Modify the periodic box size.
    //cout<<"Performing MCZBarostat..."<<endl;
    Vec3 box[3];
    context.getPeriodicBoxVectors(box[0], box[1], box[2]);
    //double volume = box[0][0]*box[1][1]*box[2][2];
    double oldzlen = box[2][2];
    double deltalen = lenScale*2*(SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber()-0.5);
    double newzlen = oldzlen+deltalen;
    //double newVolume = volume+deltaVolume;
    double lengthScale = newzlen/oldzlen;
    kernel.getAs<ApplyMCZBarostatKernel>().scaleCoordinates(context, lengthScale);
    //cout<<"apply mczbarostat done!"<<endl;
    context.getOwner().setPeriodicBoxVectors(box[0], box[1], box[2]*lengthScale);

    // Compute the energy of the modified system.
    
    double finalEnergy = context.getOwner().getState(State::Energy, false, groups).getPotentialEnergy();
    double pressure = context.getParameter(MCZBarostat::Pressure())*(AVOGADRO*1e-25); //bar*Na*nm^3, Pa*m^3/(1000*mol)=KJ/mol
    double kT = BOLTZ*context.getParameter(MCZBarostat::Temperature());
    double w = finalEnergy-initialEnergy + pressure*deltalen*box[0][0]*box[1][1] - numRealRes*kT*log(newzlen/oldzlen);
    //cout<<"compute w done!"<<endl;
    if (w > 0 && SimTKOpenMMUtilities::getUniformlyDistributedRandomNumber() > exp(-w/kT)) {
        // Reject the step.
        kernel.getAs<ApplyMCZBarostatKernel>().restoreCoordinates(context);
        context.getOwner().setPeriodicBoxVectors(box[0], box[1], box[2]);
    }
    else {
        numAccepted++;
        forcesInvalid = true;
        cout<<"Accepted ratio is "<<static_cast<float>(numAccepted)/(numAttempted+1)<<". The new box z length is "<<newzlen<<endl;
    }

    numAttempted++;
    if (numAttempted >= 10) {
        if (numAccepted < 0.25*numAttempted) {
            lenScale /= 1.1;
            numAttempted = 0;
            numAccepted = 0;
        }
        else if (numAccepted > 0.75*numAttempted) {
            lenScale = min(lenScale*1.1, newzlen*0.3);
            numAttempted = 0;
            numAccepted = 0;
        }
    }
}

map<string, double> MCZBarostatImpl::getDefaultParameters() {
    map<string, double> parameters;
    parameters[MCZBarostat::Pressure()] = getOwner().getDefaultPressure();
    parameters[MCZBarostat::Temperature()] = getOwner().getDefaultTemperature();
    return parameters;
}

vector<string> MCZBarostatImpl::getKernelNames() {
    vector<string> names;
    names.push_back(ApplyMCZBarostatKernel::Name());
    return names;
}


/* -------------------------------------------------------------------------- *
 *                               OpenMMImage                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c)                                                     *
 *         2023 ICCAS and the Authors.                                        *
 * Authors: Ruochao Wang                                                      *
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

#include "openmm/ImageLangevinIntegrator.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/ImageKernels.h"
#include <string>
#include <iostream>

using namespace OpenMM;
using std::string;
using std::vector;
using namespace std;

ImageLangevinIntegrator::ImageLangevinIntegrator(double temperature, double frictionCoeff, double stepSize, int numCells, double zmax) {
    setTemperature(temperature);
    setFriction(frictionCoeff);
    setStepSize(stepSize);
    setConstraintTolerance(1e-5);
    setRandomNumberSeed(0);
    setNumCells(numCells);
    setCellSize(zmax);
    setUseImageParticle(true);
}

// int ImageLangevinIntegrator::setImagePair(int image, int parent){
//     imageParticles.push_back(image);
//     imagePairs.emplace_back(std::make_pair(image, parent));
//     return imagePairs.size();
// }

void ImageLangevinIntegrator::initialize(ContextImpl& contextRef) {
    if(owner != NULL && &contextRef.getOwner() != owner)
        throw OpenMMException("This Integrator is already bound to a context!");

    context = &contextRef;
    owner = &contextRef.getOwner();

    ilkernel = context->getPlatform().createKernel(IntegrateImageLangevinStepKernel::Name(), contextRef);
    ilkernel.getAs<IntegrateImageLangevinStepKernel>().initialize(contextRef.getSystem(), *this);
    if(getUseImageParticle()){
        imgkernel = context->getPlatform().createKernel(ImageParticleKernel::Name(), contextRef);
        imgkernel.getAs<ImageParticleKernel>().initialize(contextRef.getSystem(), *this);
    }
}

void ImageLangevinIntegrator::cleanup() {
    ilkernel = Kernel();
    imgkernel = Kernel();
}

void ImageLangevinIntegrator::setTemperature(double temp) {
    if(temp < 0)
        throw OpenMMException("Temperature cannot be negative!");
        temperature = temp;
}

void ImageLangevinIntegrator::setFriction(double coeff) {
    if(coeff < 0)
        throw OpenMMException("Friction cannot be negative!");
    friction = coeff;
}

vector<string> ImageLangevinIntegrator::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(IntegrateImageLangevinStepKernel::Name());
    if(getUseImageParticle()){
        names.push_back(ImageParticleKernel::Name());
    }
    return names;
}

double ImageLangevinIntegrator::computeKineticEnergy() {
    return ilkernel.getAs<IntegrateImageLangevinStepKernel>().computeKineticEnergy(*context, *this);
}

void ImageLangevinIntegrator::step(int steps) {
    if(context == NULL)
        throw OpenMMException("This Integrator is not bound to a context!");

    for (int i=0; i<steps; i++)
    {
        //cout<<"Step "<<i<<endl;
        context->updateContextState();
        //cout<<"updateContextState done!"<<endl;
        context->calcForcesAndEnergy(true, false, getIntegrationForceGroups());
        //cout<<"calcForcesAndEnergy done!"<<endl;
        ilkernel.getAs<IntegrateImageLangevinStepKernel>().execute(*context, *this);
        //cout<<"IntegrateImageLangevinStepKernel excution done!"<<endl;
        if(getUseImageParticle()){
            imgkernel.getAs<ImageParticleKernel>().updateImagePositions(*context, *this);
            //cout<<"ImageParticleKernel excution done!"<<endl;
        }
    }
}
/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
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

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "openmm/internal/SlabCorrectionImpl.h"
#include "openmm/ImageKernels.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include <cmath>
#include <map>
#include <set>
#include <sstream>
#include <iostream>
#include <typeinfo>

using namespace OpenMM;
using namespace std;

SlabCorrectionImpl::SlabCorrectionImpl(const SlabCorrection& owner) : owner(owner) {
}

SlabCorrectionImpl::~SlabCorrectionImpl() {
}

void SlabCorrectionImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcSlabCorrectionKernel::Name(), context);
    kernel.getAs<CalcSlabCorrectionKernel>().initialize(context.getSystem(), owner);

    mtpImpl = NULL;
    //cout<<"The number of ForceImpls is "<<context.getForceImpls().size()<<endl;
    for(auto impl: context.getForceImpls()){
        //cout<<typeid(*impl).name()<<endl;
        mtpImpl = dynamic_cast<AmoebaMultipoleForceImpl*>(impl);
        if(mtpImpl != NULL){
            //mtpImpl->getSystemMultipoleMoments(context, multipoleMoments);
            break;
        }
    }
}

double SlabCorrectionImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    if ((groups&(1<<owner.getForceGroup())) != 0){
        // vector<double> multipoleMoments;
        // AmoebaMultipoleForceImpl* mtpImpl = NULL;
        // //cout<<"The number of ForceImpls is "<<context.getForceImpls().size()<<endl;
        // for(auto impl: context.getForceImpls()){
        //     //cout<<typeid(*impl).name()<<endl;
        //     mtpImpl = dynamic_cast<AmoebaMultipoleForceImpl*>(impl);
        //     if(mtpImpl != NULL){
        //         mtpImpl->getSystemMultipoleMoments(context, multipoleMoments);
        //         break;
        //     }
        // }
        if(mtpImpl == NULL)
            throw OpenMMException("The Context does not contain a AmoebaMultipoleForceImpl");
        else if(owner.useAmoebaDipole()){
            mtpImpl->getSystemMultipoleMoments(context, multipoleMoments); // unit is Debye
            double debye = 48.0321;
            return kernel.getAs<CalcSlabCorrectionKernel>().execute(context, includeForces, includeEnergy, multipoleMoments[3]/debye);
        }
        else
            return kernel.getAs<CalcSlabCorrectionKernel>().execute(context, includeForces, includeEnergy, 0.0);
    }  
    return 0.0;
}

std::vector<std::string> SlabCorrectionImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcSlabCorrectionKernel::Name());
    return names;
}


void SlabCorrectionImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcSlabCorrectionKernel>().copyParametersToContext(context, owner);
}

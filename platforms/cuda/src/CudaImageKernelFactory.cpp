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

#include <exception>

#include "CudaImageKernelFactory.h"
#include "CommonImageKernels.h"
#include "CudaContext.h"
#include "openmm/internal/windowsExport.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/OpenMMException.h"

using namespace OpenMM;

extern "C" OPENMM_EXPORT void registerPlatforms() {
}

extern "C" OPENMM_EXPORT void registerKernelFactories() {
    try {
        Platform& platform = Platform::getPlatformByName("CUDA");
        CudaImageKernelFactory* factory = new CudaImageKernelFactory();
        platform.registerKernelFactory(IntegrateImageLangevinStepKernel::Name(), factory);
        platform.registerKernelFactory(ImageParticleKernel::Name(), factory);
        platform.registerKernelFactory(IntegrateImageCustomStepKernel::Name(), factory);
        platform.registerKernelFactory(ApplyMCZBarostatKernel::Name(), factory);
        platform.registerKernelFactory(CalcSlabCorrectionKernel::Name(), factory);
    }
    catch (std::exception ex) {
        // Ignore
        std::cout<<"Failed to registerKernelFactory: "<<ex.what()<<std::endl;
    }
}

extern "C" OPENMM_EXPORT void registerImageCudaKernelFactories() {
    try {
        Platform::getPlatformByName("CUDA");
    }
    catch (...) {
        Platform::registerPlatform(new CudaPlatform());
    }
    registerKernelFactories();
}

KernelImpl* CudaImageKernelFactory::createKernelImpl(std::string name, const Platform& platform, ContextImpl& context) const {
    CudaContext& cu = *static_cast<CudaPlatform::PlatformData*>(context.getPlatformData())->contexts[0];
    if (name == IntegrateImageLangevinStepKernel::Name())
        return new CommonIntegrateImageLangevinStepKernel(name, platform, cu);
    if(name == ImageParticleKernel::Name())
        return new CommonImageParticleKernel(name, platform, cu);
    if (name == IntegrateImageCustomStepKernel::Name())
        return new CommonIntegrateImageCustomStepKernel(name, platform, cu);
    if (name == ApplyMCZBarostatKernel::Name())
        return new CommonApplyMCZBarostatKernel(name, platform, cu);
    if (name == CalcSlabCorrectionKernel::Name())
        return new CommonCalcSlabCorrectionKernel(name, platform, cu, context.getSystem());
    throw OpenMMException((std::string("Tried to create kernel with illegal kernel name '")+name+"'").c_str());
}

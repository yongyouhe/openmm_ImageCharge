#ifndef IMAGE_KERNELS_H_
#define IMAGE_KERNELS_H_

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
#include "openmm/ImageCustomIntegrator.h"
#include "openmm/MCZBarostat.h"
#include "openmm/Platform.h"
#include "openmm/System.h"
#include "openmm/Vec3.h"
#include <string>
#include <vector>

namespace OpenMM
{

    class IntegrateImageLangevinStepKernel : public KernelImpl {
    public:
        static std::string Name() {
            return "IntegrateImageLangevinStep";
        }
        IntegrateImageLangevinStepKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
        }
        /**
         * Initialize the kernel.
         * 
         * @param system     the System this kernel will be applied to
         * @param integrator the LangevinIntegrator this kernel will be used for
         */
        virtual void initialize(const System& system, const ImageLangevinIntegrator& integrator) = 0;
        /**
         * Execute the kernel.
         * 
         * @param context    the context in which to execute this kernel
         * @param integrator the LangevinIntegrator this kernel is being used for
         */
        virtual void execute(ContextImpl& context, const ImageLangevinIntegrator& integrator) = 0;
       /**
         * Compute the kinetic energy.
         * 
         * @param context    the context in which to execute this kernel
         * @param integrator the LangevinIntegrator this kernel is being used for
         */
        virtual double computeKineticEnergy(ContextImpl& context, const ImageLangevinIntegrator& integrator) = 0; 


    };

    /**
    * This kernel is invoked by ImageCustomIntegrator to take one time step.
    */
    class IntegrateImageCustomStepKernel : public KernelImpl {
    public:
        static std::string Name() {
            return "IntegrateImageCustomStep";
        }
        IntegrateImageCustomStepKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
        }
        /**
         * Initialize the kernel.
         * 
         * @param system     the System this kernel will be applied to
         * @param integrator the ImageCustomIntegrator this kernel will be used for
         */
        virtual void initialize(const System& system, const ImageCustomIntegrator& integrator) = 0;
        /**
         * Execute the kernel.
         * 
         * @param context    the context in which to execute this kernel
         * @param integrator the ImageCustomIntegrator this kernel is being used for
         * @param forcesAreValid if the context has been modified since the last time step, this will be
         *                       false to show that cached forces are invalid and must be recalculated.
         *                       On exit, this should specify whether the cached forces are valid at the
         *                       end of the step.
         */
        virtual void execute(ContextImpl& context, ImageCustomIntegrator& integrator, bool& forcesAreValid) = 0;
        /**
         * Compute the kinetic energy.
         * 
         * @param context    the context in which to execute this kernel
         * @param integrator the ImageCustomIntegrator this kernel is being used for
         * @param forcesAreValid if the context has been modified since the last time step, this will be
         *                       false to show that cached forces are invalid and must be recalculated.
         *                       On exit, this should specify whether the cached forces are valid at the
         *                       end of the step.
         */
        virtual double computeKineticEnergy(ContextImpl& context, ImageCustomIntegrator& integrator, bool& forcesAreValid) = 0;
        /**
         * Get the values of all global variables.
         *
         * @param context   the context in which to execute this kernel
         * @param values    on exit, this contains the values
         */
        virtual void getGlobalVariables(ContextImpl& context, std::vector<double>& values) const = 0;
        /**
         * Set the values of all global variables.
         *
         * @param context   the context in which to execute this kernel
         * @param values    a vector containing the values
         */
        virtual void setGlobalVariables(ContextImpl& context, const std::vector<double>& values) = 0;
        /**
         * Get the values of a per-DOF variable.
         *
         * @param context   the context in which to execute this kernel
         * @param variable  the index of the variable to get
         * @param values    on exit, this contains the values
         */
        virtual void getPerDofVariable(ContextImpl& context, int variable, std::vector<Vec3>& values) const = 0;
        /**
         * Set the values of a per-DOF variable.
         *
         * @param context   the context in which to execute this kernel
         * @param variable  the index of the variable to get
         * @param values    a vector containing the values
         */
        virtual void setPerDofVariable(ContextImpl& context, int variable, const std::vector<Vec3>& values) = 0;
    };

    class ApplyMCZBarostatKernel : public KernelImpl {
    public:
        static std::string Name() {
            return "ApplyMCZBarostat";
        }
        ApplyMCZBarostatKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
        }
        /**
         * Initialize the kernel.
         *
         * @param system          the System this kernel will be applied to
         * @param barostat        the MonteCarloBarostat this kernel will be used for
         * @param rigidMolecules  whether molecules should be kept rigid while scaling coordinates
         */
        virtual void initialize(const System& system, const Force& barostat, bool rigidMolecules=true) = 0;
        /**
         * Attempt a Monte Carlo step, scaling particle positions (or cluster centers) by a specified value.
         * This version scales z positions. Move the parallel plate first, then scale the positions of molecules with z axis.
         * This is called BEFORE the periodic box size is modified.  It should begin by translating each particle
         * or cluster into the first periodic box, so that coordinates will still be correct after the box size
         * is changed.
         *
         * @param context    the context in which to execute this kernel
         * @param scaleZ     the scale factor by which to multiply particle z-coordinate
         */
        virtual void scaleCoordinates(ContextImpl& context, double scaleZ) = 0;
        /**
         * Reject the most recent Monte Carlo step, restoring the particle positions to where they were before
         * scaleCoordinates() was last called.
         *
         * @param context    the context in which to execute this kernel
         */
        virtual void restoreCoordinates(ContextImpl& context) = 0;
    };

    class ImageParticleKernel : public KernelImpl {
    public:
        static std::string Name() {
            return "ImageParticle";
        }
        ImageParticleKernel(std::string name, const Platform& platform) : KernelImpl(name, platform) {
        }

        virtual void initialize(const System& system, const ImageIntegrator& integrator) = 0;
        /**
         * update positions of image partilces
         * @param context       the context in which to excute this kernel.
         * @param integrator    the ImageLangevinIntegrator this kernel is being used for.
        */
        virtual void updateImagePositions(ContextImpl& context, const ImageIntegrator& integrator) = 0;
    };

}// namespace OpenMM

#endif /*IMAGE_KERNELS_H_*/

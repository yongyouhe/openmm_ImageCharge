#ifndef COMMON_IMAGE_KERNELS_H_
#define COMMON_IMAGE_KERNELS_H_

/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014-2021 Stanford University and the Authors.      *
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

#include "openmm/ImageKernels.h"
#include "openmm/common/ComputeContext.h"
#include "openmm/common/ComputeArray.h"
#include "lepton/CompiledExpression.h"
#include "lepton/ExpressionProgram.h"
#include "openmm/internal/ImageCustomIntegratorUtilities.h"
#include "openmm/internal/CompiledExpressionSet.h"

namespace OpenMM {

/**
 * This kernel is invoked by ImageLangevinIntegrator to take one step.
 */
class CommonIntegrateImageLangevinStepKernel : public IntegrateImageLangevinStepKernel {
public:
    CommonIntegrateImageLangevinStepKernel(const std::string& name, const Platform& platform, ComputeContext& cc) :
            IntegrateImageLangevinStepKernel(name, platform), cc(cc), hasInitializedKernels(false) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system         the System this kernel will be applied to
     * @param integrator     the ImageLangevinIntegrator this kernel is being used for
     */
    void initialize(const System& system, const ImageLangevinIntegrator& integrator);
    /**
     * Execute the kernel.
     * @param constext      the context in which to execute this kernel.
     * @param integrator    the ImageLangevinIntegrator this kernel is being used for
     */
    void execute(ContextImpl& context, const ImageLangevinIntegrator& integrator);

    double computeKineticEnergy(ContextImpl& context, const ImageLangevinIntegrator& integrator);
private:
    ComputeContext& cc;
    double prevTemp, prevFriction, prevStepSize, zmax;
    bool hasInitializedKernels;
    ComputeArray params;
    ComputeKernel kernel1, kernel2;
};

/**
 * This kernel is invoked by ImageCustomIntegrator to take one time step.
 */
class CommonIntegrateImageCustomStepKernel : public IntegrateImageCustomStepKernel {
    
public:
    enum GlobalTargetType {DT, VARIABLE, PARAMETER};

    CommonIntegrateImageCustomStepKernel(std::string name, const Platform& platform, ComputeContext& cc) : IntegrateImageCustomStepKernel(name, platform), cc(cc),
            hasInitializedKernels(false), needsEnergyParamDerivs(false) {
    }
    /**
     * Initialize the kernel.
     * 
     * @param system     the System this kernel will be applied to
     * @param integrator the ImageCustomIntegrator this kernel will be used for
     */
    void initialize(const System& system, const ImageCustomIntegrator& integrator);
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
    void execute(ContextImpl& context, ImageCustomIntegrator& integrator, bool& forcesAreValid);
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
    double computeKineticEnergy(ContextImpl& context, ImageCustomIntegrator& integrator, bool& forcesAreValid);
    /**
     * Get the values of all global variables.
     *
     * @param context   the context in which to execute this kernel
     * @param values    on exit, this contains the values
     */
    void getGlobalVariables(ContextImpl& context, std::vector<double>& values) const;
    /**
     * Set the values of all global variables.
     *
     * @param context   the context in which to execute this kernel
     * @param values    a vector containing the values
     */
    void setGlobalVariables(ContextImpl& context, const std::vector<double>& values);
    /**
     * Get the values of a per-DOF variable.
     *
     * @param context   the context in which to execute this kernel
     * @param variable  the index of the variable to get
     * @param values    on exit, this contains the values
     */
    void getPerDofVariable(ContextImpl& context, int variable, std::vector<Vec3>& values) const;
    /**
     * Set the values of a per-DOF variable.
     *
     * @param context   the context in which to execute this kernel
     * @param variable  the index of the variable to get
     * @param values    a vector containing the values
     */
    void setPerDofVariable(ContextImpl& context, int variable, const std::vector<Vec3>& values);
private:
    class ReorderListener;
    class GlobalTarget;
    class DerivFunction;
    std::string createPerDofComputation(const std::string& variable, const Lepton::ParsedExpression& expr, ImageCustomIntegrator& integrator,
        const std::string& forceName, const std::string& energyName, std::vector<const TabulatedFunction*>& functions,
        std::vector<std::pair<std::string, std::string> >& functionNames);
    void prepareForComputation(ContextImpl& context, ImageCustomIntegrator& integrator, bool& forcesAreValid);
    Lepton::ExpressionTreeNode replaceDerivFunctions(const Lepton::ExpressionTreeNode& node, OpenMM::ContextImpl& context);
    void findExpressionsForDerivs(const Lepton::ExpressionTreeNode& node, std::vector<std::pair<Lepton::ExpressionTreeNode, std::string> >& variableNodes);
    void recordGlobalValue(double value, GlobalTarget target, ImageCustomIntegrator& integrator);
    void recordChangedParameters(ContextImpl& context);
    bool evaluateCondition(int step);
    ComputeContext& cc;
    double energy;
    float energyFloat;
    int numGlobalVariables, sumWorkGroupSize;
    bool hasInitializedKernels, deviceGlobalsAreCurrent, modifiesParameters, hasAnyConstraints, needsEnergyParamDerivs;
    std::vector<bool> deviceValuesAreCurrent;
    mutable std::vector<bool> localValuesAreCurrent;
    ComputeArray globalValues, sumBuffer, summedValue;
    ComputeArray uniformRandoms, randomSeed, perDofEnergyParamDerivs;
    std::vector<ComputeArray> tabulatedFunctions, perDofValues;
    std::map<int, double> savedEnergy;
    std::map<int, ComputeArray> savedForces;
    std::set<int> validSavedForces;
    mutable std::vector<std::vector<mm_float4> > localPerDofValuesFloat;
    mutable std::vector<std::vector<mm_double4> > localPerDofValuesDouble;
    std::map<std::string, double> energyParamDerivs;
    std::vector<std::string> perDofEnergyParamDerivNames;
    std::vector<double> localPerDofEnergyParamDerivs;
    std::vector<double> localGlobalValues;
    std::vector<double> initialGlobalVariables;
    std::vector<std::vector<ComputeKernel> > kernels;
    ComputeKernel randomKernel, kineticEnergyKernel, sumKineticEnergyKernel;
    std::vector<ImageCustomIntegrator::ComputationType> stepType;
    std::vector<ImageCustomIntegratorUtilities::Comparison> comparisons;
    std::vector<std::vector<Lepton::CompiledExpression> > globalExpressions;
    CompiledExpressionSet expressionSet;
    std::vector<bool> needsGlobals, needsForces, needsEnergy;
    std::vector<bool> computeBothForceAndEnergy, invalidatesForces, merged;
    std::vector<int> forceGroupFlags, blockEnd, requiredGaussian, requiredUniform;
    std::vector<int> stepEnergyVariableIndex, globalVariableIndex, parameterVariableIndex;
    int gaussianVariableIndex, uniformVariableIndex, dtVariableIndex;
    std::vector<std::string> parameterNames;
    std::vector<GlobalTarget> stepTarget;

    //double zmax;
};

class CommonIntegrateImageCustomStepKernel::GlobalTarget {
public:
    CommonIntegrateImageCustomStepKernel::GlobalTargetType type;
    int variableIndex;
    GlobalTarget() {
    }
    GlobalTarget(CommonIntegrateImageCustomStepKernel::GlobalTargetType type, int variableIndex) : type(type), variableIndex(variableIndex) {
    }
};

/**
 * This kernel is invoked by ImageLangevinIntegrator to update image particles.
 */
class CommonImageParticleKernel : public ImageParticleKernel {
public:
    CommonImageParticleKernel(const std::string & name, const Platform& platform, ComputeContext& cc) :
            ImageParticleKernel(name, platform), cc(cc), hasInitializedKernels(false){
    }

    void initialize(const System& system, const ImageIntegrator& integrator);

    void updateImagePositions(ContextImpl& context, const ImageIntegrator& integrator);

private:
    ComputeContext& cc;
    bool hasInitializedKernels;
    ComputeArray imagePairs, invAtomOrder;
    ComputeKernel kernelImage, kernelRecord;
    double zmax;
};

/**
 * This kernel is invoked by MCZBarostat to adjust the periodic box volume
 */
class CommonApplyMCZBarostatKernel : public ApplyMCZBarostatKernel {
public:
    CommonApplyMCZBarostatKernel(std::string name, const Platform& platform, ComputeContext& cc) : ApplyMCZBarostatKernel(name, platform), cc(cc),
            hasInitializedKernels(false) {
    }
    /**
     * Initialize the kernel.
     *
     * @param system     the System this kernel will be applied to
     * @param barostat   the MonteCarloBarostat this kernel will be used for
     * @param rigidMolecules  whether molecules should be kept rigid while scaling coordinates
     */
    void initialize(const System& system, const Force& barostat, bool rigidMolecules=true);
    /**
     * Attempt a Monte Carlo step, scaling particle positions (or cluster centers) by a specified value.
     * This version scales the x, y, and z positions independently.
     * This is called BEFORE the periodic box size is modified.  It should begin by translating each particle
     * or cluster into the first periodic box, so that coordinates will still be correct after the box size
     * is changed.
     *
     * @param context    the context in which to execute this kernel
     * @param scaleZ     the scale factor by which to multiply particle z-coordinate
     */
    void scaleCoordinates(ContextImpl& context, double scaleZ);
    /**
     * Reject the most recent Monte Carlo step, restoring the particle positions to where they were before
     * scaleCoordinates() was last called.
     *
     * @param context    the context in which to execute this kernel
     */
    void restoreCoordinates(ContextImpl& context);
private:
    ComputeContext& cc;
    bool hasInitializedKernels, rigidMolecules;
    int numMolecules;
    ComputeArray savedPositions, savedFloatForces, savedLongForces;
    ComputeArray moleculeAtoms;
    ComputeArray moleculeStartIndex;
    ComputeKernel kernel;
    std::vector<int> lastAtomOrder;
};

} // namespace ImagePlugin

#endif /*COMMON_IMAGE_KERNELS_H_*/

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

#include "openmm/ImageCustomIntegrator.h"
#include "openmm/Context.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/AssertionUtilities.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/ImageKernels.h"
#include "lepton/CompiledExpression.h"
#include "lepton/ParsedExpression.h"
#include "lepton/Parser.h"
#include <set>
#include <string>
#include <iostream>
#include "lepton/CustomFunction.h"
#include "lepton/ExpressionProgram.h"
#include "lepton/ExpressionTreeNode.h"
#include "lepton/Operation.h"

using namespace OpenMM;
using namespace std;

ImageCustomIntegrator::ImageCustomIntegrator(double stepSize, int numCells, double zmax) : globalsAreCurrent(true), forcesAreValid(false) {
    setStepSize(stepSize);
    setConstraintTolerance(1e-5);
    setRandomNumberSeed(0);
    //cout<<"set random seed done!"<<endl;
    setKineticEnergyExpression("m*v*v/2");
    //Lepton::ParsedExpression expr1 = Lepton::Parser::parse("m*v*v/2");
    //cout<<"parse done!"<<endl;
    //expr1.createCompiledExpression();
    //cout<<"compiled expr done!"<<endl;
    //cout<<"set kinetic expression done!"<<endl;
    setNumCells(numCells);
    setCellSize(zmax);
    setUseImageParticle(true);
    //cout<<"ImageCustomIntegrator init done!"<<endl;
}

ImageCustomIntegrator::~ImageCustomIntegrator() {
    for (auto function : functions)
        delete function.function;
}

void ImageCustomIntegrator::initialize(ContextImpl& contextRef) {
    if (owner != NULL && &contextRef.getOwner() != owner)
        throw OpenMMException("This Integrator is already bound to a context");
    vector<std::string> variableList;
    set<std::string> variableSet;
    variableList.insert(variableList.end(), globalNames.begin(), globalNames.end());
    variableList.insert(variableList.end(), perDofNames.begin(), perDofNames.end());
    for (auto& name : variableList) {
        if (variableSet.find(name) != variableSet.end())
            throw OpenMMException("The Integrator defines two variables with the same name: "+name);
        variableSet.insert(name);
        if (contextRef.getParameters().find(name) != contextRef.getParameters().end())
            throw OpenMMException("The Integrator defines a variable with the same name as a Context parameter: "+name);
    }
    set<std::string> globalTargets;
    globalTargets.insert(globalNames.begin(), globalNames.end());
    globalTargets.insert("dt");
    for (auto& param : contextRef.getParameters())
        globalTargets.insert(param.first);
    for (int i = 0; i < computations.size(); i++) {
        if (computations[i].type == ComputeGlobal && globalTargets.find(computations[i].variable) == globalTargets.end())
            throw OpenMMException("Unknown global variable: "+computations[i].variable);
    }
    context = &contextRef;
    owner = &contextRef.getOwner();
    kernel = context->getPlatform().createKernel(IntegrateImageCustomStepKernel::Name(), contextRef);
    kernel.getAs<IntegrateImageCustomStepKernel>().initialize(contextRef.getSystem(), *this);
    kernel.getAs<IntegrateImageCustomStepKernel>().setGlobalVariables(contextRef, globalValues);
    for (int i = 0; i < (int) perDofValues.size(); i++) {
        if (perDofValues[i].size() == 1)
            perDofValues[i].resize(context->getSystem().getNumParticles(), perDofValues[i][0]);
        kernel.getAs<IntegrateImageCustomStepKernel>().setPerDofVariable(contextRef, i, perDofValues[i]);
    }
    if(getUseImageParticle()){
        imgkernel = context->getPlatform().createKernel(ImageParticleKernel::Name(), contextRef);
        imgkernel.getAs<ImageParticleKernel>().initialize(contextRef.getSystem(), *this);
    }
}

void ImageCustomIntegrator::cleanup() {
    kernel = Kernel();
    imgkernel = Kernel();
}

void ImageCustomIntegrator::stateChanged(State::DataType changed) {
    forcesAreValid = false;
}

vector<string> ImageCustomIntegrator::getKernelNames() {
    vector<string> names;
    names.push_back(IntegrateImageCustomStepKernel::Name());
    if(getUseImageParticle()) {
        names.push_back(ImageParticleKernel::Name());
    }
    return names;
}

double ImageCustomIntegrator::computeKineticEnergy() {
    forcesAreValid = keNeedsForce;
    return kernel.getAs<IntegrateImageCustomStepKernel>().computeKineticEnergy(*context, *this, forcesAreValid);
}

bool ImageCustomIntegrator::kineticEnergyRequiresForce() const {
    return keNeedsForce;
}

void ImageCustomIntegrator::createCheckpoint(std::ostream& stream) const {
    for (int i = 0; i < getNumGlobalVariables(); i++) {
        double value = getGlobalVariable(i);
        stream.write((char*) &value, sizeof(double));
    }
    vector<Vec3> values;
    for (int i = 0; i < getNumPerDofVariables(); i++) {
        getPerDofVariable(i, values);
        stream.write((char*) values.data(), sizeof(Vec3)*values.size());
    }
}

void ImageCustomIntegrator::loadCheckpoint(std::istream& stream) {
    double value;
    for (int i = 0; i < getNumGlobalVariables(); i++) {
        stream.read((char*) &value, sizeof(double));
        setGlobalVariable(i, value);
    }
    // set the per-DOF for each particles
    vector<Vec3> values(context->getSystem().getNumParticles());
    for (int i = 0; i < getNumPerDofVariables(); i++) {
        stream.read((char*) values.data(), sizeof(Vec3)*values.size());
        setPerDofVariable(i, values);
    }
}

void ImageCustomIntegrator::step(int steps) {
    if (context == NULL)
        throw OpenMMException("This Integrator is not bound to a context!");
    globalsAreCurrent = false;
    for (int i = 0; i < steps; ++i) {
        kernel.getAs<IntegrateImageCustomStepKernel>().execute(*context, *this, forcesAreValid);
        if(getUseImageParticle()){
            imgkernel.getAs<ImageParticleKernel>().updateImagePositions(*context, *this);
            //cout<<"ImageParticleKernel excution done!"<<endl;
        }
    }
}

int ImageCustomIntegrator::addGlobalVariable(const string& name, double initialValue) {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    globalNames.push_back(name);
    globalValues.push_back(initialValue);
    return globalNames.size()-1;
}

const string& ImageCustomIntegrator::getGlobalVariableName(int index) const {
    ASSERT_VALID_INDEX(index, globalNames);
    return globalNames[index];
}

int ImageCustomIntegrator::addPerDofVariable(const string& name, double initialValue) {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    perDofNames.push_back(name);
    perDofValues.push_back(vector<Vec3>(1, Vec3(initialValue, initialValue, initialValue)));
    return perDofNames.size()-1;
}

const string& ImageCustomIntegrator::getPerDofVariableName(int index) const {
    ASSERT_VALID_INDEX(index, perDofNames);
    return perDofNames[index];
}

double ImageCustomIntegrator::getGlobalVariable(int index) const {
    ASSERT_VALID_INDEX(index, globalValues);
    if (owner != NULL && !globalsAreCurrent) {
        kernel.getAs<const IntegrateImageCustomStepKernel>().getGlobalVariables(*context, globalValues);
        globalsAreCurrent = true;
    }
    return globalValues[index];
}

double ImageCustomIntegrator::getGlobalVariableByName(const string& name) const {
    for (int i = 0; i < (int) globalNames.size(); i++) {
        if (name == globalNames[i]) {
            return getGlobalVariable(i);
        }
    }
    throw OpenMMException("Illegal global variable name: "+name);
}

void ImageCustomIntegrator::setGlobalVariable(int index, double value) {
    ASSERT_VALID_INDEX(index, globalValues);
    if (owner != NULL && !globalsAreCurrent) {
        kernel.getAs<IntegrateImageCustomStepKernel>().getGlobalVariables(*context, globalValues);
        globalsAreCurrent = true;
    }
    globalValues[index] = value;
    if (owner != NULL)
        kernel.getAs<IntegrateImageCustomStepKernel>().setGlobalVariables(*context, globalValues);
}

void ImageCustomIntegrator::setGlobalVariableByName(const string& name, double value) {
    for (int i = 0; i < (int) globalNames.size(); i++)
        if (name == globalNames[i]) {
            setGlobalVariable(i, value);
            return;
        }
    throw OpenMMException("Illegal global variable name: "+name);
}

void ImageCustomIntegrator::getPerDofVariable(int index, vector<Vec3>& values) const {
    ASSERT_VALID_INDEX(index, perDofValues);
    if (owner == NULL)
        values = perDofValues[index];
    else
        kernel.getAs<const IntegrateImageCustomStepKernel>().getPerDofVariable(*context, index, values);
}

void ImageCustomIntegrator::getPerDofVariableByName(const string& name,  vector<Vec3>& values) const {
    for (int i = 0; i < (int) perDofNames.size(); i++) {
        if (name == perDofNames[i]) {
            getPerDofVariable(i, values);
            return;
        }
    }
    throw OpenMMException("Illegal per-DOF variable name: "+name);
}

void ImageCustomIntegrator::setPerDofVariable(int index, const vector<Vec3>& values) {
    ASSERT_VALID_INDEX(index, perDofValues);
    if (owner != NULL && values.size() != context->getSystem().getNumParticles())
        throw OpenMMException("setPerDofVariable() called with wrong number of values");
    if (owner == NULL)
        perDofValues[index] = values;
    else
        kernel.getAs<IntegrateImageCustomStepKernel>().setPerDofVariable(*context, index, values);
}

void ImageCustomIntegrator::setPerDofVariableByName(const string& name, const vector<Vec3>& value) {
    for (int i = 0; i < (int) perDofNames.size(); i++)
        if (name == perDofNames[i]) {
            setPerDofVariable(i, value);
            return;
        }
    throw OpenMMException("Illegal per-DOF variable name: "+name);
}

int ImageCustomIntegrator::addComputeGlobal(const string& variable, const string& expression) {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    computations.push_back(ComputationInfo(ComputeGlobal, variable, expression));
    return computations.size()-1;
}

int ImageCustomIntegrator::addComputePerDof(const string& variable, const string& expression) {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    computations.push_back(ComputationInfo(ComputePerDof, variable, expression));
    return computations.size()-1;
}

int ImageCustomIntegrator::addComputeSum(const string& variable, const string& expression) {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    computations.push_back(ComputationInfo(ComputeSum, variable, expression));
    return computations.size()-1;
}

int ImageCustomIntegrator::addConstrainPositions() {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    computations.push_back(ComputationInfo(ConstrainPositions, "", ""));
    return computations.size()-1;
}

int ImageCustomIntegrator::addConstrainVelocities() {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    computations.push_back(ComputationInfo(ConstrainVelocities, "", ""));
    return computations.size()-1;
}

int ImageCustomIntegrator::addUpdateContextState() {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    computations.push_back(ComputationInfo(UpdateContextState, "", ""));
    return computations.size()-1;
}

int ImageCustomIntegrator::beginIfBlock(const string& expression) {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    computations.push_back(ComputationInfo(IfBlockStart, "", expression));
    return computations.size()-1;
}

int ImageCustomIntegrator::beginWhileBlock(const string& expression) {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    computations.push_back(ComputationInfo(WhileBlockStart, "", expression));
    return computations.size()-1;
}

int ImageCustomIntegrator::endBlock() {
    if (owner != NULL)
        throw OpenMMException("The integrator cannot be modified after it is bound to a context");
    computations.push_back(ComputationInfo(BlockEnd, "", ""));
    return computations.size()-1;
}

void ImageCustomIntegrator::getComputationStep(int index, ComputationType& type, string& variable, string& expression) const {
    ASSERT_VALID_INDEX(index, computations);
    type = computations[index].type;
    variable = computations[index].variable;
    expression = computations[index].expression;
}

int ImageCustomIntegrator::addTabulatedFunction(const std::string& name, TabulatedFunction* function) {
    functions.push_back(FunctionInfo(name, function));
    return functions.size()-1;
}

const TabulatedFunction& ImageCustomIntegrator::getTabulatedFunction(int index) const {
    ASSERT_VALID_INDEX(index, functions);
    return *functions[index].function;
}

TabulatedFunction& ImageCustomIntegrator::getTabulatedFunction(int index) {
    ASSERT_VALID_INDEX(index, functions);
    return *functions[index].function;
}

const string& ImageCustomIntegrator::getTabulatedFunctionName(int index) const {
    ASSERT_VALID_INDEX(index, functions);
    return functions[index].name;
}

const string& ImageCustomIntegrator::getKineticEnergyExpression() const {
    return kineticEnergy;
}


void ImageCustomIntegrator::setKineticEnergyExpression(const string& expression) {
    //cout<<"start set Kinetic expression..."<<endl;
    kineticEnergy = expression;
    //cout<<kineticEnergy<<endl;
    // Lepton::ParsedExpression expr1 = Lepton::Parser::parse(kineticEnergy);
    // Lepton::CompiledExpression expr = expr1.createCompiledExpression();
    Lepton::CompiledExpression expr = Lepton::Parser::parse(kineticEnergy).createCompiledExpression();
    // for(const auto& str : expr.getVariables()){
    //     cout<<str<<endl;
    // }
    //cout<<"expr done!"<<endl;
    keNeedsForce = (expr.getVariables().find("f") != expr.getVariables().end());
}

void ImageCustomIntegrator::serializeParameters(SerializationNode& node) const {
    node.setIntProperty("version", 1);
    SerializationNode& globalVariablesNode = node.createChildNode("GlobalVariables");
    for (int i = 0; i < getNumGlobalVariables(); i++)
        globalVariablesNode.setDoubleProperty(getGlobalVariableName(i), getGlobalVariable(i));
    SerializationNode& perDofVariablesNode = node.createChildNode("PerDofVariables");
    for (int i = 0; i < getNumPerDofVariables(); i++) {
        SerializationNode& perDofValuesNode = perDofVariablesNode.createChildNode(getPerDofVariableName(i));
        vector<Vec3> perDofValues;
        getPerDofVariable(i, perDofValues);
        for (int j = 0; j < perDofValues.size(); j++)
            perDofValuesNode.createChildNode("Value").setDoubleProperty("x",perDofValues[j][0]).setDoubleProperty("y",perDofValues[j][1]).setDoubleProperty("z",perDofValues[j][2]);
    }
}

void ImageCustomIntegrator::deserializeParameters(const SerializationNode& node) {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    const SerializationNode& globalVariablesNode = node.getChildNode("GlobalVariables");
    for (auto& prop : globalVariablesNode.getProperties())
        setGlobalVariableByName(prop.first, globalVariablesNode.getDoubleProperty(prop.first));
    const SerializationNode& perDofVariablesNode = node.getChildNode("PerDofVariables");
    for (auto& var : perDofVariablesNode.getChildren()) {
        vector<Vec3> perDofValues;
        for (auto& child : var.getChildren())
            perDofValues.push_back(Vec3(child.getDoubleProperty("x"), child.getDoubleProperty("y"), child.getDoubleProperty("z")));
        setPerDofVariableByName(var.getName(), perDofValues);
    }
}


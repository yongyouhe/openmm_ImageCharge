/* -------------------------------------------------------------------------- *
 *                                OpenMMImage                                 *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: RC Wang                                                     *
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

#include "openmm/serialization/ImageLangevinIntegratorProxy.h"
#include "openmm/ImageLangevinIntegrator.h"
#include "openmm/serialization/SerializationNode.h"
#include <sstream>

using namespace OpenMM;
using namespace std;
//
ImageLangevinIntegratorProxy::ImageLangevinIntegratorProxy() : SerializationProxy("ImageLangevinIntegrator") {
}

void ImageLangevinIntegratorProxy::serialize(const void* object, SerializationNode& node) const {
    node.setIntProperty("version", 1);
    const ImageLangevinIntegrator& integrator = *reinterpret_cast<const ImageLangevinIntegrator*>(object);
    node.setDoubleProperty("stepSize", integrator.getStepSize());
    node.setDoubleProperty("constraintTolerance", integrator.getConstraintTolerance());
    node.setDoubleProperty("temperature", integrator.getTemperature());
    node.setDoubleProperty("friction", integrator.getFriction());
    node.setIntProperty("randomSeed", integrator.getRandomNumberSeed());
    node.setIntProperty("numCells", integrator.getNumCells());
    node.setDoubleProperty("cellsize", integrator.getCellSize());
    node.setBoolProperty("useImageParticle", integrator.getUseImageParticle());
    SerializationNode& imagePairs = node.createChildNode("ImagePairs");
    for(auto& pair : integrator.getImagePairs()){
        int image = pair.first;
        int parent = pair.second;
        imagePairs.createChildNode("Pair").setIntProperty("image", image).setIntProperty("parent", parent);
    }
}

void* ImageLangevinIntegratorProxy::deserialize(const SerializationNode& node) const {
    if (node.getIntProperty("version") != 1)
        throw OpenMMException("Unsupported version number");
    ImageLangevinIntegrator *integrator = new ImageLangevinIntegrator(node.getDoubleProperty("temperature"),
            node.getDoubleProperty("frction"), node.getDoubleProperty("stepSize"), 
            node.getIntProperty("numCells"), node.getDoubleProperty("cellsize"));
    integrator->setConstraintTolerance(node.getDoubleProperty("constraintTolerance"));
    integrator->setRandomNumberSeed(node.getIntProperty("randomSeed"));
    integrator->setUseImageParticle(node.getBoolProperty("useImageParticle"));
    const SerializationNode& imagePairs = node.getChildNode("ImagePairs");
    for(auto& pair : imagePairs.getChildren())
        integrator->setImagePair(pair.getIntProperty("image"), pair.getIntProperty("parent"));
    
    return integrator;
}

#ifndef OPENMM_SLABCORRECTION_H_
#define OPENMM_SLABCORRECTION_H_

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

#include "openmm/Context.h"
#include "openmm/Force.h"
#include <vector>
#include "openmm/internal/windowsExport.h"

namespace OpenMM {

/**
 * This class implements an slab correction of two parallel plates system. (Yeh JCP, 1999)
 */

class OPENMM_EXPORT SlabCorrection : public Force {
public:

    SlabCorrection(bool applytoAll=true, bool useAmoebaDip=true);
    /**
     * Get whether apply the slab correction to all atoms.
    */
    bool getApplytoAll() const {
        return applytoAll;
    }
    /**
     * Set whether apply the slab correction to all atoms.
    */
    void setApplytoAll(bool apply) {
        applytoAll = apply;
    }
    /**
     * Add a particle to the correction.
     * @param index     the index of added particle
     * @return          the number of particles applied the slab correction
     */
    int addParticles(int index) {
        particlesCorr.push_back(index);
        return particlesCorr.size();
    }
    /**
     * Get particles which are applied to correction
     */
    const std::vector<int>& getParticlesCorr() const {
        return particlesCorr;
    }
    int getNumParticlesCorr() const {
        return particlesCorr.size();
    }
    /**
     * Get whether to use the dipoles in Amoeba force field (including permenant and induced dipoles).
    */
    bool useAmoebaDipole() const {
        return useAmoebaDip;
    }
    /**
     * Set whether use the dipoles in Amoeba force field.
    */
    void setUseAmoebaDipole(bool use) {
        useAmoebaDip = use;
    }
    /**
     * Update the per-bond parameters in a Context to match those stored in this Force object.  This method provides
     * an efficient method to update certain parameters in an existing Context without needing to reinitialize it.
     * Simply call setBondParameters() to modify this object's parameters, then call updateParametersInState()
     * to copy them over to the Context.
     * 
     * The only information this method updates is the values of per-bond parameters.  The set of particles involved
     * in a bond cannot be changed, nor can new bonds be added.
     */
    void updateParametersInContext(Context& context);
    /**
     * Returns whether or not this force makes use of periodic boundary
     * conditions.
     *
     * @returns true if force uses PBC and false otherwise
     */
    bool usesPeriodicBoundaryConditions() const {
        return false;
    }

protected:
    ForceImpl* createImpl() const;
private:
    bool applytoAll, useAmoebaDip;
    std::vector<int> particlesCorr;
};

}
#endif

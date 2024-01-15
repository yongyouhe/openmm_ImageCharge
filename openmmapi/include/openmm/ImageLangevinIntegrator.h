#ifndef OPENMM_IMAGELANGEVININTEGRATOR_H_
#define OPENMM_IMAGELANGEVININTEGRATOR_H_

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

#include "openmm/Integrator.h"
#include "openmm/ImageIntegrator.h"
#include "openmm/Kernel.h"
#include "openmm/internal/windowsExport.h"
#include <algorithm>
#include <map>

namespace OpenMM{

class OPENMM_EXPORT ImageLangevinIntegrator : public Integrator, public ImageIntegrator {
    
    public:
        /**
        * Create a ImageLangevinIntegrator.
        * 
        * @param temperature    the temperature of the heat bath (in Kelvin)
        * @param frictionCoeff  the friction coefficient which couples the system to the heat bath (in inverse picoseconds)
        * @param stepSize       the step size with which to integrate the system (in picoseconds)
        * @param numCells       the number of cells including real and image cells (default = 2, indicating single real+image pair)
        * @param zmax           the size of unit cell in z-dimension
        */
        ImageLangevinIntegrator(double temperature, double frictionCoeff, double stepSize, int numCells=2, double zmax=-1);
        /**
         * Get the temperature of the heat bath (in Kelvin).
         * @return the temperature of the heat bath, measured in Kelvin
         */
        double getTemperature() const {
            return temperature;
        }
        /**
         * Set the temperature of the heat bath (in Kelvin).
         * @param temp  the temperature of the heat bath, measured in Kelvin.
         */
        void setTemperature(double temp);
        /**
         * Get the friction coefficient which determines how strongly the system is coupled to the heat bath (in inverse ps)
         * @return the friction coefficient, measured in 1/ps
         */
        double getFriction() const {
            return friction;
        }
        /**
         * Set the friction coefficient which determined how strongly the system is coupled to the heat bath (in inverse ps)
         * @param coeff  the friction coefficient, measured in 1/ps
         */
        void setFriction(double coeff);
        /**
         * Get the random number seed. See setRandomNumberSeed() for details.
         */
        int getRandomNumberSeed() const {
            return randomNumberSeed;
        }
        /**
         * Set the random number seed.  The precise meaning of this parameter is undefined, and is left up
         * to each Platform to interpret in an appropriate way.  It is guaranteed that if two simulations
         * are run with different random number seeds, the sequence of random forces will be different.  On
         * the other hand, no guarantees are made about the behavior of simulations that use the same seed.
         * In particular, Platforms are permitted to use non-deterministic algorithms which produce different
         * results on successive runs, even if those runs were initialized identically.
         *
         * If seed is set to 0 (which is the default value assigned), a unique seed is chosen when a Context
         * is created from this Force. This is done to ensure that each Context receives unique random seeds
         * without you needing to set them explicitly.
         */
        void setRandomNumberSeed(int seed) {
            randomNumberSeed = seed;
        }
        /**
         * Advance a simulation through time by taking a series of time steps.
         * 
         * @param steps   the number of time steps to take
         */
        void step(int steps);
        // /**
        //  * Get the number of real+image cells.
        //  * @return the number of cells.
        //  */
        // int getNumCells() const {
        //     return numCells;
        // }
        // /**
        //  * Set the number of real+image cells.
        //  * @param numCells the number of real+image cells (must be multiple of 2, value larger than 2 indicates multiple real+image pair)
        //  */
        // void setNumCells(int cells) {
        //     numCells = cells;
        // }
        // /**
        //  * Get the size of cell in z-dimension
        //  */
        // double getCellSize() const {
        //     return cellsize;
        // }
        // /**
        //  * Set the cell size in z-dimension
        //  */
        // void setCellSize(double zmax) {
        //     cellsize = zmax;
        // }
        // /**
        //  * Choose whether to use image particles.
        //  */
        // void setUseImageParticle(bool use) {
        //     useImageParticle = use;
        // }
        // /**
        //  * Get whether to use image particle.
        //  */
        // const bool& getUseImageParticle() const {
        //     return useImageParticle;
        // }
        // /**
        //  * Set a particle as the image of another particle.
        //  * @param image     the index of image particle
        //  * @param parent    the index of parent of the image particle
        //  * @return          the number of image particles in the system
        //  */
        // int setImagePair(int image, int parent);
        // /**
        //  * Get all image pairs
        //  */
        // const std::vector<std::pair<int,int> >& getImagePairs() const{
        //     return imagePairs;
        // }
        // /**
        //  * Check whether the particle i is a image particle
        //  * @return  bool
        //  */
        // bool isParticleImage(int i) const{
        //     return std::find(imageParticles.begin(), imageParticles.end(), i) != imageParticles.end(); 
        // }
    protected:
        /**
         * This will be called by the Context when it is created.  It informs the Integrator
         * of what context it will be integrating, and gives it a chance to do any necessary initialization.
         * It will also get called again if the application calls reinitialize() on the Context.
         */
        void initialize(ContextImpl& context);
        /**
         * This will be called by the Context when it is destroyed to let the Integrator do any necessary
         * cleanup.  It will also get called again if the application calls reinitialize() on the Context.
         */
        void cleanup();
        /**
         * Get the names of all Kernels used by this Integrator.
         */
        std::vector<std::string> getKernelNames();
        /**
         * Compute the kinetic energy of the system at the current time.
         */
        double computeKineticEnergy();
        /**
         * Get the time interval by which velocities are offset from positions. This is used to adjust velocities 
         * when setVelocitiesToTemperature() is called on a Context.
         */
        double getVelocityTimeOffset() const {
            return getStepSize()/2;
        }

    private:
        int randomNumberSeed;
        double temperature, friction;
        Kernel ilkernel, imgkernel;
        //int numCells;
        //double cellsize;
        //bool useImageParticle;

        //std::vector<std::pair<int,int> > imagePairs;
        //std::vector<int> imageParticles;

};

}

#endif
#ifndef OPENMM_IMAGEINTEGRATOR_H_
#define OPENMM_IMAGEINTEGRATOR_H_

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
#include "openmm/Kernel.h"
#include "openmm/internal/windowsExport.h"
#include <algorithm>
#include <map>

namespace OpenMM{

class OPENMM_EXPORT ImageIntegrator {
    
    public:
        ImageIntegrator(int numCells=2, double zmax=-1){
            setNumCells(numCells);
            setCellSize(zmax);
            useImageParticle = true;
        }

        /**
         * Get the number of real+image cells.
         * @return the number of cells.
         */
        int getNumCells() const {
            return numCells;
        }
        /**
         * Set the number of real+image cells.
         * @param numCells the number of real+image cells (must be multiple of 2, value larger than 2 indicates multiple real+image pair)
         */
        void setNumCells(int cells) {
            numCells = cells;
        }
        /**
         * Get the size of cell in z-dimension
         */
        double getCellSize() const {
            return cellsize;
        }
        /**
         * Set the cell size in z-dimension
         */
        void setCellSize(double zmax) {
            cellsize = zmax;
        }
        /**
         * Choose whether to use image particles.
         */
        void setUseImageParticle(bool use) {
            useImageParticle = use;
        }
        /**
         * Get whether to use image particle.
         */
        const bool& getUseImageParticle() const {
            return useImageParticle;
        }
        /**
         * Set a particle as the image of another particle.
         * @param image     the index of image particle
         * @param parent    the index of parent of the image particle
         * @return          the number of image particles in the system
         */
        int setImagePair(int image, int parent) {
            imageParticles.push_back(image);
            imagePairs.emplace_back(std::make_pair(image, parent));
            return imagePairs.size();
        }
        /**
         * Get all image pairs
         */
        const std::vector<std::pair<int,int> >& getImagePairs() const{
            return imagePairs;
        }
        /**
         * Check whether the particle i is a image particle
         * @return  bool
         */
        bool isParticleImage(int i) const{
            return std::find(imageParticles.begin(), imageParticles.end(), i) != imageParticles.end(); 
        }

    private:
        int numCells;
        double cellsize;
        bool useImageParticle;

        std::vector<std::pair<int,int> > imagePairs;
        std::vector<int> imageParticles;

};

}

#endif
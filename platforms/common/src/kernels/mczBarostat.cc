/**
 * Scale the particle positions with each axis independent.
 */

KERNEL void scalePositions(float scaleZ, int numMolecules, real4 periodicBoxSize,
        real4 invPeriodicBoxSize, real4 periodicBoxVecX, real4 periodicBoxVecY, real4 periodicBoxVecZ, GLOBAL real4* RESTRICT posq, 
        GLOBAL mixed4* RESTRICT velm, GLOBAL const int* RESTRICT moleculeAtoms, GLOBAL const int* RESTRICT moleculeStartIndex) {
    for (int index = GLOBAL_ID; index < numMolecules; index += GLOBAL_SIZE) {
        int first = moleculeStartIndex[index];
        int last = moleculeStartIndex[index+1];
        int numAtoms = last-first;

        if(numAtoms==1){
            real4 pos = posq[moleculeAtoms[first]];
            if(velm[moleculeAtoms[first]].w==0.0 && pos.z>periodicBoxVecZ.z*0.1){
                pos.z *= scaleZ;
                posq[moleculeAtoms[first]] = pos;
            }
        }
        else{
            // Find the center of each molecule.

            real centerZ = 0;
            for (int atom = first; atom < last; atom++) {
                real4 pos = posq[moleculeAtoms[atom]];
                centerZ += pos.z;
            }
            real invNumAtoms = RECIP((real) numAtoms);
            centerZ *= invNumAtoms;

            // Now scale the position of the molecule center.

            real delta;
            delta = centerZ*(scaleZ-1);
            for (int atom = first; atom < last; atom++) {
                real4 pos = posq[moleculeAtoms[atom]];
                pos.z += delta;
                posq[moleculeAtoms[atom]] = pos;
            }
        }

    }
}

/**
 * Add slab correction to the selected particles.
 */

KERNEL void addSlabCorrection(GLOBAL mm_long* forces, int bufferSize, bool applytoAll, GLOBAL real4* RESTRICT posq,
        int* RESTRICT invAtomOrder, GLOBAL const int* particlesCorr, bool useAmoebaDip, real muz, real fscale, 
        mm_long* RESTRICT sumQZ) {
    for (int index = GLOBAL_ID; index < NUM_PARTICLECORR; index += GLOBAL_SIZE) {
        if(useAmoebaDip == true) {
            if(applytoAll == true) {
                //forces[index] += mm_long (posq.w*muz*fscale*0x100000000);
                forces[index+PADDED_NUM_ATOMS*2] += (mm_long) (posq[index].w*muz*fscale*0x100000000);
                //ATOMIC_ADD(&forces[index+PADDED_NUM_ATOMS*2], (mm_ulong) ((mm_long) (posq[index].w*muz*fscale*0x100000000)));
            }
            else {
                int particle = particlesCorr[index];
                int iparticle = invAtomOrder[particle]+PADDED_NUM_ATOMS*2;
                forces[iparticle] += (mm_long) (posq[invAtomOrder[particle]].w*muz*fscale*0x100000000);
            }
        }
        else {
            if(applytoAll == true) {
                int iparticle = index+PADDED_NUM_ATOMS*2;
                forces[iparticle] += (mm_long) ((*sumQZ)*posq[index].w*fscale);
            }
            else {
                int particle = particlesCorr[index];
                int iparticle = invAtomOrder[particle]+PADDED_NUM_ATOMS*2;
                forces[iparticle] += (mm_long) ((*sumQZ)*posq[invAtomOrder[particle]].w*fscale);
            }
        }
    }
}

/**
 * Compute the sum of q*z
*/
KERNEL void computeQZ(int* RESTRICT invAtomOrder, bool applytoAll, GLOBAL const int* particlesCorr, mm_ulong* RESTRICT sumQZ, 
        GLOBAL real4* RESTRICT posq) {
    *sumQZ = 0;
    for (int index = GLOBAL_ID; index < NUM_PARTICLECORR; index += GLOBAL_SIZE) {
        if(applytoAll == true) {
            ATOMIC_ADD(sumQZ, (mm_ulong) ((mm_long) (posq[index].z*posq[index].w*0x100000000)));
        }
        else {
            int particle = particlesCorr[index];
            int iparticle = invAtomOrder[particle];
            ATOMIC_ADD(sumQZ, (mm_ulong) ((mm_long) (posq[iparticle].z*posq[iparticle].w*0x100000000)));
        }

    }
}

/**
 * Record the atom index
 */
KERNEL void recordAtomIndexes(int numAtoms, int* RESTRICT atomIndex, int* RESTRICT invAtomOrder) {

    int index = GLOBAL_ID;
    while (index < numAtoms) {
        invAtomOrder[atomIndex[index]] = index; //invAtomOrder里atomindex对应index
        index += GLOBAL_SIZE;
    }
}

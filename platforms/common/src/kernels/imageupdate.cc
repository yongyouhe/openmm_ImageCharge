
/**
 * Update the positions of image particles.
 */
KERNEL void updateImageParticlePositions(GLOBAL real4* RESTRICT posq, GLOBAL real4* RESTRICT posqCorrection, 
        GLOBAL const int2* RESTRICT imagePairs, mixed zmax, int* RESTRICT invAtomOrder) {
    
    for(int i = GLOBAL_ID; i < NUM_IMAGES; i += GLOBAL_SIZE) {
        int2 pair = imagePairs[i];
        int indexImg = invAtomOrder[pair.x];
        int indexPar = invAtomOrder[pair.y];
        
        // if(posq[indexPar].z > zmax*2){
        //     printf("The previous parent particle coord is %f %f\n", posq[indexPar].z, posqCorrection[indexPar]);
        // }
        if(posq[indexPar].w != -posq[indexPar].w) {//charge of particle is not 0.
#ifdef USE_MIXED_PRECISION
            real4 pos1 = posq[indexPar];
            real4 pos2 = posqCorrection[indexPar];
            mixed4 pos = make_mixed4(pos1.x+(mixed)pos2.x, pos1.y+(mixed)pos2.y, pos1.z+(mixed)pos2.z, pos1.w);
#else
            real4 pos = posq[indexPar];
#endif
            pos.z = -pos.z +2*zmax;
            pos.w = posq[indexImg].w;
#ifdef USE_MIXED_PRECISION
            posq[indexImg] = make_real4((real) pos.x, (real) pos.y, (real) pos.z, (real) pos.w);
            posqCorrection[indexImg] = make_real4(pos.x-(real) pos.x, pos.y-(real) pos.y, pos.z-(real) pos.z, 0);
#else
            posq[indexImg] = pos;
#endif
        }

        // if(posq[indexPar].z > zmax*2){
        //     printf("i = %d, indexImg = %d, indexPar = %d, posq[img] = %f %f %f %f, posq[par] = %f %f %f %f\n", 
        //             i, indexImg, indexPar, posq[indexImg].x, posq[indexImg].y, posq[indexImg].z, posq[indexImg].w,
        //             posq[indexPar].x, posq[indexPar].y, posq[indexPar].z, posq[indexPar].w);
        // }

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
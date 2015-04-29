// CUDA_ART_getProj.cu
// many parts of the code might be improved (speed?? memory usage?) by using "index-pointers" ; always the same computations for different indiex-combinations ; don't know if the cuda compiler compresses this to safe memory for commands!

// in proj already raySum/rayWeight is kept!

__global__ void proj(const float * R, const int * volSize, const int * projSize, const float * vol, const float * recOrig, const float * origProj, float * proj) {
    // declare some useful varibles
    int count = 0;
    float det_pos_lab[3];
    float det_pos_sample[3];
    float ray_dir[3];
    short ray_sig[3];	
    float boundary_hit[3][2];
    float boundary_hit_pos[3][2][3];
    float current_k;
    int start_index[3];
    float start_pos[3];
    float t[3];
    float raySum = 0.0;
    //int internalProjSize[2];
    //internalProjSize[0] = (int) ((projSize[0] + 1 - gridTrace[0])/2.0);
    //internalProjSize[1] = (int) ((projSize[1] + 1 - gridTrace[1])/2.0);
    int noThreads = gridDim.x*gridDim.y*blockDim.x*blockDim.y*blockDim.z;
    int threadInd = (blockIdx.y*gridDim.x + blockIdx.x)*(blockDim.x*blockDim.y*blockDim.z) + blockDim.x*blockDim.y*threadIdx.z + blockDim.x*threadIdx.y + threadIdx.x; // just go through all indexed threads going to be started
    
    //while((count*noThreads+threadInd)<(internalProjSize[0]*internalProjSize[1])) {
    while((count*noThreads+threadInd)<(projSize[0]*projSize[1])) {
        // map the threads on the specified voxel indices (by gridTrace)
        int ind = threadInd+count*noThreads;
        //int indY = gridTrace[1] + 2*((int)(ind / internalProjSize[0]));
        //int indX = gridTrace[0] + 2*(ind - ((int)(ind / internalProjSize[0]))*internalProjSize[0]);
        int indY = ind / projSize[0];
        int indX = ind - indY*projSize[0];

        // adding 0.00001 to positions and directions in the case of even volSize but odd projection size leads to artifacts at the volume borders!
        det_pos_lab[0] = 0.0+0.00001;  // +0.00001 to avoid divisions by zero and other singularities of the raytracing algorithm --> better solution possible?!? (but shouldnot be a problem..)
        det_pos_lab[1] = indX+0.5-origProj[0]+0.00001;
        det_pos_lab[2] = indY+0.5-origProj[1]+0.00001;
        det_pos_sample[0] = R[0]*det_pos_lab[0] + R[3]*det_pos_lab[1] + R[6]*det_pos_lab[2] + recOrig[0] - volSize[0]/2.0;
        det_pos_sample[1] = R[1]*det_pos_lab[0] + R[4]*det_pos_lab[1] + R[7]*det_pos_lab[2] + recOrig[1] - volSize[1]/2.0;
        det_pos_sample[2] = R[2]*det_pos_lab[0] + R[5]*det_pos_lab[1] + R[8]*det_pos_lab[2] + recOrig[2] - volSize[2]/2.0;
        ray_dir[0] = R[0]; // might be improved by directly using R instead of an additional variable ray_dir // mustn't it be "ray_dir[i] = -R[i]" ?!?!?
        ray_dir[1] = R[1];
        ray_dir[2] = R[2];
        if (ray_dir[0]<0) ray_sig[0] = -1;
        if (ray_dir[0]>0) ray_sig[0] = 1;
        if (ray_dir[0]==0) {ray_sig[0] = 1; ray_dir[0] = 0.00001;} // is there a more elegant and safer (more general) way to avoid these indefinit expressions (divisions by zero)?!?
        if (ray_dir[1]<0) ray_sig[1] = -1;
        if (ray_dir[1]>0) ray_sig[1] = 1;
        if (ray_dir[1]==0) {ray_sig[1] = 1; ray_dir[1] = 0.00001;}
        if (ray_dir[2]<0) ray_sig[2] = -1;
        if (ray_dir[2]>0) ray_sig[2] = 1;
        if (ray_dir[2]==0) {ray_sig[2] = 1; ray_dir[2] = 0.00001;}

        // calculate intersections with the volume boundary
        // boundary_hit[0][0] = (volSize[0]/2.0-det_pos_sample[0]) / ray_dir[0];
        // boundary_hit[0][1] = (-volSize[0]/2.0-det_pos_sample[0]) / ray_dir[0];
        // boundary_hit[1][0] = (volSize[1]/2.0-det_pos_sample[1]) / ray_dir[1];
        // boundary_hit[1][1] = (-volSize[1]/2.0-det_pos_sample[1]) / ray_dir[1];
        // boundary_hit[2][0] = (volSize[2]/2.0-det_pos_sample[2]) / ray_dir[2];
        // boundary_hit[2][1] = (-volSize[2]/2.0-det_pos_sample[2]) / ray_dir[2];
        boundary_hit[0][0] = (volSize[0]/2.0-det_pos_sample[0]) / ray_dir[0];
        boundary_hit[0][1] = (-volSize[0]/2.0-det_pos_sample[0]) / ray_dir[0];
        boundary_hit[1][0] = (volSize[1]/2.0-det_pos_sample[1]) / ray_dir[1];
        boundary_hit[1][1] = (-volSize[1]/2.0-det_pos_sample[1]) / ray_dir[1];
        boundary_hit[2][0] = (volSize[2]/2.0-det_pos_sample[2]) / ray_dir[2];
        boundary_hit[2][1] = (-volSize[2]/2.0-det_pos_sample[2]) / ray_dir[2];

        // for the following optimization by using a loop possible? Or maybe put it all in a separate kernel and safe boundary matrix!
        boundary_hit_pos[0][0][0] = det_pos_sample[0] + ray_dir[0] * boundary_hit[0][0]; // might be optimized by only calculating the missing components and only until 2 possibilities are found!
        boundary_hit_pos[0][0][1] = det_pos_sample[1] + ray_dir[1] * boundary_hit[0][0];
        boundary_hit_pos[0][0][2] = det_pos_sample[2] + ray_dir[2] * boundary_hit[0][0];
        boundary_hit_pos[0][1][0] = det_pos_sample[0] + ray_dir[0] * boundary_hit[0][1]; // might be optimized by only calculating the missing components and only until 2 possibilities are found!
        boundary_hit_pos[0][1][1] = det_pos_sample[1] + ray_dir[1] * boundary_hit[0][1];
        boundary_hit_pos[0][1][2] = det_pos_sample[2] + ray_dir[2] * boundary_hit[0][1];
        boundary_hit_pos[1][0][0] = det_pos_sample[0] + ray_dir[0] * boundary_hit[1][0]; // might be optimized by only calculating the missing components and only until 2 possibilities are found!
        boundary_hit_pos[1][0][1] = det_pos_sample[1] + ray_dir[1] * boundary_hit[1][0];
        boundary_hit_pos[1][0][2] = det_pos_sample[2] + ray_dir[2] * boundary_hit[1][0];
        boundary_hit_pos[1][1][0] = det_pos_sample[0] + ray_dir[0] * boundary_hit[1][1]; // might be optimized by only calculating the missing components and only until 2 possibilities are found!
        boundary_hit_pos[1][1][1] = det_pos_sample[1] + ray_dir[1] * boundary_hit[1][1];
        boundary_hit_pos[1][1][2] = det_pos_sample[2] + ray_dir[2] * boundary_hit[1][1];
        boundary_hit_pos[2][0][0] = det_pos_sample[0] + ray_dir[0] * boundary_hit[2][0]; // might be optimized by only calculating the missing components and only until 2 possibilities are found!
        boundary_hit_pos[2][0][1] = det_pos_sample[1] + ray_dir[1] * boundary_hit[2][0];
        boundary_hit_pos[2][0][2] = det_pos_sample[2] + ray_dir[2] * boundary_hit[2][0];
        boundary_hit_pos[2][1][0] = det_pos_sample[0] + ray_dir[0] * boundary_hit[2][1]; // might be optimized by only calculating the missing components and only until 2 possibilities are found!
        boundary_hit_pos[2][1][1] = det_pos_sample[1] + ray_dir[1] * boundary_hit[2][1];
        boundary_hit_pos[2][1][2] = det_pos_sample[2] + ray_dir[2] * boundary_hit[2][1];

        current_k = 10000000000.0;
        short hit_count = 0;
        count++;

        // this part I do not like, yet: seems too complicated..
        if ( boundary_hit_pos[0][0][1]>-volSize[1]/2.0 && boundary_hit_pos[0][0][1]<volSize[1]/2.0 && boundary_hit_pos[0][0][2]>-volSize[2]/2.0 && boundary_hit_pos[0][0][2]<volSize[2]/2.0) {
            hit_count++;
            if (boundary_hit[0][0]<current_k) {
                current_k = boundary_hit[0][0];
                start_index[0] = volSize[0]-1;
                start_index[1] = (int) (boundary_hit_pos[0][0][1]+volSize[1]/2.0);
                start_index[2] = (int) (boundary_hit_pos[0][0][2]+volSize[2]/2.0);
                start_pos[0] = boundary_hit_pos[0][0][0];
                start_pos[1] = boundary_hit_pos[0][0][1];
                start_pos[2] = boundary_hit_pos[0][0][2];
                t[0] = ((float) ray_sig[0]) / ray_dir[0];
                // do it the easy way:
                if (ray_sig[1]>0) t[1] = (1.0 / ray_dir[1]) * (1.0 - ( start_pos[1] + volSize[1]/2.0 - (int) (start_pos[1] + volSize[1]/2.0) ) );
                else t[1] = (-1.0 / ray_dir[1]) * ( start_pos[1] + volSize[1]/2.0 - (int) (start_pos[1] + volSize[1]/2.0) );
                if (ray_sig[2]>0) t[2] = (1.0 / ray_dir[2]) * (1.0 - ( start_pos[2] + volSize[2]/2.0 - (int) (start_pos[2] + volSize[2]/2.0) ) );
                else t[2] = (-1.0 / ray_dir[2]) * ( start_pos[2] + volSize[2]/2.0 - (int) (start_pos[2] + volSize[2]/2.0) );
                }
            }
        if ( boundary_hit_pos[0][1][1]>-volSize[1]/2.0 && boundary_hit_pos[0][1][1]<volSize[1]/2.0 && boundary_hit_pos[0][1][2]>-volSize[2]/2.0 && boundary_hit_pos[0][1][2]<volSize[2]/2.0 ) {
            hit_count++;
            if (boundary_hit[0][1]<current_k) {
                current_k = boundary_hit[0][1];
                start_index[0] = 0;
                start_index[1] = (int) (boundary_hit_pos[0][1][1]+volSize[1]/2.0);
                start_index[2] = (int) (boundary_hit_pos[0][1][2]+volSize[2]/2.0);
                start_pos[0] = boundary_hit_pos[0][1][0];
                start_pos[1] = boundary_hit_pos[0][1][1];
                start_pos[2] = boundary_hit_pos[0][1][2];
                t[0] = ((float) ray_sig[0]) / ray_dir[0];
                if (ray_sig[1]>0) t[1] = (1.0 / ray_dir[1]) * (1.0 - ( start_pos[1] + volSize[1]/2.0 - (int) (start_pos[1] + volSize[1]/2.0) ) );
                else t[1] = (-1.0 / ray_dir[1]) * ( start_pos[1] + volSize[1]/2.0 - (int) (start_pos[1] + volSize[1]/2.0) );
                if (ray_sig[2]>0) t[2] = (1.0 / ray_dir[2]) * (1.0 - ( start_pos[2] + volSize[2]/2.0 - (int) (start_pos[2] + volSize[2]/2.0) ) );
                else t[2] = (-1.0 / ray_dir[2]) * ( start_pos[2] + volSize[2]/2.0 - (int) (start_pos[2] + volSize[2]/2.0) );
                }
            }
        if ( boundary_hit_pos[1][0][0]>-volSize[0]/2.0 && boundary_hit_pos[1][0][0]<volSize[0]/2.0 && boundary_hit_pos[1][0][2]>-volSize[2]/2.0 && boundary_hit_pos[1][0][2]<volSize[2]/2.0 ) {
            hit_count++;
            if (boundary_hit[1][0]<current_k) {
                current_k = boundary_hit[1][0];
                start_index[1] = volSize[1]-1;
                start_index[0] = (int) (boundary_hit_pos[1][0][0]+volSize[0]/2.0);
                start_index[2] = (int) (boundary_hit_pos[1][0][2]+volSize[2]/2.0);
                start_pos[0] = boundary_hit_pos[1][0][0];
                start_pos[1] = boundary_hit_pos[1][0][1];
                start_pos[2] = boundary_hit_pos[1][0][2];
                t[1] = ((float) ray_sig[1]) / ray_dir[1];
                if (ray_sig[0]>0) t[0] = (1.0 / ray_dir[0]) * (1.0 - ( start_pos[0] + volSize[0]/2.0 - (int) (start_pos[0] + volSize[0]/2.0) ) );
                else t[0] = (-1.0 / ray_dir[0]) * ( start_pos[0] + volSize[0]/2.0 - (int) (start_pos[0] + volSize[0]/2.0) );
                if (ray_sig[2]>0) t[2] = (1.0 / ray_dir[2]) * (1.0 - ( start_pos[2] + volSize[2]/2.0 - (int) (start_pos[2] + volSize[2]/2.0) ) );
                else t[2] = (-1.0 / ray_dir[2]) * ( start_pos[2] + volSize[2]/2.0 - (int) (start_pos[2] + volSize[2]/2.0) );
                }
            }
        if ( boundary_hit_pos[1][1][0]>-volSize[0]/2.0 && boundary_hit_pos[1][1][0]<volSize[0]/2.0 && boundary_hit_pos[1][1][2]>-volSize[2]/2.0 && boundary_hit_pos[1][1][2]<volSize[2]/2.0 ) {
            hit_count++;
            if (boundary_hit[1][1]<current_k) {
                current_k = boundary_hit[1][1];
                start_index[1] = 0;
                start_index[0] = (int) (boundary_hit_pos[1][1][0]+volSize[0]/2.0);
                start_index[2] = (int) (boundary_hit_pos[1][1][2]+volSize[2]/2.0);
                start_pos[0] = boundary_hit_pos[1][1][0];
                start_pos[1] = boundary_hit_pos[1][1][1];
                start_pos[2] = boundary_hit_pos[1][1][2];
                t[1] = ((float) ray_sig[1]) / ray_dir[1];
                if (ray_sig[0]>0) t[0] = (1.0 / ray_dir[0]) * (1.0 - ( start_pos[0] + volSize[0]/2.0 - (int) (start_pos[0] + volSize[0]/2.0) ) );
                else t[0] = (-1.0 / ray_dir[0]) * ( start_pos[0] + volSize[0]/2.0 - (int) (start_pos[0] + volSize[0]/2.0) );
                if (ray_sig[2]>0) t[2] = (1.0 / ray_dir[2]) * (1.0 - ( start_pos[2] + volSize[2]/2.0 - (int) (start_pos[2] + volSize[2]/2.0) ) );
                else t[2] = (-1.0 / ray_dir[2]) * ( start_pos[2] + volSize[2]/2.0 - (int) (start_pos[2] + volSize[2]/2.0) );
                }
            }
        if ( boundary_hit_pos[2][0][0]>-volSize[0]/2.0 && boundary_hit_pos[2][0][0]<volSize[0]/2.0 && boundary_hit_pos[2][0][1]>-volSize[1]/2.0 && boundary_hit_pos[2][0][1]<volSize[1]/2.0 ) {
            hit_count++;
            if (boundary_hit[2][0]<current_k) {
                current_k = boundary_hit[2][0];
                start_index[2] = volSize[2]-1;
                start_index[1] = (int) (boundary_hit_pos[2][0][1]+volSize[1]/2.0);
                start_index[0] = (int) (boundary_hit_pos[2][0][0]+volSize[0]/2.0);
                start_pos[0] = boundary_hit_pos[2][0][0];
                start_pos[1] = boundary_hit_pos[2][0][1];
                start_pos[2] = boundary_hit_pos[2][0][2];
                t[2] = ((float) ray_sig[2]) / ray_dir[2];
                if (ray_sig[0]>0) t[0] = (1.0 / ray_dir[0]) * (1.0 - ( start_pos[0] + volSize[0]/2.0 - (int) (start_pos[0] + volSize[0]/2.0) ) );
                else t[0] = (-1.0 / ray_dir[0]) * ( start_pos[0] + volSize[0]/2.0 - (int) (start_pos[0] + volSize[0]/2.0) );
                if (ray_sig[1]>0) t[1] = (1.0 / ray_dir[1]) * (1.0 - ( start_pos[1] + volSize[1]/2.0 - (int) (start_pos[1] + volSize[1]/2.0) ) );
                else t[1] = (-1.0 / ray_dir[1]) * ( start_pos[1] + volSize[1]/2.0 - (int) (start_pos[1] + volSize[1]/2.0) );
                }
            }
        if ( boundary_hit_pos[2][1][0]>-volSize[0]/2.0 && boundary_hit_pos[2][1][0]<volSize[0]/2.0 && boundary_hit_pos[2][1][1]>-volSize[1]/2.0 && boundary_hit_pos[2][1][1]<volSize[1]/2.0 ) {
            hit_count++;
            if (boundary_hit[2][1]<current_k) {
                current_k = boundary_hit[2][1];
                start_index[2] = 0;
                start_index[1] = (int) (boundary_hit_pos[2][1][1]+volSize[1]/2.0);
                start_index[0] = (int) (boundary_hit_pos[2][1][1]+volSize[2]/2.0);
                start_pos[0] = boundary_hit_pos[2][1][0];
                start_pos[1] = boundary_hit_pos[2][1][1];
                start_pos[2] = boundary_hit_pos[2][1][2];
                t[2] = ((float) ray_sig[2]) / ray_dir[2];
                if (ray_sig[0]>0) t[0] = (1.0 / ray_dir[0]) * (1.0 - ( start_pos[0] + volSize[0]/2.0 - (int) (start_pos[0] + volSize[0]/2.0) ) );
                else t[0] = (-1.0 / ray_dir[0]) * ( start_pos[0] + volSize[0]/2.0 - (int) (start_pos[0] + volSize[0]/2.0) );
                if (ray_sig[1]>0) t[1] = (1.0 / ray_dir[1]) * (1.0 - ( start_pos[1] + volSize[1]/2.0 - (int) (start_pos[1] + volSize[1]/2.0) ) );
                else t[1] = (-1.0 / ray_dir[1]) * ( start_pos[1] + volSize[1]/2.0 - (int) (start_pos[1] + volSize[1]/2.0) );
                }
            }
        if ( hit_count < 2 ) continue;

        float ray_length = 0.0;
        int X = start_index[0];
        int Y = start_index[1];
        int Z = start_index[2];

        raySum = 0.0;

        // the voxel ray tracing looks quite nice and efficient:
        while (X>=0 && X<volSize[0] && Y>=0 && Y<volSize[1] && Z>=0 && Z<volSize[2]) {
            if (t[0]<t[1]) {
                if (t[2]<t[0]) {
                    // 2=min:
                    t[0] -= t[2]; t[1] -= t[2];
                    raySum += t[2]*vol[Z*volSize[0]*volSize[1]+Y*volSize[0]+X];
                    //proj[indY*projSize[0]+indX] += t[2]*vol[Z*volSize[0]*volSize[1]+Y*volSize[0]+X];
                    //vol[Z*volSize[0]*volSize[1]+Y*volSize[0]+X] += t[2]*proj[indY*projSize[0]+indX];
                    ray_length += t[2];
                    t[2] = ((float) ray_sig[2]) / ray_dir[2];
                    Z += ray_sig[2];
                } else {
                    // 0=min:
                    t[1] -= t[0]; t[2] -= t[0];
                    raySum += t[0]*vol[Z*volSize[0]*volSize[1]+Y*volSize[0]+X];
                    //proj[indY*projSize[0]+indX] += t[0]*vol[Z*volSize[0]*volSize[1]+Y*volSize[0]+X];
                    //vol[Z*volSize[0]*volSize[1]+Y*volSize[0]+X] += t[0]*proj[indY*projSize[0]+indX];
                    ray_length += t[0];
                    t[0] = ((float) ray_sig[0]) / ray_dir[0];
                    X += ray_sig[0];
                }
            } else {
                if (t[2]<t[1]) {
                    // 2=min:
                    t[0] -= t[2]; t[1] -= t[2];
                    raySum += t[2]*vol[Z*volSize[0]*volSize[1]+Y*volSize[0]+X];
                    //proj[indY*projSize[0]+indX] += t[2]*vol[Z*volSize[0]*volSize[1]+Y*volSize[0]+X];
                    //vol[Z*volSize[0]*volSize[1]+Y*volSize[0]+X] += t[2]*proj[indY*projSize[0]+indX];
                    ray_length += t[2];
                    t[2] = ((float) ray_sig[2]) / ray_dir[2];
                    Z += ray_sig[2];
                } else {
                    // 1=min:
                    t[0] -= t[1]; t[2] -= t[1];
                    raySum += t[1]*vol[Z*volSize[0]*volSize[1]+Y*volSize[0]+X];
                    //proj[indY*projSize[0]+indX] += t[1]*vol[Z*volSize[0]*volSize[1]+Y*volSize[0]+X];
                    //vol[Z*volSize[0]*volSize[1]+Y*volSize[0]+X] += t[1]*proj[indY*projSize[0]+indX];
                    ray_length += t[1];
                    t[1] = ((float) ray_sig[1]) / ray_dir[1];
                    Y += ray_sig[1];
                }
            }
        }
        proj[indY*projSize[0]+indX] -= raySum;
        proj[indY*projSize[0]+indX] /= ray_length; 
    }
}
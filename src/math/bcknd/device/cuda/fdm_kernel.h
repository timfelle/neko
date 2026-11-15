#ifndef __MATH_FDM_KERNEL_H__
#define __MATH_FDM_KERNEL_H__
/*
 Copyright (c) 2021-2022, The Neko Authors
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.

   * Redistributions in binary form must reproduce the above
     copyright notice, this list of conditions and the following
     disclaimer in the documentation and/or other materials provided
     with the distribution.

   * Neither the name of the authors nor the names of its
     contributors may be used to endorse or promote products derived
     from this software without specific prior written permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE.
*/

template< typename T, const int NL >
__global__ void fdm_do_fast_kernel(T  * __restrict__ e,
                                   T * __restrict__ r,
                                   T * __restrict__ s,
                                   T * __restrict__ d) {
  __shared__ T shwork[NL*NL*NL];
  __shared__ T shwork2[NL*NL*NL];
  __shared__ T A[NL*NL];
  __shared__ T Bt[NL*NL];
  __shared__ T Ct[NL*NL];

  const int idx = threadIdx.x;
  const int str = blockDim.x;
  const int el = blockIdx.x;
  if( idx < NL*NL){
     A[idx] = s[idx+NL*NL+el*NL*NL*3*2];
    Bt[idx] = s[idx+2*NL*NL+el*NL*NL*3*2];
    Ct[idx] = s[idx+2*2*NL*NL+el*NL*NL*3*2];
  }    
  __syncthreads();

  for (int ii = idx; ii< NL*NL*NL; ii += str){
    T tmp = 0.0;
    int j = ii/NL;
    int i = ii - j*NL;
    for( int l = 0; l < NL; l++){
      tmp += A[i+l*NL]*r[l+NL*j+el*NL*NL*NL];
    }
    shwork[ii] = tmp;
  }
  __syncthreads();
  for (int ijk = idx; ijk< NL*NL*NL; ijk += str){
    const int jk = ijk / NL;
    const int i = ijk - jk * NL;
    const int k = jk / NL;
    const int j = jk - k * NL;
    T tmp = 0.0;
    const int ik2 = i + k*NL*NL; 
    for( int l = 0; l < NL; l++){
      tmp += Bt[l+j*NL]*shwork[l*NL+ik2];
    }
    shwork2[ijk] = tmp;
  }
  __syncthreads();
  for (int ijk = idx; ijk< NL*NL*NL; ijk += str){
    const int jk = ijk / NL;
    const int i = ijk - jk * NL;
    const int k = jk / NL;
    const int j = jk - k * NL;
    T tmp = 0.0;
    const int ij2 = i + j*NL; 
    for( int l = 0; l < NL; l++){
      tmp += Ct[l+k*NL]*shwork2[ij2 + l*NL*NL];
    }
    r[ijk+el*NL*NL*NL] = tmp*d[ijk+el*NL*NL*NL];
  }
  __syncthreads();
  if( idx < NL*NL){
     A[idx] = s[idx+el*NL*NL*3*2];
    Bt[idx] = s[idx+NL*NL+2*NL*NL+el*NL*NL*3*2];
    Ct[idx] = s[idx+NL*NL+2*2*NL*NL+el*NL*NL*3*2];
  }  
  __syncthreads();

  for (int ii = idx; ii< NL*NL*NL; ii += str){
    T tmp = 0.0;
    int j = ii/NL;
    int i = ii - j*NL;
    for( int l = 0; l < NL; l++){
      tmp += A[i+l*NL]*r[l+NL*j+el*NL*NL*NL];
    }
    shwork[ii] = tmp;
  }
  __syncthreads();
  for (int ijk = idx; ijk< NL*NL*NL; ijk += str){
    const int jk = ijk / NL;
    const int i = ijk - jk * NL;
    const int k = jk / NL;
    const int j = jk - k * NL;
    T tmp = 0.0;
    const int ik2 = i + k*NL*NL; 
    for( int l = 0; l < NL; l++){
      tmp += Bt[l+j*NL]*shwork[l*NL+ik2];
    }
    shwork2[ijk] = tmp;
  }
  __syncthreads();
  for (int ijk = idx; ijk< NL*NL*NL; ijk += str){
    const int jk = ijk / NL;
    const int i = ijk - jk * NL;
    const int k = jk / NL;
    const int j = jk - k * NL;
    T tmp = 0.0;
    const int ij2 = i + j*NL; 
    for( int l = 0; l < NL; l++){
      tmp += Ct[l+k*NL]*shwork2[ij2 + l*NL*NL];
    }
    e[ijk+el*NL*NL*NL] = tmp;
  }

}

/**
 * Contractions along r and s of one k plane, for fdm_do_fast_kernel_kstep
 *
 * Every thread of the (NL, NL) plane passes its value @a u of the plane and
 * gets back its value of the contracted plane. Both barriers are inside, so
 * the shared planes can be reused by the next call straight away.
 */
template< typename T, const int NL >
__device__ __forceinline__ T fdm_kstep_rs(const T * A, const T * Bt,
                                          T * shu, T * shw, const T u,
                                          const int i, const int j) {
  shu[i+j*NL] = u;
  __syncthreads();
  T tmp = 0.0;
#pragma unroll
  for( int l = 0; l < NL; l++){
    tmp += A[i+l*NL]*shu[l+j*NL];
  }
  shw[i+j*NL] = tmp;
  __syncthreads();
  tmp = 0.0;
#pragma unroll
  for( int l = 0; l < NL; l++){
    tmp += Bt[l+j*NL]*shw[i+l*NL];
  }
  return tmp;
}

/**
 * Fast diagonalisation solve for large NL
 *
 * fdm_do_fast_kernel keeps two whole NL^3 work arrays in shared memory,
 * which with its three matrices is more than the 48 kB a block can have
 * statically from NL = 15 in double precision. A block here is one (NL, NL)
 * thread plane per element, as in the kstep operator kernels. The
 * contractions along r and s are done one k plane at a time through two
 * NL^2 planes of shared memory, and the contraction along t in registers,
 * along each thread's column. Every sum runs in the same order as in
 * fdm_do_fast_kernel, and r is likewise left holding the scaled result of
 * the first solve.
 */
template< typename T, const int NL >
__global__ void __launch_bounds__(NL*NL)
fdm_do_fast_kernel_kstep(T * __restrict__ e,
                         T * __restrict__ r,
                         const T * __restrict__ s,
                         const T * __restrict__ d) {
  __shared__ T A[NL*NL];
  __shared__ T Bt[NL*NL];
  __shared__ T Ct[NL*NL];
  __shared__ T shu[NL*NL];
  __shared__ T shw[NL*NL];

  const int i = threadIdx.x;
  const int j = threadIdx.y;
  const int ij = i + j*NL;
  const int ele = blockIdx.x*NL*NL*NL;
  const int els = blockIdx.x*NL*NL*3*2;

  /* This thread's column after the r and s contractions, and of the scaled
     result of the first solve, which is the input of the second */
  T rs[NL];
  T res[NL];

  A[ij] = s[ij+NL*NL+els];
  Bt[ij] = s[ij+2*NL*NL+els];
  Ct[ij] = s[ij+2*2*NL*NL+els];

#pragma unroll
  for (int k = 0; k < NL; k++) {
    rs[k] = fdm_kstep_rs<T, NL>(A, Bt, shu, shw, r[ij+k*NL*NL+ele], i, j);
  }

#pragma unroll
  for (int k = 0; k < NL; k++) {
    T tmp = 0.0;
#pragma unroll
    for( int l = 0; l < NL; l++){
      tmp += Ct[l+k*NL]*rs[l];
    }
    res[k] = tmp*d[ij+k*NL*NL+ele];
    r[ij+k*NL*NL+ele] = res[k];
  }

  /* Every thread is done with the matrices of the first solve */
  __syncthreads();
  A[ij] = s[ij+els];
  Bt[ij] = s[ij+NL*NL+2*NL*NL+els];
  Ct[ij] = s[ij+NL*NL+2*2*NL*NL+els];

#pragma unroll
  for (int k = 0; k < NL; k++) {
    rs[k] = fdm_kstep_rs<T, NL>(A, Bt, shu, shw, res[k], i, j);
  }

#pragma unroll
  for (int k = 0; k < NL; k++) {
    T tmp = 0.0;
#pragma unroll
    for( int l = 0; l < NL; l++){
      tmp += Ct[l+k*NL]*rs[l];
    }
    e[ij+k*NL*NL+ele] = tmp;
  }
}


#endif // __MATH_FDM_KERNEL_H__

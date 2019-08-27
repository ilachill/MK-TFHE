#ifndef MKTFHESAMPLES_H
#define MKTFHESAMPLES_H

#include "tfhe_core.h"



// MK LWE sample (a_1, ..., a_k, b)
struct MKLweSample {
	Torus32* a; //-- the parties*n coefs of the mask
    Torus32 b;  //
   	double current_variance; //-- average noise of the sample
   	const int32_t parties;
   	const int32_t n;

#ifdef __cplusplus
   MKLweSample(const LweParams* LWEparams, const MKTFHEParams* MKparams);
   ~MKLweSample();
   MKLweSample(const MKLweSample&)=delete;
   MKLweSample& operator=(const MKLweSample&)=delete;
#endif
};


// alloc 
EXPORT MKLweSample* alloc_MKLweSample();
EXPORT MKLweSample* alloc_MKLweSample_array(int32_t nbelts);
//free memory space 
EXPORT void free_MKLweSample(MKLweSample* ptr);
EXPORT void free_MKLweSample_array(int32_t nbelts, MKLweSample* ptr);
// init
EXPORT void init_MKLweSample(MKLweSample* obj, const LweParams* LWEparams, const MKTFHEParams* MKparams);
EXPORT void init_MKLweSample_array(int32_t nbelts, MKLweSample* obj, const LweParams* LWEparams, 
        const MKTFHEParams* MKparams);
// destroys the structure
EXPORT void destroy_MKLweSample(MKLweSample* obj);
EXPORT void destroy_MKLweSample_array(int32_t nbelts, MKLweSample* obj);
// new = alloc + init
EXPORT MKLweSample* new_MKLweSample(const LweParams* LWEparams, const MKTFHEParams* MKparams);
EXPORT MKLweSample* new_MKLweSample_array(int32_t nbelts, const LweParams* LWEparams, const MKTFHEParams* MKparams);
// delete = destroy + free
EXPORT void delete_MKLweSample(MKLweSample* obj);
EXPORT void delete_MKLweSample_array(int32_t nbelts, MKLweSample* obj);












// MK RLWE sample (a_1, ..., a_k, b)
struct MKTLweSample {
    TorusPolynomial *a; ///< array of length parties+1: mask + right term
    TorusPolynomial *b; ///< alias of a[parties] to get the right term
    double current_variance; ///< avg variance of the sample
    const int32_t parties;
    const int32_t N;


#ifdef __cplusplus
    MKTLweSample(const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
    ~MKTLweSample();
    MKTLweSample(const MKTLweSample &) = delete;
    void operator=(const MKTLweSample &) = delete;
#endif
};


// alloc 
EXPORT MKTLweSample* alloc_MKTLweSample();
EXPORT MKTLweSample* alloc_MKTLweSample_array(int32_t nbelts);
//free memory space 
EXPORT void free_MKTLweSample(MKTLweSample* ptr);
EXPORT void free_MKTLweSample_array(int32_t nbelts, MKTLweSample* ptr);
// init
EXPORT void init_MKTLweSample(MKTLweSample* obj, const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
EXPORT void init_MKTLweSample_array(int32_t nbelts, MKTLweSample* obj, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams);
// destroys the structure
EXPORT void destroy_MKTLweSample(MKTLweSample* obj);
EXPORT void destroy_MKTLweSample_array(int32_t nbelts, MKTLweSample* obj);
// new = alloc + init
EXPORT MKTLweSample* new_MKTLweSample(const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
EXPORT MKTLweSample* new_MKTLweSample_array(int32_t nbelts, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams);
// delete = destroy + free
EXPORT void delete_MKTLweSample(MKTLweSample* obj);
EXPORT void delete_MKTLweSample_array(int32_t nbelts, MKTLweSample* obj);

















// MK RLWE sample FFT (a_1, ..., a_k, b)
struct MKTLweSampleFFT {
    LagrangeHalfCPolynomial *a; ///< array of length parties+1: mask + right term
    LagrangeHalfCPolynomial *b; ///< alias of a[parties] to get the right term
    double current_variance; ///< avg variance of the sample
    const int32_t parties; 

#ifdef __cplusplus
    MKTLweSampleFFT(const TLweParams *params, const MKTFHEParams* MKparams, LagrangeHalfCPolynomial *arr, 
		double current_variance);
    ~MKTLweSampleFFT();
    MKTLweSampleFFT(const MKTLweSampleFFT &) = delete;
    void operator=(const MKTLweSampleFFT &) = delete;
#endif
};


// alloc 
EXPORT MKTLweSampleFFT* alloc_MKTLweSampleFFT();
EXPORT MKTLweSampleFFT* alloc_MKTLweSampleFFT_array(int32_t nbelts);
//free memory space 
EXPORT void free_MKTLweSampleFFT(MKTLweSampleFFT* ptr);
EXPORT void free_MKTLweSampleFFT_array(int32_t nbelts, MKTLweSampleFFT* ptr);
// init
EXPORT void init_MKTLweSampleFFT(MKTLweSampleFFT* obj, const TLweParams* RLWEparams, const MKTFHEParams* MKparams, 
        LagrangeHalfCPolynomial *arr, double current_variance);
EXPORT void init_MKTLweSampleFFT_array(int32_t nbelts, MKTLweSampleFFT* obj, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams, LagrangeHalfCPolynomial *arr, double current_variance);
// destroys the structure
EXPORT void destroy_MKTLweSampleFFT(MKTLweSampleFFT* obj);
EXPORT void destroy_MKTLweSampleFFT_array(int32_t nbelts, MKTLweSampleFFT* obj);
// new = alloc + init
EXPORT MKTLweSampleFFT* new_MKTLweSampleFFT(const TLweParams* RLWEparams, const MKTFHEParams* MKparams, 
        LagrangeHalfCPolynomial *arr, double current_variance);
EXPORT MKTLweSampleFFT* new_MKTLweSampleFFT_array(int32_t nbelts, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams, LagrangeHalfCPolynomial *arr, double current_variance);
// delete = destroy + free
EXPORT void delete_MKTLweSampleFFT(MKTLweSampleFFT* obj);
EXPORT void delete_MKTLweSampleFFT_array(int32_t nbelts, MKTLweSampleFFT* obj);




























































/* **************************************************************************
***************************** VERSION 2 *************************************
************************************************************************** */



// MK RGSW UniEnc sample (d,F)=(d,f0,f1)
struct MKTGswUESample_v2 {
    TorusPolynomial *d; ///< array of length 3*dg
    TorusPolynomial *f0; ///< alias of d[dg]
    TorusPolynomial *f1; ///< alias of d[2*dg]
    int32_t party; // party 
    double current_variance; ///< avg variance of the sample
    const int32_t dg;
    const int32_t N;

#ifdef __cplusplus
    MKTGswUESample_v2(const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
    ~MKTGswUESample_v2();
    MKTGswUESample_v2(const MKTGswUESample_v2 &) = delete;
    void operator=(const MKTGswUESample_v2 &) = delete;
#endif
};

// alloc
EXPORT MKTGswUESample_v2* alloc_MKTGswUESample_v2();
EXPORT MKTGswUESample_v2* alloc_MKTGswUESample_v2_array(int32_t nbelts);
//free memory space 
EXPORT void free_MKTGswUESample_v2(MKTGswUESample_v2* ptr);
EXPORT void free_MKTGswUESample_v2_array(int32_t nbelts, MKTGswUESample_v2* ptr);
// initialize the structure
EXPORT void init_MKTGswUESample_v2(MKTGswUESample_v2* obj, const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
EXPORT void init_MKTGswUESample_v2_array(int32_t nbelts, MKTGswUESample_v2* obj, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams);
//destroys the structure
EXPORT void destroy_MKTGswUESample_v2(MKTGswUESample_v2* obj);
EXPORT void destroy_MKTGswUESample_v2_array(int32_t nbelts, MKTGswUESample_v2* obj);
// new = alloc + init
EXPORT MKTGswUESample_v2* new_MKTGswUESample_v2(const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
EXPORT MKTGswUESample_v2* new_MKTGswUESample_v2_array(int32_t nbelts, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams);
// delete = destroy + free
EXPORT void delete_MKTGswUESample_v2(MKTGswUESample_v2* obj);
EXPORT void delete_MKTGswUESample_v2_array(int32_t nbelts, MKTGswUESample_v2* obj);










// MK RGSW UniEnc sample FFT (d,F)=(d1,f0,f1)
struct MKTGswUESampleFFT_v2 {
    LagrangeHalfCPolynomial *d; ///< array of length 3*dg
    LagrangeHalfCPolynomial *f0; ///< alias of d[dg]
    LagrangeHalfCPolynomial *f1; ///< alias of d[2*dg]
    int32_t party; // party 
    double current_variance; ///< avg variance of the sample
    const int32_t dg; 

#ifdef __cplusplus
    MKTGswUESampleFFT_v2(const TLweParams* RLWEparams, const MKTFHEParams* MKparams, 
        LagrangeHalfCPolynomial *arr, int32_t party, double current_variance);
    ~MKTGswUESampleFFT_v2();
    MKTGswUESampleFFT_v2(const MKTGswUESampleFFT_v2 &) = delete;
    void operator=(const MKTGswUESampleFFT_v2 &) = delete;
#endif
};

// alloc
EXPORT MKTGswUESampleFFT_v2* alloc_MKTGswUESampleFFT_v2();
EXPORT MKTGswUESampleFFT_v2* alloc_MKTGswUESampleFFT_v2_array(int32_t nbelts);
//free memory space 
EXPORT void free_MKTGswUESampleFFT_v2(MKTGswUESampleFFT_v2* ptr);
EXPORT void free_MKTGswUESampleFFT_v2_array(int32_t nbelts, MKTGswUESampleFFT_v2* ptr);
// initialize the structure
EXPORT void init_MKTGswUESampleFFT_v2(MKTGswUESampleFFT_v2* obj, const TLweParams* RLWEparams, const MKTFHEParams* MKparams, 
        int32_t party, double current_variance);
EXPORT void init_MKTGswUESampleFFT_v2_array(int32_t nbelts, MKTGswUESampleFFT_v2* obj, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams, int32_t party, double current_variance);
//destroys the structure
EXPORT void destroy_MKTGswUESampleFFT_v2(MKTGswUESampleFFT_v2* obj);
EXPORT void destroy_MKTGswUESampleFFT_v2_array(int32_t nbelts, MKTGswUESampleFFT_v2* obj);
// new = alloc + init
EXPORT MKTGswUESampleFFT_v2* new_MKTGswUESampleFFT_v2(const TLweParams* RLWEparams, const MKTFHEParams* MKparams, 
        int32_t party, double current_variance);
EXPORT MKTGswUESampleFFT_v2* new_MKTGswUESampleFFT_v2_array(int32_t nbelts, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams, int32_t party, double current_variance);
// delete = destroy + free
EXPORT void delete_MKTGswUESampleFFT_v2(MKTGswUESampleFFT_v2* obj);
EXPORT void delete_MKTGswUESampleFFT_v2_array(int32_t nbelts, MKTGswUESampleFFT_v2* obj);










// MK RGSW Expanded sample: party i D_i=(x_0, ..., x_{parties-1}, x_parties + d_i, y_0, ..., d_i+y_i, ..., y_perties, d_i)
struct MKTGswExpSample_v2 {
    TorusPolynomial *x; ///< array of length (2*(parties+1)+1)*dg
    TorusPolynomial *y;
    TorusPolynomial *d;
    int32_t party; // party (from 0 to parties-1)
    double current_variance; ///< avg variance of the sample
    const int32_t parties; 
    const int32_t dg;
    const int32_t N;

#ifdef __cplusplus
    MKTGswExpSample_v2(const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
    ~MKTGswExpSample_v2();
    MKTGswExpSample_v2(const MKTGswExpSample_v2 &) = delete;
    void operator=(const MKTGswExpSample_v2 &) = delete;
#endif
};

// alloc
EXPORT MKTGswExpSample_v2* alloc_MKTGswExpSample_v2() ;
EXPORT MKTGswExpSample_v2* alloc_MKTGswExpSample_v2_array(int32_t nbelts);
//free memory space 
EXPORT void free_MKTGswExpSample_v2(MKTGswExpSample_v2* ptr);
EXPORT void free_MKTGswExpSample_v2_array(int32_t nbelts, MKTGswExpSample_v2* ptr);
// initialize the structure
EXPORT void init_MKTGswExpSample_v2(MKTGswExpSample_v2* obj, const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
EXPORT void init_MKTGswExpSample_v2_array(int32_t nbelts, MKTGswExpSample_v2* obj, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams);
//destroys the structure
EXPORT void destroy_MKTGswExpSample_v2(MKTGswExpSample_v2* obj);
EXPORT void destroy_MKTGswExpSample_v2_array(int32_t nbelts, MKTGswExpSample_v2* obj);
// new = alloc + init
EXPORT MKTGswExpSample_v2* new_MKTGswExpSample_v2(const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
EXPORT MKTGswExpSample_v2* new_MKTGswExpSample_v2_array(int32_t nbelts, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams);
// delete = destroy + free
EXPORT void delete_MKTGswExpSample_v2(MKTGswExpSample_v2* obj);
EXPORT void delete_MKTGswExpSample_v2_array(int32_t nbelts, MKTGswExpSample_v2* obj);









// MK RGSW Expanded sample FFT: party i D_i=(x_0, ..., x_{parties-1}, x_parties + d_i, y_0, ..., d_i+y_i, ..., y_perties, d_i)
struct MKTGswExpSampleFFT_v2 {
    LagrangeHalfCPolynomial *x; ///< array of length (2*(parties+1)+1)*dg
    LagrangeHalfCPolynomial *y;
    LagrangeHalfCPolynomial *d;
    int32_t party; // party (from 0 to parties-1)
    double current_variance; ///< avg variance of the sample
    const int32_t parties;
    const int32_t dg;
    
#ifdef __cplusplus
    MKTGswExpSampleFFT_v2(const TLweParams* RLWEparams, const MKTFHEParams* MKparams, 
        LagrangeHalfCPolynomial *arr, double current_variance);
    ~MKTGswExpSampleFFT_v2();
    MKTGswExpSampleFFT_v2(const MKTGswExpSampleFFT_v2 &) = delete;
    void operator=(const MKTGswExpSampleFFT_v2 &) = delete;
#endif
};

// alloc
EXPORT MKTGswExpSampleFFT_v2* alloc_MKTGswExpSampleFFT_v2();
EXPORT MKTGswExpSampleFFT_v2* alloc_MKTGswExpSampleFFT_v2_array(int32_t nbelts);
//free memory space 
EXPORT void free_MKTGswExpSampleFFT_v2(MKTGswExpSampleFFT_v2* ptr);
EXPORT void free_MKTGswExpSampleFFT_v2_array(int32_t nbelts, MKTGswExpSampleFFT_v2* ptr);
// initialize the structure
EXPORT void init_MKTGswExpSampleFFT_v2(MKTGswExpSampleFFT_v2* obj, const TLweParams* RLWEparams, const MKTFHEParams* MKparams, 
        double current_variance);
EXPORT void init_MKTGswExpSampleFFT_v2_array(int32_t nbelts, MKTGswExpSampleFFT_v2* obj, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams, double current_variance);
//destroys the structure
EXPORT void destroy_MKTGswExpSampleFFT_v2(MKTGswExpSampleFFT_v2* obj);
EXPORT void destroy_MKTGswExpSampleFFT_v2_array(int32_t nbelts, MKTGswExpSampleFFT_v2* obj);
// new = alloc + init
EXPORT MKTGswExpSampleFFT_v2* new_MKTGswExpSampleFFT_v2(const TLweParams* RLWEparams, const MKTFHEParams* MKparams, 
        double current_variance);
EXPORT MKTGswExpSampleFFT_v2* new_MKTGswExpSampleFFT_v2_array(int32_t nbelts, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams, double current_variance);
// delete = destroy + free
EXPORT void delete_MKTGswExpSampleFFT_v2(MKTGswExpSampleFFT_v2* obj);
EXPORT void delete_MKTGswExpSampleFFT_v2_array(int32_t nbelts, MKTGswExpSampleFFT_v2* obj);













#endif //MKTFHESAMPLES_H


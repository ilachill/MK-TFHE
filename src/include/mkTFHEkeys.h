#ifndef MKTFHEKEYS_H
#define MKTFHEKEYS_H


#include "tfhe_core.h"


struct MKLweKey {
   const LweParams* LWEparams;
   const MKTFHEParams* MKparams;

   LweKey* key; // LWE secret keys for all the parties

#ifdef __cplusplus   
   MKLweKey(const LweParams* LWEparams, const MKTFHEParams* MKparams);
   ~MKLweKey();
   MKLweKey(const MKLweKey&) = delete; //forbidden 
   MKLweKey* operator=(const MKLweKey&) = delete; //forbidden
#endif
};



// allocate memory space 
EXPORT MKLweKey* alloc_MKLweKey();
EXPORT MKLweKey* alloc_MKLweKey_array(int32_t nbelts);
// free memory space 
EXPORT void free_MKLweKey(MKLweKey* ptr);
EXPORT void free_MKLweKey_array(int32_t nbelts, MKLweKey* ptr);
// initialize the structure
EXPORT void init_MKLweKey(MKLweKey* obj, const LweParams* LWEparams, const MKTFHEParams* MKparams);
EXPORT void init_MKLweKey_array(int32_t nbelts, MKLweKey* obj, const LweParams* LWEparams, 
        const MKTFHEParams* MKparams);
// destroys the structure
EXPORT void destroy_MKLweKey(MKLweKey* obj);
EXPORT void destroy_MKLweKey_array(int32_t nbelts, MKLweKey* obj);
// new = alloc + init
EXPORT MKLweKey* new_MKLweKey(const LweParams* LWEparams, const MKTFHEParams* MKparams);
EXPORT MKLweKey* new_MKLweKey_array(int32_t nbelts, const LweParams* LWEparams, const MKTFHEParams* MKparams);
// delete = destroy + free
EXPORT void delete_MKLweKey(MKLweKey* obj);
EXPORT void delete_MKLweKey_array(int32_t nbelts, MKLweKey* obj);



















struct MKRLweKey {
    const TLweParams* RLWEparams;
    const MKTFHEParams* MKparams;

    TLweKey* key; // RLWE secret keys for all the parties
    TorusPolynomial* Pkey; // RLWE public keys for all the parties

#ifdef __cplusplus
    MKRLweKey(const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
    ~MKRLweKey();
    MKRLweKey(const MKRLweKey &) = delete;
    MKRLweKey* operator=(const MKRLweKey &) = delete;
#endif
};



// allocate memory space 
EXPORT MKRLweKey* alloc_MKRLweKey();
EXPORT MKRLweKey* alloc_MKRLweKey_array(int32_t nbelts);
// free memory space 
EXPORT void free_MKRLweKey(MKRLweKey* ptr);
EXPORT void free_MKRLweKey_array(int32_t nbelts, MKRLweKey* ptr);
// initialize the structure
EXPORT void init_MKRLweKey(MKRLweKey* obj, const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
EXPORT void init_MKRLweKey_array(int32_t nbelts, MKRLweKey* obj, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams);
// destroys the structure
EXPORT void destroy_MKRLweKey(MKRLweKey* obj);
EXPORT void destroy_MKRLweKey_array(int32_t nbelts, MKRLweKey* obj);
// new = alloc + init
EXPORT MKRLweKey* new_MKRLweKey(const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
EXPORT MKRLweKey* new_MKRLweKey_array(int32_t nbelts, const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
// delete = destroy + free
EXPORT void delete_MKRLweKey(MKRLweKey* obj);
EXPORT void delete_MKRLweKey_array(int32_t nbelts, MKRLweKey* obj);


























/* *******************************************************
*************** Key Switching Key ************************
******************************************************* */

struct MKLweKeySwitchKey {
    int32_t n_in;                 // length input key
    int32_t n_out;                // length output key
    int32_t parties;              // number of parties
    int32_t Bksbit;               // KS basebit
    int32_t Bks;                  // KS base
    int32_t dks;                  // KS lenght
    const MKTFHEParams* params;
    MKLweSample* ks0_raw;         // vector of size parties*n_in*dks*Bks
    MKLweSample** ks1_raw;        
    MKLweSample*** ks2_raw;       
    MKLweSample**** ks;           

#ifdef __cplusplus
    MKLweKeySwitchKey(int32_t n_in, const MKTFHEParams* params, MKLweSample* ks0_raw);
    ~MKLweKeySwitchKey();
    MKLweKeySwitchKey(const MKLweKeySwitchKey&) = delete;
    void operator=(const MKLweKeySwitchKey&) = delete;
#endif
};


// alloc 
EXPORT MKLweKeySwitchKey* alloc_MKLweKeySwitchKey();
EXPORT MKLweKeySwitchKey* alloc_MKLweKeySwitchKey_array(int32_t nbelts);
// free memory space 
EXPORT void free_MKLweKeySwitchKey(MKLweKeySwitchKey* ptr);
EXPORT void free_MKLweKeySwitchKey_array(int32_t nbelts, MKLweKeySwitchKey* ptr);
// initialize the structure
EXPORT void init_MKLweKeySwitchKey(MKLweKeySwitchKey* obj, int32_t n_in, const LweParams* LWEparams, 
        const MKTFHEParams* params);
EXPORT void init_MKLweKeySwitchKey_array(int32_t nbelts, MKLweKeySwitchKey* obj, int32_t n_in, 
        const LweParams* LWEparams, const MKTFHEParams* params);
// destroy 
EXPORT void destroy_MKLweKeySwitchKey(MKLweKeySwitchKey* obj);
EXPORT void destroy_MKLweKeySwitchKey_array(int32_t nbelts, MKLweKeySwitchKey* obj);
// new = alloc + init 
EXPORT MKLweKeySwitchKey* new_MKLweKeySwitchKey(int32_t n_in, const LweParams* LWEparams, const MKTFHEParams* params);
EXPORT MKLweKeySwitchKey* new_MKLweKeySwitchKey_array(int32_t nbelts, int32_t n_in, const LweParams* LWEparams, 
        const MKTFHEParams* params);
// delete = destroy + free 
EXPORT void delete_MKLweKeySwitchKey(MKLweKeySwitchKey* obj);
EXPORT void delete_MKLweKeySwitchKey_array(int32_t nbelts, MKLweKeySwitchKey* obj);
















/* *******************************************************
*************** Bootstrapping Key v2 *********************
******************************************************* */


struct MKLweBootstrappingKey_v2{
    const MKTFHEParams* MKparams;
    MKTGswUESample_v2* bk;
    LweKeySwitchKey* ks; //MKLweKeySwitchKey* ks;

#ifdef __cplusplus
   MKLweBootstrappingKey_v2(const MKTFHEParams* MKparams, MKTGswUESample_v2* bk, 
        LweKeySwitchKey* ks);
    ~MKLweBootstrappingKey_v2();
    MKLweBootstrappingKey_v2(const MKLweBootstrappingKey_v2&) = delete;
    void operator=(const MKLweBootstrappingKey_v2&) = delete;
  
#endif
};



// alloc
EXPORT MKLweBootstrappingKey_v2 *alloc_MKLweBootstrappingKey_v2();
EXPORT MKLweBootstrappingKey_v2 *alloc_MKLweBootstrappingKey_v2_array(int32_t nbelts);
// free memory space 
EXPORT void free_MKLweBootstrappingKey_v2(MKLweBootstrappingKey_v2 *ptr);
EXPORT void free_MKLweBootstrappingKey_v2_array(int32_t nbelts, MKLweBootstrappingKey_v2 *ptr);
//initialize the structure
// in mkTFHEkeygen.h
// init_MKLweBootstrappingKey_v2(MKLweBootstrappingKey_v2 *obj, const int32_t n_in, 
//        const LweParams* LWEparams, const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
EXPORT void init_MKLweBootstrappingKey_v2_array(int32_t nbelts, MKLweBootstrappingKey_v2 *obj,  
        const LweParams* LWEparams, const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
// destroys the structure
// in mkTFHEkeygen.h
// destroy_MKLweBootstrappingKey_v2(MKLweBootstrappingKey_v2 *obj);
EXPORT void destroy_MKLweBootstrappingKey_v2_array(int32_t nbelts, MKLweBootstrappingKey_v2 *obj);
// new = alloc + init
EXPORT MKLweBootstrappingKey_v2 *new_MKLweBootstrappingKey_v2(const LweParams* LWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
EXPORT MKLweBootstrappingKey_v2 *new_MKLweBootstrappingKey_v2_array(int32_t nbelts,
        const LweParams* LWEparams, const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
// delete = destroy + free
EXPORT void delete_MKLweBootstrappingKey_v2(MKLweBootstrappingKey_v2 *obj);
EXPORT void delete_MKLweBootstrappingKey_v2_array(int32_t nbelts, MKLweBootstrappingKey_v2 *obj);











/*
 * MKLweBootstrappingKey is converted to a BootstrappingKeyFFT
 */
struct MKLweBootstrappingKeyFFT_v2 {
    const MKTFHEParams* MKparams; 
    MKTGswUESampleFFT_v2* bkFFT;
    LweKeySwitchKey* ks; //const MKLweKeySwitchKey* ks;

#ifdef __cplusplus
   MKLweBootstrappingKeyFFT_v2(const MKTFHEParams* MKparams, 
        MKTGswUESampleFFT_v2* bkFFT, LweKeySwitchKey* ks);
    ~MKLweBootstrappingKeyFFT_v2();
    MKLweBootstrappingKeyFFT_v2(const MKLweBootstrappingKeyFFT_v2&) = delete;
    void operator=(const MKLweBootstrappingKeyFFT_v2&) = delete;
  
#endif
};


// alloc
EXPORT MKLweBootstrappingKeyFFT_v2 *alloc_MKLweBootstrappingKeyFFT_v2();
EXPORT MKLweBootstrappingKeyFFT_v2 *alloc_MKLweBootstrappingKeyFFT_v2_array(int32_t nbelts);
// free memory space 
EXPORT void free_MKLweBootstrappingKeyFFT_v2(MKLweBootstrappingKeyFFT_v2 *ptr);
EXPORT void free_MKLweBootstrappingKeyFFT_v2_array(int32_t nbelts, MKLweBootstrappingKeyFFT_v2 *ptr);
//initialize the structure
// in mkTFHEkeygen.h
// EXPORT void init_MKLweBootstrappingKeyFFT_v2(MKLweBootstrappingKeyFFT_v2 *obj, const MKLweBootstrappingKey *bk,
//   const LweParams* LWEparams, const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
EXPORT void init_MKLweBootstrappingKeyFFT_v2_array(int32_t nbelts, MKLweBootstrappingKeyFFT_v2 *obj, 
        const MKLweBootstrappingKey_v2 *bk, const LweParams* LWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
// destroys the structure
// in mkTFHEkeygen.h
// EXPORT void destroy_MKLweBootstrappingKeyFFT_v2(MKLweBootstrappingKeyFFT_v2 *obj);
EXPORT void destroy_MKLweBootstrappingKeyFFT_v2_array(int32_t nbelts, MKLweBootstrappingKeyFFT_v2 *obj);
// new = alloc + init
EXPORT MKLweBootstrappingKeyFFT_v2 *new_MKLweBootstrappingKeyFFT_v2(const MKLweBootstrappingKey_v2 *bk,  
        const LweParams* LWEparams, const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
EXPORT MKLweBootstrappingKeyFFT_v2 *new_MKLweBootstrappingKeyFFT_v2_array(int32_t nbelts, const MKLweBootstrappingKey_v2 *bk,  
        const LweParams* LWEparams, const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
// delete = destroy + free
EXPORT void delete_MKLweBootstrappingKeyFFT_v2(MKLweBootstrappingKeyFFT_v2 *obj);
EXPORT void delete_MKLweBootstrappingKeyFFT_v2_array(int32_t nbelts, MKLweBootstrappingKeyFFT_v2 *obj);






#endif //MKTFHEKEYS_H



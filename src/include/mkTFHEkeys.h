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
*************** Bootstrapping Key ************************
******************************************************* */


struct MKLweBootstrappingKey{
    const MKTFHEParams* MKparams;
    MKTGswExpSample* bk;
    LweKeySwitchKey* ks; //MKLweKeySwitchKey* ks;

#ifdef __cplusplus
   MKLweBootstrappingKey(const MKTFHEParams* MKparams, MKTGswExpSample* bk, 
        LweKeySwitchKey* ks);
    ~MKLweBootstrappingKey();
    MKLweBootstrappingKey(const MKLweBootstrappingKey&) = delete;
    void operator=(const MKLweBootstrappingKey&) = delete;
  
#endif
};



// alloc
EXPORT MKLweBootstrappingKey *alloc_MKLweBootstrappingKey();
EXPORT MKLweBootstrappingKey *alloc_MKLweBootstrappingKey_array(int32_t nbelts);
// free memory space 
EXPORT void free_MKLweBootstrappingKey(MKLweBootstrappingKey *ptr);
EXPORT void free_MKLweBootstrappingKey_array(int32_t nbelts, MKLweBootstrappingKey *ptr);
//initialize the structure
// in mkTFHEkeygen.h
// init_MKLweBootstrappingKey(MKLweBootstrappingKey *obj, const int32_t n_in, 
//        const LweParams* LWEparams, const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
EXPORT void init_MKLweBootstrappingKey_array(int32_t nbelts, MKLweBootstrappingKey *obj,  
        const LweParams* LWEparams, const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
// destroys the structure
// in mkTFHEkeygen.h
// destroy_MKLweBootstrappingKey(MKLweBootstrappingKey *obj);
EXPORT void destroy_MKLweBootstrappingKey_array(int32_t nbelts, MKLweBootstrappingKey *obj);
// new = alloc + init
EXPORT MKLweBootstrappingKey *new_MKLweBootstrappingKey(const LweParams* LWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
EXPORT MKLweBootstrappingKey *new_MKLweBootstrappingKey_array(int32_t nbelts,
        const LweParams* LWEparams, const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
// delete = destroy + free
EXPORT void delete_MKLweBootstrappingKey(MKLweBootstrappingKey *obj);
EXPORT void delete_MKLweBootstrappingKey_array(int32_t nbelts, MKLweBootstrappingKey *obj);



























/*
 * MKLweBootstrappingKey is converted to a BootstrappingKeyFFT
 */
struct MKLweBootstrappingKeyFFT {
    const MKTFHEParams* MKparams; 
    MKTGswExpSampleFFT* bkFFT;
    LweKeySwitchKey* ks; //const MKLweKeySwitchKey* ks;

#ifdef __cplusplus
   MKLweBootstrappingKeyFFT(const MKTFHEParams* MKparams, 
        MKTGswExpSampleFFT* bkFFT, LweKeySwitchKey* ks);
    ~MKLweBootstrappingKeyFFT();
    MKLweBootstrappingKeyFFT(const MKLweBootstrappingKeyFFT&) = delete;
    void operator=(const MKLweBootstrappingKeyFFT&) = delete;
  
#endif
};


// alloc
EXPORT MKLweBootstrappingKeyFFT *alloc_MKLweBootstrappingKeyFFT();
EXPORT MKLweBootstrappingKeyFFT *alloc_MKLweBootstrappingKeyFFT_array(int32_t nbelts);
// free memory space 
EXPORT void free_MKLweBootstrappingKeyFFT(MKLweBootstrappingKeyFFT *ptr);
EXPORT void free_MKLweBootstrappingKeyFFT_array(int32_t nbelts, MKLweBootstrappingKeyFFT *ptr);
//initialize the structure
// in mkTFHEkeygen.h
// EXPORT void init_MKLweBootstrappingKeyFFT(MKLweBootstrappingKeyFFT *obj, const MKLweBootstrappingKey *bk,
//   const LweParams* LWEparams, const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
EXPORT void init_MKLweBootstrappingKeyFFT_array(int32_t nbelts, MKLweBootstrappingKeyFFT *obj, 
        const MKLweBootstrappingKey *bk, const LweParams* LWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
// destroys the structure
// in mkTFHEkeygen.h
// EXPORT void destroy_MKLweBootstrappingKeyFFT(MKLweBootstrappingKeyFFT *obj);
EXPORT void destroy_MKLweBootstrappingKeyFFT_array(int32_t nbelts, MKLweBootstrappingKeyFFT *obj);
// new = alloc + init
EXPORT MKLweBootstrappingKeyFFT *new_MKLweBootstrappingKeyFFT(const MKLweBootstrappingKey *bk,  
        const LweParams* LWEparams, const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
EXPORT MKLweBootstrappingKeyFFT *new_MKLweBootstrappingKeyFFT_array(int32_t nbelts, const MKLweBootstrappingKey *bk,  
        const LweParams* LWEparams, const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
// delete = destroy + free
EXPORT void delete_MKLweBootstrappingKeyFFT(MKLweBootstrappingKeyFFT *obj);
EXPORT void delete_MKLweBootstrappingKeyFFT_array(int32_t nbelts, MKLweBootstrappingKeyFFT *obj);



#endif //MKTFHEKEYS_H



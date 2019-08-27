#ifndef MKTFHEKEYGEN_H
#define MKTFHEKEYGEN_H

#include "tfhe_core.h"


/**
 * This functions generate keys for the given parameters.
 * The keys for the result must be allocated and initialized
 * (this means that the parameters are already in the result)
 */

// MKLwe 
// key generation for every party 
EXPORT void MKLweKeyGen(MKLweKey* result);

// MKRLwe
// key generation for every party
// secret and public keys
EXPORT void MKRLweKeyGen(MKRLweKey *result);



//extractions Ring Lwe -> Lwe (extracted)
EXPORT void MKtLweExtractKey(MKLweKey* LWEkey, const MKRLweKey* RLWEkey); 







/* *******************************************************
*************** Key Switching Key ************************
******************************************************* */


/*
Create the key switching key: normalize the error in the beginning
 * chose a random vector of gaussian noises (same size as ks) 
 * recenter the noises 
 * generate the ks by creating noiseless encryprions and then add the noise
*/
EXPORT void MKlweCreateKeySwitchKey(MKLweKeySwitchKey* result, const MKLweKey* in_key, const MKLweKey* out_key,
        const MKTFHEParams* MKparams);




















/* *******************************************************
*************** Bootstrapping Key v2 *********************
******************************************************* */


EXPORT void init_MKLweBootstrappingKey_v2(MKLweBootstrappingKey_v2 *obj,
        const LweParams* LWEparams, const TLweParams* RLWEparams, const MKTFHEParams* MKparams);
EXPORT void destroy_MKLweBootstrappingKey_v2(MKLweBootstrappingKey_v2 *obj);


EXPORT void MKlweCreateBootstrappingKey_v2(MKLweBootstrappingKey_v2* result, const MKLweKey* LWEkey, 
        const MKRLweKey* RLWEkey, const MKLweKey* extractedLWEkey, const LweParams *extractedLWEparams,
        const LweParams *LWEparams, const TLweParams *RLWEparams, const MKTFHEParams* MKparams);

// FFT
EXPORT void init_MKLweBootstrappingKeyFFT_v2(MKLweBootstrappingKeyFFT_v2 *obj, 
    const MKLweBootstrappingKey_v2 *bk, const LweParams* LWEparams, const TLweParams* RLWEparams, 
    const MKTFHEParams* MKparams);
EXPORT void destroy_MKLweBootstrappingKeyFFT_v2(MKLweBootstrappingKeyFFT_v2 *obj);


#endif //MKTFHEKEYGEN_H

#ifndef MKTFHEFUNCTIONS_H
#define MKTFHEFUNCTIONS_H


#include "tfhe_core.h"
#include "lwekey.h"
#include "lweparams.h"
#include "lwesamples.h"
#include "tlwe.h"


#include "mkTFHEparams.h"
#include "mkTFHEkeys.h"
#include "mkTFHEsamples.h"



/* *************************
******** MK-LWE ************
************************* */

// b = \sum <a_i,s_i> + m + e
// m comes with a scaling factor 
EXPORT void MKlweSymEncrypt(MKLweSample* result, Torus32 message, double alpha, const MKLweKey* key);

EXPORT void MKlweSymEncryptWithExternalNoise(MKLweSample* result, Torus32 message, double noise, double alpha, 
        const MKLweKey* key);

/**
 * This function computes the phase of sample by using key : phi = b - \sum <a_i,s_i>
 */
EXPORT Torus32 MKlwePhase(const MKLweSample* sample, const MKLweKey* key);

/**
 * This function computes the decryption of sample by using key
 * The constant Msize indicates the message space and is used to approximate the phase
 */
EXPORT Torus32 MKlweSymDecrypt(const MKLweSample* sample, const MKLweKey* key, const int32_t Msize);

/** result = (0, ..., 0, mu) */
EXPORT void MKlweNoiselessTrivial(MKLweSample* result, Torus32 mu, const MKTFHEParams* params);


/** result = result - sample */
EXPORT void MKlweSubTo(MKLweSample* result, const MKLweSample* sample, const MKTFHEParams* MKparams);

/** result = sample */
EXPORT void MKlweCopy(MKLweSample* result, const MKLweSample* sample, const MKTFHEParams* params);













/* *************************
******** MK-RLWE ***********
************************* */


// b = \sum <a_i,s_i> + m + e
// m comes with a scaling factor 
EXPORT void MKtLweSymEncrypt(MKTLweSample *result, TorusPolynomial *message, double alpha, const MKRLweKey *key);

// b = \sum <a_i,s_i> + m + e
// m constant message, comes with a scaling factor 
EXPORT void MKtLweSymEncryptT(MKTLweSample *result, Torus32 message, double alpha, const MKRLweKey *key);

/** result = (0, ..., 0,mu) */
EXPORT void MKtLweNoiselessTrivial(MKTLweSample *result, const TorusPolynomial *mu, const MKTFHEParams *MKparams);


/**
 * This function computes the phase of sample by using key : phi = b - \sum <a_i,s_i>
 */
EXPORT void MKtLwePhase(TorusPolynomial *phase, const MKTLweSample *sample, const MKRLweKey *key);


/**
 * This function computes the decryption of sample by using key
 * The constant Msize indicates the message space and is used to approximate the phase
 */
EXPORT void MKtLweSymDecrypt(TorusPolynomial *result, const MKTLweSample *sample, const MKRLweKey *key, int32_t Msize);


/**
 * This function computes the decryption of sample by using key
 * The constant Msize indicates the message space and is used to approximate the phase
 * The message is a constant torus element
 */
EXPORT Torus32 MKtLweSymDecryptT(const MKTLweSample *sample, const MKRLweKey *key, int32_t Msize);

/** result = sample */
EXPORT void MKtLweCopy(MKTLweSample *result, const MKTLweSample *sample, const MKTFHEParams *MKparams);



// external multiplication of ACC by X^ai-1
EXPORT void MKtLweMulByXaiMinusOne(MKTLweSample *result, int32_t ai, const MKTLweSample *ACC, 
        const MKTFHEParams *MKparams);


/** result = result + sample */
EXPORT void MKtLweAddTo(MKTLweSample *result, const MKTLweSample *sample, const MKTFHEParams *MKparams);



// EXTRACT
EXPORT void MKtLweExtractMKLweSampleIndex(MKLweSample* result, const MKTLweSample* x, const int32_t index, 
        const MKTFHEParams* MKparams) ;
// extract index 0
EXPORT void MKtLweExtractMKLweSample(MKLweSample* result, const MKTLweSample* x, const MKTFHEParams* MKparams);



























/* *************************
******** MK-RGSW ***********
************************* */

/* Uni-Encrypt */
// Encrypt a integer polynomial
EXPORT void MKTGswUniEncrypt(MKTGswUESample *result, IntPolynomial *message, int32_t party, double alpha, const MKRLweKey *key);

// Encrypt an integer value
EXPORT void MKTGswUniEncryptI(MKTGswUESample *result, int32_t message, int32_t party, double alpha, const MKRLweKey *key);


// decrypt (decrypts every line of c and d)
EXPORT void MKtGswSymDecrypt(TorusPolynomial *result, const MKTGswUESample *sample, const MKRLweKey *key);




// same function as tGswTorus32PolynomialDecompH, without the assembly
// (t_0, ..., t_N-1) -> (I_0, ...,I_dg-1)
// decomp_g(t_i) = (I0,j, ..., Idg-1,j)
EXPORT void MKtGswTorus32PolynomialDecompG(IntPolynomial *result, const TorusPolynomial *sample, 
        const MKTFHEParams *params);

EXPORT void MKtGswTorus32PolynomialDecompGassembly(IntPolynomial *result, const TorusPolynomial *sample, 
        const MKTFHEParams *params);


/* EXPAND */
// (C,D,F) = (c0,c1,d0,d1,f0,f1) -> C = (c0, d0+x1, ..., d0, ..., d0+xk, c1, y1, ..., d1, ..., yk)
EXPORT void MKTGswExpand(MKTGswExpSample *result, MKTGswUESample *sample, const MKRLweKey *key, 
		const MKTFHEParams* MKparams);
/**
 * This function is used to verify that the expansion is done correctly
 * C = (y1, ..., d1, ..., yk, c1, d0+x1, ..., d0, ..., d0+xk, c0)
 * The constant Msize indicates the message space and is used to approximate the phase
 */
// result is an array composed by (parties+1)*dg torus polynomials 
EXPORT void MKtGswEXPSymDecrypt(TorusPolynomial *result, MKTGswExpSample *sample, const MKRLweKey *key);


/* EXPAND */
// (C,D,F) = (c0,c1,d0,d1,f0,f1) -> C = (y1, ..., d1, ..., yk, c1, d0+x1, ..., d0, ..., d0+xk, c0)
// sample UE --> resultFFT expand
EXPORT void MKTGswExpandFFT(MKTGswExpSampleFFT *resultFFT, MKTGswUESample *sample, const MKRLweKey *key, 
        const TLweParams* RLWEparams, const MKTFHEParams* MKparams);















/* ********************************************************************************
*************************** EXTERNAL PRODUCT **************************************
******************************************************************************** */

// c' = G^{-1}(c)*C, with C = (y1, ..., d1, ..., yk, c1, d0+x1, ..., d0, ..., d0+xk, c0) 
EXPORT void MKtGswExpExternMulToMKTLwe(MKTLweSample* result, MKTLweSample* sample, const MKTGswExpSample* sampleExp, 
        const MKTFHEParams* MKparams);
// External product FFT
// Result not FFT
EXPORT void MKtGswExpExternMulToMKTLweFFT(MKTLweSample* result, MKTLweSample* sample, 
        const MKTGswExpSampleFFT* sampleExpFFT, const TLweParams* RLWEparams, const MKTFHEParams* MKparams);










/* ********************************************************************************
****************************** KEY SWITCHING **************************************
******************************************************************************** */

EXPORT void MKlweKeySwitch(MKLweSample* result, const LweKeySwitchKey* ks, const MKLweSample* sample, 
        const LweParams* LWEparams, const MKTFHEParams* MKparams);











/* ********************************************************************************
****************************** BOOTSTRAPPING **************************************
******************************************************************************** */

// Blind rotate
EXPORT void MKtfhe_blindRotate(MKTLweSample *accum, const MKTGswExpSample *bk, const int32_t *bara, 
    const TLweParams* RLWEparams, const MKTFHEParams *MKparams);
// Blind rotate
EXPORT void MKtfhe_blindRotateFFT(MKTLweSample *accum, const MKTGswExpSampleFFT *bkFFT, const int32_t *bara, 
    const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey);



EXPORT void MKtfhe_blindRotateAndExtract(MKLweSample *result,
                                       const TorusPolynomial *v,
                                       const MKTGswExpSample *bk,
                                       const int32_t barb,
                                       const int32_t *bara,
                                       const TLweParams* RLWEparams, 
                                       const MKTFHEParams *MKparams);
EXPORT void MKtfhe_blindRotateAndExtractFFT(MKLweSample *result,
                                       const TorusPolynomial *v,
                                       const MKTGswExpSampleFFT *bkFFT,
                                       const int32_t barb,
                                       const int32_t *bara,
                                       const TLweParams* RLWEparams, 
                                       const MKTFHEParams *MKparams,
                                       const MKRLweKey *MKrlwekey); 






// MK bootstrapping without keyswitching
EXPORT void MKtfhe_bootstrap_woKS(MKLweSample *result, const MKLweBootstrappingKey *bk, 
        Torus32 mu, const MKLweSample *x, const TLweParams* RLWEparams, const MKTFHEParams *MKparams);
EXPORT void MKtfhe_bootstrap_woKSFFT(MKLweSample *result, const MKLweBootstrappingKeyFFT *bk, 
        Torus32 mu, const MKLweSample *x, const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey);








// MK bootstrapping 
EXPORT void MKtfhe_bootstrap(MKLweSample *result, const MKLweBootstrappingKey *bk, Torus32 mu, 
        const MKLweSample *x, const LweParams* LWEparams, const LweParams* extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams);
EXPORT void MKtfhe_bootstrapFFT(MKLweSample *result, const MKLweBootstrappingKeyFFT *bkFFT, Torus32 mu, 
        const MKLweSample *x, const LweParams* LWEparams, const LweParams* extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey); 










// Encrypt and decrypt for gate bootstrap
/** encrypts a boolean */
EXPORT void MKbootsSymEncrypt(MKLweSample *result, int32_t message, const MKLweKey* key);
/** decrypts a boolean */
EXPORT int32_t MKbootsSymDecrypt(const MKLweSample *sample, const MKLweKey* key);








// MK Bootstrapped NAND (no FFT)
EXPORT void MKbootsNAND(MKLweSample *result, const MKLweSample *ca, const MKLweSample *cb, 
        const MKLweBootstrappingKey *bk, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams);
EXPORT void MKbootsNAND_FFT(MKLweSample *result, const MKLweSample *ca, const MKLweSample *cb, 
        const MKLweBootstrappingKeyFFT *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey);
















/* **************************************************************************
***************************** VERSION 2 *************************************
************************************************************************** */




/* *************************
******** MK-RGSW ***********
************************* */


/* Uni-Encrypt */
// Encrypt a integer polynomial as (d, F) = (d, f0, f1)
EXPORT void MKTGswUniEncrypt_v2(MKTGswUESample_v2 *result, IntPolynomial *message, int32_t party, double alpha, const MKRLweKey *key);
EXPORT void MKTGswUniEncryptI_v2(MKTGswUESample_v2 *result, int32_t message, int32_t party, double alpha, const MKRLweKey *key);
// TO BE CODED 
EXPORT void MKtGswSymDecrypt_v2(TorusPolynomial *result, const MKTGswUESample_v2 *sample, const MKRLweKey *key);







/* EXPAND */
// (d,F) = (d,f0,f1) -> D_i=(x_0, ..., x_{parties-1}, x_parties + d_i, y_0, ..., d_i+y_i, ..., y_perties, d_i)
EXPORT void MKTGswExpand_v2(MKTGswExpSample_v2 *result, MKTGswUESample_v2 *sample, const MKRLweKey *key, const MKTFHEParams* MKparams);
// TO BE CODED 
EXPORT void MKtGswEXPSymDecrypt_v2(TorusPolynomial *result, MKTGswExpSample_v2 *sample, const MKRLweKey *key);
/* EXPAND */
// (d,F) = (d,f0,f1) -> D_i=(x_0, ..., x_{parties-1}, x_parties + d_i, y_0, ..., d_i+y_i, ..., y_perties, d_i)
// sample UE --> resultFFT expand
EXPORT void MKTGswExpandFFT_v2(MKTGswExpSampleFFT_v2 *resultFFT, MKTGswUESample_v2 *sample, const MKRLweKey *key, 
        const TLweParams* RLWEparams, const MKTFHEParams* MKparams);








#endif //MKTFHEFUNCTIONS_H

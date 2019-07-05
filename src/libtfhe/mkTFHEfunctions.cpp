#include <cstdlib>
#include <iostream>
#include <random>
#include <cassert>
#include "tfhe_generic_templates.h"
#include "tfhe_core.h"
#include "numeric_functions.h"
#include "lweparams.h"
#include "lwekey.h"
#include "lwe-functions.h"
#include "lwesamples.h"
#include "lwekeyswitch.h"
#include "tlwe_functions.h"
#include "polynomials_arithmetic.h"
#include "lagrangehalfc_arithmetic.h"




#include "mkTFHEparams.h"
#include "mkTFHEkeys.h"
#include "mkTFHEsamples.h"
#include "mkTFHEfunctions.h"


using namespace std;


/* Encryption and Decryption for the MK samples */


/* *************************
******** MK-LWE ************
************************* */

// b = \sum <a_i,s_i> + m + e
// m comes with a scaling factor 
EXPORT void MKlweSymEncrypt(MKLweSample* result, Torus32 message, double alpha, const MKLweKey* key){
    const int32_t n = key->LWEparams->n;
    const int32_t parties = key->MKparams->parties;
    
    result->b = gaussian32(message, alpha); 

    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            result->a[i*n +j] = uniformTorus32_distrib(generator);
            result->b += result->a[i*n +j]*key->key[i].key[j];
        } 
    }
    
    result->current_variance = alpha*alpha;
}



/* 
 * This function encrypts a message by using key and a given noise value
*/
EXPORT void MKlweSymEncryptWithExternalNoise(MKLweSample* result, Torus32 message, double noise, double alpha, 
        const MKLweKey* key)
{
    const int32_t n = key->LWEparams->n;
    const int32_t parties = key->MKparams->parties;
    
    result->b = message + dtot32(noise); 

    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            result->a[i*n +j] = uniformTorus32_distrib(generator);
            result->b += result->a[i*n +j]*key->key[i].key[j];
        } 
    }
    
    result->current_variance = alpha*alpha;
}






/**
 * This function computes the phase of sample by using key : phi = b - \sum <a_i,s_i>
 */
EXPORT Torus32 MKlwePhase(const MKLweSample* sample, const MKLweKey* key){
    const int32_t n = key->LWEparams->n;
    const int32_t parties = key->MKparams->parties;
    
    Torus32 axs = 0;

    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            axs += sample->a[i*n +j]*key->key[i].key[j];
        } 
    }

    return sample->b - axs;
}




/**
 * This function computes the decryption of sample by using key
 * The constant Msize indicates the message space and is used to approximate the phase
 */
EXPORT Torus32 MKlweSymDecrypt(const MKLweSample* sample, const MKLweKey* key, const int32_t Msize){
    Torus32 phi;

    phi = MKlwePhase(sample, key);
    return approxPhase(phi, Msize);
}






/** result = (0, ..., 0, mu) */
EXPORT void MKlweNoiselessTrivial(MKLweSample* result, Torus32 mu, const MKTFHEParams* params){
    const int32_t parties = params->parties;
    const int32_t n = params->n;

    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            result->a[i*n +j] = 0;
        } 
    }

    result->b = mu;
    result->current_variance = 0.0;
}




/** result = result - sample */
EXPORT void MKlweSubTo(MKLweSample* result, const MKLweSample* sample, const MKTFHEParams* MKparams){
    const int32_t n = MKparams->n;
    const int32_t parties = MKparams->parties;

    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            result->a[i*n+j] -= sample->a[i*n+j];
        }
    }
    
    result->b -= sample->b;

    result->current_variance += sample->current_variance; 
}


/** result = sample */
EXPORT void MKlweCopy(MKLweSample* result, const MKLweSample* sample, const MKTFHEParams* params){
    const int32_t n = params->n;
    const int32_t parties = params->parties;

    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            result->a[i*n+j] = sample->a[i*n+j];
        }
    }
    
    result->b = sample->b;

    result->current_variance = sample->current_variance;
}






















/* *************************
******** MK-RLWE ***********
************************* */


// b = \sum <a_i,s_i> + m + e
// m comes with a scaling factor 
EXPORT void MKtLweSymEncrypt(MKTLweSample *result, TorusPolynomial *message, double alpha, const MKRLweKey *key) {
    const int32_t N = key->RLWEparams->N;
    const int32_t parties = key->MKparams->parties;

    for (int j = 0; j < N; ++j)
    {
        result->b->coefsT[j] = gaussian32(0, alpha);
        result->b->coefsT[j] += message->coefsT[j];
    }

    for (int i = 0; i < parties; ++i)
    {
        torusPolynomialUniform(&result->a[i]);
        torusPolynomialAddMulR(result->b, key->key[i].key, &result->a[i]);
    }

    result->current_variance = alpha * alpha;
}




// b = \sum <a_i,s_i> + m + e
// m constant message, comes with a scaling factor 
EXPORT void MKtLweSymEncryptT(MKTLweSample *result, Torus32 message, double alpha, const MKRLweKey *key) {
    const int32_t N = key->RLWEparams->N;
    const int32_t parties = key->MKparams->parties;

    for (int j = 0; j < N; ++j)
    {
        result->b->coefsT[j] = gaussian32(0, alpha);
    }
    result->b->coefsT[0] += message;

    for (int i = 0; i < parties; ++i)
    {
        torusPolynomialUniform(&result->a[i]);
        torusPolynomialAddMulR(result->b, key->key[i].key, &result->a[i]);
    }

    result->current_variance = alpha * alpha;
}









/** result = (0, ..., 0,mu) */
EXPORT void MKtLweNoiselessTrivial(MKTLweSample *result, const TorusPolynomial *mu, const MKTFHEParams *MKparams) {
    const int32_t parties = MKparams->parties;

    for (int i = 0; i < parties; ++i)
    {
        torusPolynomialClear(&result->a[i]);
    }

    torusPolynomialCopy(result->b, mu);

    result->current_variance = 0.0;
}













/**
 * This function computes the phase of sample by using key : phi = b - \sum <a_i,s_i>
 */
EXPORT void MKtLwePhase(TorusPolynomial *phase, const MKTLweSample *sample, const MKRLweKey *key) {
    const int32_t parties = key->MKparams->parties;

    torusPolynomialCopy(phase, sample->b); // phi = b

    for (int i = 0; i < parties; ++i)
    {
        torusPolynomialSubMulR(phase, key->key[i].key, &sample->a[i]);
    }
}








/**
 * This function computes the decryption of sample by using key
 * The constant Msize indicates the message space and is used to approximate the phase
 */
EXPORT void MKtLweSymDecrypt(TorusPolynomial *result, const MKTLweSample *sample, const MKRLweKey *key, int32_t Msize) {
    MKtLwePhase(result, sample, key);
    tLweApproxPhase(result, result, Msize, key->RLWEparams->N);
}


/**
 * This function computes the decryption of sample by using key
 * The constant Msize indicates the message space and is used to approximate the phase
 * The message is a constant torus element
 */
EXPORT Torus32 MKtLweSymDecryptT(const MKTLweSample *sample, const MKRLweKey *key, int32_t Msize) {
    TorusPolynomial *phase = new_TorusPolynomial(key->RLWEparams->N);

    MKtLwePhase(phase, sample, key);
    Torus32 result = approxPhase(phase->coefsT[0], Msize);

    delete_TorusPolynomial(phase);
    return result;
}









/** result = sample */
EXPORT void MKtLweCopy(MKTLweSample *result, const MKTLweSample *sample, const MKTFHEParams *MKparams) {
    const int32_t parties = MKparams->parties;
    const int32_t N = MKparams->N;

    for (int32_t i = 0; i <= parties; ++i)
    {
        for (int32_t j = 0; j < N; ++j)
        {
            result->a[i].coefsT[j] = sample->a[i].coefsT[j];
        }
    }

    result->current_variance = sample->current_variance;
}








// external multiplication of ACC by X^ai-1
EXPORT void MKtLweMulByXaiMinusOne(MKTLweSample *result, int32_t ai, const MKTLweSample *ACC, 
        const MKTFHEParams *MKparams) 
{
    const int32_t parties = MKparams->parties;

    for (int i = 0; i <= parties; ++i)
    {
        torusPolynomialMulByXaiMinusOne(&result->a[i], ai, &ACC->a[i]);
    }
}






/** result = result + sample */
EXPORT void MKtLweAddTo(MKTLweSample *result, const MKTLweSample *sample, const MKTFHEParams *MKparams) {
    const int32_t parties = MKparams->parties;

    for (int i = 0; i < parties; ++i)
    {
        torusPolynomialAddTo(&result->a[i], &sample->a[i]);
    }
    torusPolynomialAddTo(result->b, sample->b);

    result->current_variance += sample->current_variance;
}





// EXTRACT
EXPORT void MKtLweExtractMKLweSampleIndex(MKLweSample* result, const MKTLweSample* x, const int32_t index, 
        const MKTFHEParams* MKparams) 
{
    const int32_t parties = MKparams->parties;
    const int32_t N = MKparams->N;

    assert(MKparams->n_extract == N);

    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j <= index; ++j)
        {
            result->a[i*N+j] = x->a[i].coefsT[index-j];
        }
        for (int j = index+1; j < N; ++j)
        {
            result->a[i*N+j] = -x->a[i].coefsT[N+index-j];
        }
    }

    result->b = x->b->coefsT[index];
}
// extract index 0
EXPORT void MKtLweExtractMKLweSample(MKLweSample* result, const MKTLweSample* x, const MKTFHEParams* MKparams) {
    MKtLweExtractMKLweSampleIndex(result, x, 0, MKparams);
}






















































/* *************************
******** MK-RGSW ***********
************************* */


/* Uni-Encrypt */
// Encrypt a integer polynomial
EXPORT void MKTGswUniEncrypt(MKTGswUESample *result, IntPolynomial *message, int32_t party, double alpha, const MKRLweKey *key) {
    const int32_t N = key->RLWEparams->N;
    const int32_t dg = key->MKparams->dg;
    const int32_t parties = key->MKparams->parties;

    // generate r, the randomness
    uniform_int_distribution<int32_t> distribution(0, 1);
    IntPolynomial* r = new_IntPolynomial(N);
    for (int j = 0; j < N; ++j){
        r->coefs[j] = distribution(generator);
    }
            

    // C = (c0,c1) \in T^2dg, with c0 = s_party*c1 +e_c + m*g
    for (int i = 0; i < dg; ++i)
    {
        // c1  = c+dg
        torusPolynomialUniform(&result->c[dg + i]); 

        // c0 = e_c[i] + m*g[i]
        for (int j = 0; j < N; ++j)
        {
            result->c[i].coefsT[j] = gaussian32(0, alpha); // e_c
            result->c[i].coefsT[j] += message->coefs[j] * key->MKparams->g[i]; // m*g[i]
        }

        torusPolynomialAddMulR(&result->c[i], key->key[party].key, &result->c[dg + i]);       
    }


    // D = (d0, d1) = r*[Pkey_party | Pkey_parties] + [E0|E1] + [0|m*g] \in T^2dg
    for (int i = 0; i < dg; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            // d1  = d+dg
            result->d[dg + i].coefsT[j] = gaussian32(0, alpha); // E1
            // d1 = E1 + m*g[i]
            result->d[dg + i].coefsT[j] += message->coefs[j] * key->MKparams->g[i]; // m*g[i]

            // d0 = E0
            result->d[i].coefsT[j] = gaussian32(0, alpha); // E0
        }

        // d1 = r*Pkey_parties[i] + E1 + m*g[i]
        torusPolynomialAddMulR(&result->d[dg + i], r, &key->Pkey[parties*dg + i]); 
        // d0 = r*Pkey_parties[i] + E0
        torusPolynomialAddMulR(&result->d[i], r, &key->Pkey[party*dg + i]);       
    }


    // F = (f0,f1) \in T^2dg, with f0 = s_party*f1 +e_f + r*g
    for (int i = 0; i < dg; ++i)
    {
        // f1  = f+dg
        torusPolynomialUniform(&result->f[dg + i]); 

        // f0 = e_f[i] + r*g[i]
        for (int j = 0; j < N; ++j)
        {
            result->f[i].coefsT[j] = gaussian32(0, alpha); // e_f
            result->f[i].coefsT[j] += r->coefs[j] * key->MKparams->g[i]; // r*g[i]
        }

        torusPolynomialAddMulR(&result->f[i], key->key[party].key, &result->f[dg + i]);       
    }
    
   

    result->current_variance = alpha * alpha;
    delete_IntPolynomial(r);
}







// Encrypt an integer value
EXPORT void MKTGswUniEncryptI(MKTGswUESample *result, int32_t message, int32_t party, double alpha, const MKRLweKey *key) {
    // cout << "alpha for uniencrypt = " << alpha << endl;
    const int32_t N = key->RLWEparams->N;
    const int32_t dg = key->MKparams->dg;
    const int32_t parties = key->MKparams->parties;

    // generate r, the randomness
    uniform_int_distribution<int32_t> distribution(0, 1);
    IntPolynomial* r = new_IntPolynomial(N);

    
    // cout << "r = "; // just for verification of f part
    for (int j = 0; j < N; ++j){
        //r ->coefs[j] = 0;
        r->coefs[j] = distribution(generator);
        // cout << r->coefs[j] << ", "; // just for verification of f part 
    }
    // cout << endl; // just for verification of f part 


    // C = (c0,c1) \in T^2dg, with c0 = s_party*c1 + e_c + m*g
    for (int i = 0; i < dg; ++i)
    {
        // c1  = c+dg
        torusPolynomialUniform(&result->c[dg + i]); 

        // c0 = e_c[i] + m*g[i]
        for (int j = 0; j < N; ++j)
        {
            result->c[i].coefsT[j] = gaussian32(0, alpha); // e_c
        }
        result->c[i].coefsT[0] += message * key->MKparams->g[i]; // m*g[i]
        // c0 = s_party*c1 + e_c + m*g
        torusPolynomialAddMulR(&result->c[i], key->key[party].key, &result->c[dg + i]);       
    }


    // D = (d0, d1) = r*[Pkey_party | Pkey_parties] + [E0|E1] + [0|m*g] \in T^2dg
    for (int i = 0; i < dg; ++i)
    {
        //cout << "E:";
        for (int j = 0; j < N; ++j)
        {
            // d1  = d+dg
            result->d[dg + i].coefsT[j] = gaussian32(0, alpha); // E1
            // d0 = E0
            result->d[i].coefsT[j] = gaussian32(0, alpha); // E0
            // cout << result->d[i].coefsT[j] << ", "; 
        }
        //cout << endl;
        // d1 = E1 + m*g[i]
        result->d[dg + i].coefsT[0] += message * key->MKparams->g[i]; // m*g[i]

        // d1 = r*Pkey_parties[i] + E1 + m*g[i] 
        torusPolynomialAddMulR(&result->d[dg + i], r, &key->Pkey[dg*parties + i]); 
        //cout << "d1 vec = ";
        //for(int j = 0; j < N; j++){
        //    cout << result->d[dg+i].coefsT[j] << ", "; 
        // }
        //cout << endl; 
        // d0 = r*Pkey_party[i] + E0
        torusPolynomialAddMulR(&result->d[i], r, &key->Pkey[dg*party + i]);   
        //cout << "d0 vec = ";
        //for(int j = 0; j < N; j++){
        //    cout << result->d[i].coefsT[j] << ", "; 
        //}
        //cout << endl;     
    }


    // F = (f0,f1) \in T^2dg, with f0 = s_party*f1 + e_f + r*g
    for (int i = 0; i < dg; ++i)
    {
        // f1  = f+dg
        torusPolynomialUniform(&result->f[dg + i]); 

        // f0 = e_f[i] + r*g[i]
        for (int j = 0; j < N; ++j)
        {
            result->f[i].coefsT[j] = gaussian32(0, alpha); // e_f
            result->f[i].coefsT[j] += r->coefs[j] * key->MKparams->g[i]; // r*g[i]
        }
        // f0 = s_party*f1 + e_f + r*g
        torusPolynomialAddMulR(&result->f[i], key->key[party].key, &result->f[dg + i]);       
    }
    
   

    result->current_variance = alpha * alpha;
    delete_IntPolynomial(r);
}







/**
 * This function computes the decryption (actually the phase) of sample by using key
 * The constant Msize indicates the message space and is used to approximate the phase
 */
// result is an array composed by 3*dg torus polynomials 
EXPORT void MKtGswSymDecrypt(TorusPolynomial *result, const MKTGswUESample *sample, const MKRLweKey *key) {
    const int32_t dg = key->MKparams->dg;
    const int32_t party = sample->party;
    

    for (int j = 0; j < dg; ++j)
    {
        // c part
        torusPolynomialCopy(&result[j], &sample->c[j]); // phi = c0[j] 
        // phi = c0 - c1*s_party
        torusPolynomialSubMulR(&result[j], key->key[party].key, &sample->c[dg+j]);

        // d part
        // cout << "d0 in decrypt = ";
        // for (int h = 0; h < key->MKparams->N; h++){
        //     cout << sample->d[j].coefsT[h] << ",";
        // }
        // cout << endl;
        torusPolynomialCopy(&result[dg+j], &sample->d[j]); // phi = d0[j] 
        // phi = - (d0 - d1*s_party)
        torusPolynomialSubMulR2(&result[dg+j], key->key[party].key, &sample->d[dg+j]);
        // cout << "d1 in decrypt = ";
        // for (int h = 0; h < key->MKparams->N; h++){
        //     cout << sample->d[dg+j].coefsT[h] << ",";
        // }
        // cout << endl;

        // f part 
        torusPolynomialCopy(&result[2*dg + j], &sample->f[j]); // phi = f0[j] 
        // phi = f0 - f1*s_party
        torusPolynomialSubMulR(&result[2*dg + j], key->key[party].key, &sample->f[dg+j]);
    }

}


































// same function as tGswTorus32PolynomialDecompH, without the assembly
// (t_0, ..., t_N-1) -> (I_0, ...,I_dg-1)
// decomp_g(t_j) = (I0,j, ..., Idg-1,j)
EXPORT void MKtGswTorus32PolynomialDecompG(IntPolynomial *result, const TorusPolynomial *sample, 
        const MKTFHEParams *params) 
{
    const int32_t N = params->N;
    const int32_t dg = params->dg;
    const int32_t Bgbit = params->Bgbit;

    uint32_t *buf = (uint32_t *) sample->coefsT;

    const uint32_t maskMod = params->maskMod; // Bg - 1
    const int32_t halfBg = params->halfBg; // Bg / 2
    const uint32_t offset = params->offset;


    //First, add offset to everyone
    for (int32_t j = 0; j < N; ++j) buf[j] += offset;

    //then, do the decomposition 
    for (int32_t p = 0; p < dg; ++p) {
        const int32_t decal = (32 - (p + 1) * Bgbit);
        int32_t *res_p = result[p].coefs;

        for (int32_t j = 0; j < N; ++j) {
            uint32_t temp1 = (buf[j] >> decal) & maskMod;
            res_p[j] = temp1 - halfBg;
        }
    }

    //finally, remove offset to everyone
    for (int32_t j = 0; j < N; ++j) buf[j] -= offset;

}




EXPORT void MKtGswTorus32PolynomialDecompGassembly(IntPolynomial *result, const TorusPolynomial *sample, 
        const MKTFHEParams *params)
{
    const int32_t N = params->N;
    const int32_t dg = params->dg;
    const int32_t Bgbit = params->Bgbit;
    uint32_t *buf = (uint32_t *) sample->coefsT;

//#define __AVX2__ //(to test)
#ifndef __AVX2__
    const uint32_t maskMod = params->maskMod;
    const int32_t halfBg = params->halfBg;
    const uint32_t offset = params->offset;
#else
    const uint32_t* maskMod_addr = &params->maskMod;
    const int32_t* halfBg_addr = &params->halfBg;
    const uint32_t* offset_addr = &params->offset;
#endif





    //First, add offset to everyone
#ifndef __AVX2__
    for (int32_t j = 0; j < N; ++j) buf[j] += offset;
#else
    {
    const uint32_t* sit = buf;
    const uint32_t* send = buf+N;
    __asm__ __volatile__ (
        "vpbroadcastd (%2),%%ymm0\n"
        "0:\n"
        "vmovdqu (%0),%%ymm3\n"
        "vpaddd %%ymm0,%%ymm3,%%ymm3\n" // add offset
        "vmovdqu %%ymm3,(%0)\n"
        "addq $32,%0\n"
        "cmpq %1,%0\n"
        "jb 0b\n"
        : "=r"(sit),"=r"(send),"=r"(offset_addr)
        :  "0"(sit), "1"(send), "2"(offset_addr)
        : "%ymm0","%ymm3","memory"
        );
    }
#endif





    //then, do the decomposition (in parallel)
    for (int32_t p = 0; p < dg; ++p) {
        const int32_t decal = (32 - (p + 1) * Bgbit);

#ifndef __AVX2__
        int32_t *res_p = result[p].coefs;
        for (int32_t j = 0; j < N; ++j) {
            uint32_t temp1 = (buf[j] >> decal) & maskMod;
            res_p[j] = temp1 - halfBg;
        }
#else
        int32_t* dst = result[p].coefs;
        const uint32_t* sit = buf;
        const uint32_t* send = buf+N;
        const int32_t* decal_addr = &decal;
        __asm__ __volatile__ (
            "vpbroadcastd (%4),%%ymm0\n"
            "vpbroadcastd (%5),%%ymm1\n"
            "vmovd (%3),%%xmm2\n"
            "1:\n"
            "vmovdqu (%1),%%ymm3\n"
            "VPSRLD %%xmm2,%%ymm3,%%ymm3\n" // shift by decal
            "VPAND %%ymm1,%%ymm3,%%ymm3\n"  // and maskMod
            "VPSUBD %%ymm0,%%ymm3,%%ymm3\n" // sub halfBg
            "vmovdqu %%ymm3,(%0)\n"
            "addq $32,%0\n"
            "addq $32,%1\n"
            "cmpq %2,%1\n"
            "jb 1b\n"
            : "=r"(dst),"=r"(sit),"=r"(send),"=r"(decal_addr),"=r"(halfBg_addr),"=r"(maskMod_addr)
            :  "0"(dst), "1"(sit), "2"(send), "3"(decal_addr), "4"(halfBg_addr) ,"5"(maskMod_addr)
            : "%ymm0","%ymm1","%ymm2","%ymm3","memory"
            );
        /* // verify that the assembly block was ok
        int32_t* res_p = result[p].coefs;
        for (int32_t j = 0; j < N; ++j)
        {
            uint32_t temp1 = (buf[j] >> decal) & maskMod;
            if (res_p[j] != int32_t(temp1 - halfBg)) {
            fprintf(stderr, "j=%d,buf[j]=%u,decal=%u,mask=%u,halfbg=%d,res_p[j]=%d\n",j,buf[j],decal,maskMod,halfBg,res_p[j]);
            abort();
            }
        }*/
#endif
    }






    //finally, remove offset to everyone
#ifndef __AVX2__
    for (int32_t j = 0; j < N; ++j) buf[j] -= offset;
#else
    {
    const uint32_t* sit = buf;
    const uint32_t* send = buf+N;
    __asm__ __volatile__ (
        "vpbroadcastd (%2),%%ymm0\n"
        "2:\n"
        "vmovdqu (%0),%%ymm3\n"
        "vpsubd %%ymm0,%%ymm3,%%ymm3\n" // add offset
        "vmovdqu %%ymm3,(%0)\n"
        "addq $32,%0\n"
        "cmpq %1,%0\n"
        "jb 2b\n"
        "vzeroall\n"
        : "=r"(sit),"=r"(send),"=r"(offset_addr)
        :  "0"(sit), "1"(send), "2"(offset_addr)
        : "%ymm0","%ymm3","memory"
        );
    }
#endif
}











/* EXPAND */
// (C,D,F) = (c0,c1,d0,d1,f0,f1) -> C = (y1, ..., d1, ..., yk, c1, d0+x1, ..., d0, ..., d0+xk, c0)
EXPORT void MKTGswExpand(MKTGswExpSample *result, MKTGswUESample *sample, const MKRLweKey *key, const MKTFHEParams* MKparams) {
    const int32_t N = key->RLWEparams->N;
    const int32_t dg = key->MKparams->dg;
    const int32_t party = sample->party;
    const int32_t parties = key->MKparams->parties;


    // INITIALIZE
    // C = (0, ..., d1, ..., 0, c1, d0, ..., d0, ..., d0, c0)
    for (int j = 0; j < dg; ++j)
    {
        // Initialize y: (0)*dg ... (parties-1)*dg as 0
        // d1 will stay only in place (party)*dg
        for (int i = 0; i <= parties; ++i)
        {
            torusPolynomialClearN(&result->y[i*dg + j], N);
        }

        // d1 in place y[party*dg]
        torusPolynomialCopyN(&result->y[party*dg + j], &sample->d[dg + j], N);

        // c1 in place c1
        torusPolynomialCopyN(&result->c1[j], &sample->c[dg + j], N);
        
        // Initialize x: 0*dg ... (parties-1)*dg as d0
        // d0 will stay only in place (party)*dg
        for (int i = 0; i < parties; ++i)
        {
            torusPolynomialCopyN(&result->x[i*dg + j], &sample->d[j], N);
        }

        // c0 in place c0
        torusPolynomialCopyN(&result->c0[j], &sample->c[j], N);
    }



    TorusPolynomial* X = new_TorusPolynomial(N);
    TorusPolynomial* Y = new_TorusPolynomial(N);
    TorusPolynomial* b_temp = new_TorusPolynomial(N);
    IntPolynomial* u = new_IntPolynomial_array(dg, N);

    // 0 ... party-1
    for (int i = 0; i < party; ++i) 
    {
        for (int j = 0; j < dg; ++j)
        {
            // b_temp = b_i[j] - b_party[j] = Pkey[i*dg + j] - Pkey[party*dg + j]
            torusPolynomialSubN(b_temp, &key->Pkey[i*dg + j], &key->Pkey[party*dg + j], N);
            // g^{-1}(b_temp) = [u_0, ...,u_dg-1] intPolynomials
            MKtGswTorus32PolynomialDecompGassembly(u, b_temp, MKparams);

            // X=0 and Y=0
            torusPolynomialClearN(X, N);
            torusPolynomialClearN(Y, N);
            for (int l = 0; l < dg; ++l)
            {
                // X = xi[j] = <g^{-1}(b_temp), f0>
                torusPolynomialAddMulRFFTN(X, &u[l], &sample->f[l], N);
                // Y = yi[j] = <g^{-1}(b_temp), f1> 
                torusPolynomialAddMulRFFTN(Y, &u[l], &sample->f[dg + l], N);          
            }
            
            // xi = d0 + xi
            torusPolynomialAddTo1(&result->x[i*dg + j], X); // N = X->N
            // yi = 0 + yi
            torusPolynomialAddTo1(&result->y[i*dg + j], Y); // N = Y->N
        }   
    }
    // party+1 ... parties-1
    for (int i = party+1; i < parties; ++i) 
    {
        for (int j = 0; j < dg; ++j)
        {
            // b_temp = b_i[j] - b_party[j] = Pkey[i*dg + j] - Pkey[party*dg + j]
            torusPolynomialSubN(b_temp, &key->Pkey[i*dg + j], &key->Pkey[party*dg + j], N);
            // g^{-1}(b_temp) = [u_0, ...,u_dg-1] intPolynomials
            MKtGswTorus32PolynomialDecompGassembly(u, b_temp, MKparams);

            // X=0 and Y=0
            torusPolynomialClearN(X, N);
            torusPolynomialClearN(Y, N);
            for (int l = 0; l < dg; ++l)
            {
                // X = xi[j] = <g^{-1}(b_temp), f0>
                torusPolynomialAddMulRFFTN(X, &u[l], &sample->f[l], N);
                // Y = yi[j] = <g^{-1}(b_temp), f1> 
                torusPolynomialAddMulRFFTN(Y, &u[l], &sample->f[dg + l], N);          
            }
            
            // xi = d0 + xi
            torusPolynomialAddTo1(&result->x[i*dg + j], X); // N = X->N
            // yi = 0 + yi
            torusPolynomialAddTo1(&result->y[i*dg + j], Y); // N = Y->N
        }    
    }


    result->party = sample->party;
    // TODO: fix this
    result->current_variance = sample->current_variance;
    
    delete_TorusPolynomial(X);
    delete_TorusPolynomial(Y);
    delete_TorusPolynomial(b_temp);
    delete_IntPolynomial_array(dg, u);
}







/**
 * This function is used to verify that the expansion is done correctly
 * C = (y1, ..., d1, ..., yk, c1, d0+x1, ..., d0, ..., d0+xk, c0) (k parties)
 * The constant Msize indicates the message space and is used to approximate the phase
 */
// result is an array composed by (parties+1)*dg torus polynomials, containing the phases 
EXPORT void MKtGswEXPSymDecrypt(TorusPolynomial *result, MKTGswExpSample *sample, const MKRLweKey *key) {
    const int32_t parties = key->MKparams->parties;
    const int32_t dg = key->MKparams->dg;
    const int32_t party = sample->party;


    //cout << "party id for GSW decryption :" << party << endl;


    // for i= 0,...,parties-1, i!=party
    // phi_i = -(d0 + xi - d1*s_i - yi*s_party) = E_i + mu*g*s_i
    for (int i = 0; i < party; ++i)
    {
        for (int j = 0; j < dg; ++j)
        {
            // phi_i = d0 + xi 
            torusPolynomialCopy(&result[i*dg + j], &sample->x[i*dg + j]);
            // phi_i = -(d0 + xi) + d1*s_i
            torusPolynomialSubMulR2(&result[i*dg + j], key->key[i].key, &sample->y[party*dg + j]);
            // phi_i = -(d0 + xi) + d1*s_i + yi*s_party
            torusPolynomialAddMulR(&result[i*dg + j], key->key[party].key, &sample->y[i*dg + j]);
        }
    }
    for (int i = party + 1; i < parties; ++i)
    {
        for (int j = 0; j < dg; ++j)
        {
            // phi_i = d0 + xi 
            torusPolynomialCopy(&result[i*dg + j], &sample->x[i*dg + j]);
            // phi_i = -(d0 + xi) + d1*s_i
            torusPolynomialSubMulR2(&result[i*dg + j], key->key[i].key, &sample->y[party*dg + j]);
            // phi_i = -(d0 + xi) + d1*s_i + yi*s_party
            torusPolynomialAddMulR(&result[i*dg + j], key->key[party].key, &sample->y[i*dg + j]);
        }
    }

    // phi_party = -(d0 - d1*s_party) = E_party + mu*g*s_party
    for (int j = 0; j < dg; ++j)
    {
        // phi_party = d0
        torusPolynomialCopy(&result[party*dg + j], &sample->x[party*dg + j]);
        // phi_party = -(d0 - d1*s_party)
        torusPolynomialSubMulR2(&result[party*dg + j], key->key[party].key, &sample->y[party*dg + j]); 
    }



    // phi_parties = phi_b = c0 - c1*s_party = E_b + mu*g
    for (int j = 0; j < dg; ++j)
    {
        // phi_b = c0
        torusPolynomialCopy(&result[parties*dg + j], &sample->c0[j]);
        // phi_b = c0 - c1*s_party
        torusPolynomialSubMulR(&result[parties*dg + j], key->key[party].key, &sample->c1[j]); 
    }


}




















/* EXPAND */
// (C,D,F) = (c0,c1,d0,d1,f0,f1) -> C = (y1, ..., d1, ..., yk, c1, d0+x1, ..., d0, ..., d0+xk, c0)
// sample UE --> resultFFT expand
EXPORT void MKTGswExpandFFT(MKTGswExpSampleFFT *resultFFT, MKTGswUESample *sample, const MKRLweKey *key, 
        const TLweParams* RLWEparams, const MKTFHEParams* MKparams) 
{
    const int32_t N = key->RLWEparams->N;
    const int32_t dg = key->MKparams->dg;
    const int32_t party = sample->party;
    const int32_t parties = key->MKparams->parties;

    LagrangeHalfCPolynomial* tempFFT = new_LagrangeHalfCPolynomial(N);
    

    // sample -> sampleFFT
    LagrangeHalfCPolynomial *arr = new_LagrangeHalfCPolynomial_array(6*dg, N);
    for (int i = 0; i < 6*dg; ++i)
    {
        TorusPolynomial_ifft(&arr[i], &sample->c[i]);
    }
    MKTGswUESampleFFT* sampleFFT = new_MKTGswUESampleFFT(RLWEparams, MKparams, arr, sample->current_variance);



    // INITIALIZE
    // C = (0, ..., d1, ..., 0, c1, d0, ..., d0, ..., d0, c0)
    for (int j = 0; j < dg; ++j)
    {
        // Initialize y: (0)*dg ... (parties-1)*dg as 0
        // d1 will stay only in place (party)*dg
        for (int i = 0; i < parties; ++i)
        {
            LagrangeHalfCPolynomialClear(&resultFFT->y[i*dg + j]);
        }

        // d1 in place y[party*dg]
        LagrangeHalfCPolynomialCopy(&resultFFT->y[party*dg + j], &sampleFFT->d[dg + j]);
        
        // c1 in place c1
        LagrangeHalfCPolynomialCopy(&resultFFT->c1[j], &sampleFFT->c[dg + j]);
        
        // Initialize x: 0*dg ... (parties-1)*dg as d0
        // d0 will stay only in place (party)*dg
        for (int i = 0; i < parties; ++i)
        {
            LagrangeHalfCPolynomialCopy(&resultFFT->x[i*dg + j], &sampleFFT->d[j]);
        }

        // c0 in place c0
        LagrangeHalfCPolynomialCopy(&resultFFT->c0[j], &sampleFFT->c[j]);
    }





    LagrangeHalfCPolynomial* X = new_LagrangeHalfCPolynomial(N);
    LagrangeHalfCPolynomial* Y = new_LagrangeHalfCPolynomial(N);
    TorusPolynomial* b_temp = new_TorusPolynomial(N);
    IntPolynomial* u = new_IntPolynomial_array(dg, N);
    LagrangeHalfCPolynomial *uFFT = new_LagrangeHalfCPolynomial_array(dg, N); //fft version

    // 0 ... party-1
    for (int i = 0; i < party; ++i) 
    {
        for (int j = 0; j < dg; ++j)
        {
            // b_temp = b_i[j] - b_party[j] = Pkey[i*dg + j] - Pkey[party*dg + j]
            torusPolynomialSubN(b_temp, &key->Pkey[i*dg + j], &key->Pkey[party*dg + j], N);
            // g^{-1}(b_temp) = [u_0, ...,u_dg-1] intPolynomials
            MKtGswTorus32PolynomialDecompGassembly(u, b_temp, MKparams);
            for (int p = 0; p < dg; ++p){
                IntPolynomial_ifft(&uFFT[p], &u[p]); // FFT
            }

            // X=0 and Y=0
            LagrangeHalfCPolynomialClear(X);
            LagrangeHalfCPolynomialClear(Y);
            for (int l = 0; l < dg; ++l)
            {
                // X = xi[j] = <g^{-1}(b_temp), f0>
                LagrangeHalfCPolynomialMul(tempFFT, &uFFT[l], &sampleFFT->f[l]);
                LagrangeHalfCPolynomialAddTo(X, tempFFT);
                // Y = yi[j] = <g^{-1}(b_temp), f1> 
                LagrangeHalfCPolynomialMul(tempFFT, &uFFT[l], &sampleFFT->f[dg + l]);
                LagrangeHalfCPolynomialAddTo(Y, tempFFT);         
            }
            
            // xi = d0 + xi
            LagrangeHalfCPolynomialAddTo(&resultFFT->x[i*dg + j], X);
            // yi = 0 + yi
            LagrangeHalfCPolynomialAddTo(&resultFFT->y[i*dg + j], Y);
        }   
    }

    // party+1 ... parties-1
    for (int i = party+1; i < parties; ++i) 
    {
        for (int j = 0; j < dg; ++j)
        {
            // b_temp = b_i[j] - b_party[j] = Pkey[i*dg + j] - Pkey[party*dg + j]
            torusPolynomialSubN(b_temp, &key->Pkey[i*dg + j], &key->Pkey[party*dg + j], N);
            // g^{-1}(b_temp) = [u_0, ...,u_dg-1] intPolynomials
            MKtGswTorus32PolynomialDecompGassembly(u, b_temp, MKparams);
            for (int p = 0; p < dg; ++p){
                IntPolynomial_ifft(&uFFT[p], &u[p]); // FFT
            }

            // X=0 and Y=0
            LagrangeHalfCPolynomialClear(X);
            LagrangeHalfCPolynomialClear(Y);
            for (int l = 0; l < dg; ++l)
            {
                // X = xi[j] = <g^{-1}(b_temp), f0>
                LagrangeHalfCPolynomialMul(tempFFT, &uFFT[l], &sampleFFT->f[l]);
                LagrangeHalfCPolynomialAddTo(X, tempFFT);
                // Y = yi[j] = <g^{-1}(b_temp), f1> 
                LagrangeHalfCPolynomialMul(tempFFT, &uFFT[l], &sampleFFT->f[dg + l]);
                LagrangeHalfCPolynomialAddTo(Y, tempFFT);         
            }
            
            // xi = d0 + xi
            LagrangeHalfCPolynomialAddTo(&resultFFT->x[i*dg + j], X);
            // yi = 0 + yi
            LagrangeHalfCPolynomialAddTo(&resultFFT->y[i*dg + j], Y);
        }   
    }





    resultFFT->party = sample->party;
    // TODO: fix this
    resultFFT->current_variance = sample->current_variance;
    


    // delete 
    delete_LagrangeHalfCPolynomial_array(dg, uFFT);
    delete_IntPolynomial_array(dg, u);
    delete_TorusPolynomial(b_temp);
    delete_LagrangeHalfCPolynomial(Y);
    delete_LagrangeHalfCPolynomial(X);
    delete_MKTGswUESampleFFT(sampleFFT);
    delete_LagrangeHalfCPolynomial_array(6*dg, arr);
    delete_LagrangeHalfCPolynomial(tempFFT);
   
}

























/* ********************************************************************************
*************************** EXTERNAL PRODUCT **************************************
******************************************************************************** */


// c' = G^{-1}(c)*C, with C = (y1, ..., d1, ..., yk, c1, d0+x1, ..., d0, ..., d0+xk, c0) 
EXPORT void MKtGswExpExternMulToMKTLwe(MKTLweSample* result, MKTLweSample* sample, const MKTGswExpSample* sampleExp, 
        const MKTFHEParams* MKparams)
{
    const int32_t N = MKparams->N;
    const int32_t dg = MKparams->dg;
    const int32_t party = sampleExp->party;
    const int32_t parties = MKparams->parties;

    // DECOMPOSE sample
    // u[i] = g^{-1}(a_i), u[parties] = g^{-1}(b), 
    IntPolynomial* u = new_IntPolynomial_array((parties+1)*dg, N);
    for (int i = 0; i <= parties; ++i)
    {
        MKtGswTorus32PolynomialDecompGassembly(&u[i*dg], &sample->a[i], MKparams);
    }


    // c'_i = g^{-1}(a_i)*d1, i<parties, i!=party
    for (int i = 0; i < party; ++i)
    {
        torusPolynomialClearN(&result->a[i], N);
        for (int l = 0; l < dg; ++l)
        {
            torusPolynomialAddMulRFFTN(&result->a[i], &u[i*dg+l], &sampleExp->y[party*dg + l], N);
        }
    }
    for (int i = party+1; i < parties; ++i)
    {
        torusPolynomialClearN(&result->a[i], N);
        for (int l = 0; l < dg; ++l)
        {
            torusPolynomialAddMulRFFTN(&result->a[i], &u[i*dg+l], &sampleExp->y[party*dg + l], N);
        }
    }



    // c'_party = \sum g^{-1}(a_i)*yi + g^{-1}(b)*c1
    torusPolynomialClearN(&result->a[party], N);
    for (int i = 0; i < parties; ++i)
    {
        for (int l = 0; l < dg; ++l)
        {
            torusPolynomialAddMulRFFTN(&result->a[party], &u[i*dg+l], &sampleExp->y[dg*i + l], N);
        }
    }
    for (int l = 0; l < dg; ++l)
    {
        torusPolynomialAddMulRFFTN(&result->a[party], &u[parties*dg+l], &sampleExp->c1[l], N);
    }


    // c'_parties = \sum g^{-1}(a_i)*xi + g^{-1}(b)*c0
    torusPolynomialClearN(&result->a[parties], N);
    for (int i = 0; i < parties; ++i)
    {
        for (int l = 0; l < dg; ++l)
        {
            torusPolynomialAddMulRFFTN(&result->a[parties], &u[i*dg+l], &sampleExp->x[dg*i + l], N);
        }
    }
    for (int l = 0; l < dg; ++l)
    {
        torusPolynomialAddMulRFFTN(&result->a[parties], &u[parties*dg+l], &sampleExp->c0[l], N);
    }


    // TODO current_variance

    delete_IntPolynomial_array((parties+1)*dg, u);
}









// External product FFT

// result += poly1*poly2
EXPORT void MulFFTAndAddTo(TorusPolynomial* result, const LagrangeHalfCPolynomial* poly1, 
        const LagrangeHalfCPolynomial* poly2, const int32_t N)
{
    LagrangeHalfCPolynomial* tempFFT = new_LagrangeHalfCPolynomial(N);
    TorusPolynomial* temp = new_TorusPolynomial(N);

    LagrangeHalfCPolynomialMul(tempFFT, poly1, poly2);
    TorusPolynomial_fft(temp, tempFFT);
    torusPolynomialAddTo(result, temp);

    delete_TorusPolynomial(temp);
    delete_LagrangeHalfCPolynomial(tempFFT);
} 

// result -= poly1*poly2
EXPORT void MulFFTAndSubTo(TorusPolynomial* result, const LagrangeHalfCPolynomial* poly1, 
        const LagrangeHalfCPolynomial* poly2, const int32_t N)
{
    LagrangeHalfCPolynomial* tempFFT = new_LagrangeHalfCPolynomial(N);
    TorusPolynomial* temp = new_TorusPolynomial(N);

    LagrangeHalfCPolynomialMul(tempFFT, poly1, poly2);
    TorusPolynomial_fft(temp, tempFFT);
    torusPolynomialSubTo(result, temp);

    delete_TorusPolynomial(temp);
    delete_LagrangeHalfCPolynomial(tempFFT);
} 

// Result not FFT
EXPORT void MKtGswExpExternMulToMKTLweFFT(MKTLweSample* result, MKTLweSample* sample, 
        const MKTGswExpSampleFFT* sampleExpFFT, const TLweParams* RLWEparams, const MKTFHEParams* MKparams)
{
    const int32_t N = MKparams->N;
    const int32_t dg = MKparams->dg;
    const int32_t party = sampleExpFFT->party;
    const int32_t parties = MKparams->parties;
    const int32_t parties1dg = (parties+1)*dg;

    // DECOMPOSE sample and convert it to FFT
    // u[i*dg] = g^{-1}(a_i), u[parties*dg] = g^{-1}(b), 
    IntPolynomial* u = new_IntPolynomial_array(parties1dg, N);
    LagrangeHalfCPolynomial *uFFT = new_LagrangeHalfCPolynomial_array(parties1dg, N); //fft version

    for (int i = 0; i <= parties; ++i){
        MKtGswTorus32PolynomialDecompGassembly(&u[i*dg], &sample->a[i], MKparams);
    }
    for (int p = 0; p < parties1dg; ++p){
        IntPolynomial_ifft(&uFFT[p], &u[p]); // FFT
    }







    
    // c'_i = g^{-1}(a_i)*d1, i<parties, i!=party
    for (int i = 0; i < party; ++i)
    {
        torusPolynomialClearN(&result->a[i], N);
        for (int l = 0; l < dg; ++l)
        {
            MulFFTAndAddTo(&result->a[i], &uFFT[i*dg+l], &sampleExpFFT->y[party*dg + l], N);
        }
    }
    for (int i = party+1; i < parties; ++i)
    {
        torusPolynomialClearN(&result->a[i], N);
        for (int l = 0; l < dg; ++l)
        {
            MulFFTAndAddTo(&result->a[i], &uFFT[i*dg+l], &sampleExpFFT->y[party*dg + l], N);
        }
    }




    // c'_party = \sum g^{-1}(a_i)*yi + g^{-1}(b)*c1
    torusPolynomialClearN(&result->a[party], N);
    for (int i = 0; i < parties; ++i)
    {
        for (int l = 0; l < dg; ++l)
        {
            MulFFTAndAddTo(&result->a[party], &uFFT[i*dg+l], &sampleExpFFT->y[dg*i + l], N);
        }
    }
    for (int l = 0; l < dg; ++l)
    {
        MulFFTAndAddTo(&result->a[party], &uFFT[parties*dg+l], &sampleExpFFT->c1[l], N);
    }




    // c'_parties = \sum g^{-1}(a_i)*xi + g^{-1}(b)*c0
    torusPolynomialClearN(&result->a[parties], N);
    for (int i = 0; i < parties; ++i)
    {
        for (int l = 0; l < dg; ++l)
        {
            MulFFTAndAddTo(&result->a[parties], &uFFT[i*dg+l], &sampleExpFFT->x[dg*i + l], N);
        }
    }
    for (int l = 0; l < dg; ++l)
    {
        MulFFTAndAddTo(&result->a[parties], &uFFT[parties*dg+l], &sampleExpFFT->c0[l], N);
    }



    // TODO current_variance


    // delete 
    delete_LagrangeHalfCPolynomial_array(parties1dg, uFFT);
    delete_IntPolynomial_array(parties1dg, u);
}

























/* ********************************************************************************
****************************** KEY SWITCHING **************************************
******************************************************************************** */

EXPORT void MKlweKeySwitch(MKLweSample* result, const LweKeySwitchKey* ks, const MKLweSample* sample, 
        const LweParams* LWEparams, const MKTFHEParams* MKparams)
{
    const int32_t n_extract = MKparams->n_extract;
    const int32_t Bksbit = MKparams->Bksbit;
    const int32_t dks = MKparams->dks;
    const int32_t parties = MKparams->parties;
    const int32_t Bks = 1 << Bksbit;
    const int32_t prec_offset = 1 << (32-(1+Bksbit*dks)); //precision
    const int32_t mask = Bks-1;
    const int32_t n = LWEparams->n;

    LweSample* temp = new_LweSample(LWEparams);
 
    // result = (b, 0,...,0)
    MKlweNoiselessTrivial(result, sample->b, MKparams);

    for (int p = 0; p < parties; ++p)
    {
        // temp = (0,0)
        lweClear(temp, LWEparams);

        // temp = (a', b')
        for (int i = 0; i < n_extract; ++i)
        {
            const uint32_t aibar = sample->a[p*n_extract + i] + prec_offset;

            for (int j = 0; j < dks; ++j)
            {
                const uint32_t aij = (aibar >> (32-(j+1)*Bksbit)) & mask;
                if(aij != 0) 
                {
                    lweSubTo(temp, &ks[p].ks[i][j][aij], LWEparams);
                }
            }
        }

        // result = (b + \sum b', a1', ..., ak')
        result->b += temp->b;
        for (int i = 0; i < n; ++i)
        {
            result->a[p*n +i] = temp->a[i]; 
        }        
    }

    /*
    for (int p = 0; p < parties; ++p)
    {
        for (int i = 0; i < n_in; ++i)
        {
            const uint32_t aibar = sample->a[p*n_in + i] + prec_offset;

            for (int j = 0; j < dks; ++j)
            {
                const uint32_t aij = (aibar >> (32-(j+1)*Bksbit)) & mask;
                if(aij != 0) 
                {
                    MKlweSubTo(result, &ks->ks[p][i][j][aij], MKparams);
                }
            }
        }
    }
    */

}





































/* ********************************************************************************
****************************** BOOTSTRAPPING **************************************
******************************************************************************** */

// MUX -> rotate
void MKtfhe_MuxRotate(MKTLweSample *result, MKTLweSample *accum, const MKTGswExpSample *bki, 
    const int32_t barai, const TLweParams* RLWEparams, const MKTFHEParams* MKparams) 
{
    MKTLweSample *temp_result = new_MKTLweSample(RLWEparams, MKparams);

    // ACC = BKi*[(X^barai-1)*ACC]+ACC
    // temp = (X^barai-1)*ACC
    MKtLweMulByXaiMinusOne(temp_result, barai, accum, MKparams);
    // temp *= BKi
    MKtGswExpExternMulToMKTLwe(result, temp_result, bki, MKparams);
    // ACC += temp
    MKtLweAddTo(result, accum, MKparams);

    delete_MKTLweSample(temp_result);
}


// MUX -> rotate
void MKtfhe_MuxRotateFFT(MKTLweSample *result, MKTLweSample *accum, const MKTGswExpSampleFFT *bkiFFT, 
    const int32_t barai, const TLweParams* RLWEparams, const MKTFHEParams* MKparams) 
{
    MKTLweSample *temp_result = new_MKTLweSample(RLWEparams, MKparams);

    // ACC = BKi*[(X^barai-1)*ACC]+ACC
    // temp = (X^barai-1)*ACC
    MKtLweMulByXaiMinusOne(temp_result, barai, accum, MKparams);
    // temp *= BKi
    MKtGswExpExternMulToMKTLweFFT(result, temp_result, bkiFFT, RLWEparams, MKparams);
    // ACC += temp
    MKtLweAddTo(result, accum, MKparams);

    delete_MKTLweSample(temp_result);
}











// Blind rotate
EXPORT void MKtfhe_blindRotate(MKTLweSample *accum, const MKTGswExpSample *bk, const int32_t *bara, 
    const TLweParams* RLWEparams, const MKTFHEParams *MKparams) 
{
    const int32_t parties = MKparams->parties;
    const int32_t n = MKparams->n;

    MKTLweSample *temp = new_MKTLweSample(RLWEparams, MKparams);

    // MKTLweSample *temp1 = accum;
    MKTLweSample *temp1 = new_MKTLweSample(RLWEparams, MKparams);
    MKtLweCopy(temp1, accum, MKparams);    


    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            const int32_t baraij = bara[n*i+j];

            if (baraij == 0) continue; //indeed, this is an easy case!

            MKtfhe_MuxRotate(temp, temp1, bk + (n*i+j), baraij, RLWEparams, MKparams);
            swap(temp, temp1);
        }
    }

    if (temp1 != accum) {
        MKtLweCopy(accum, temp1, MKparams);
    }


    delete_MKTLweSample(temp1);
    delete_MKTLweSample(temp);
}





// Blind rotate
//  const MKRLweKey *MKrlwekey ONLY TO TEST
EXPORT void MKtfhe_blindRotateFFT(MKTLweSample *accum, const MKTGswExpSampleFFT *bkFFT, const int32_t *bara, 
    const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey) 
{
    const int32_t parties = MKparams->parties;
    const int32_t n = MKparams->n;

    MKTLweSample *temp = new_MKTLweSample(RLWEparams, MKparams);

    // MKTLweSample *temp1 = accum;
    MKTLweSample *temp1 = new_MKTLweSample(RLWEparams, MKparams);
    MKtLweCopy(temp1, accum, MKparams); 

/*
    const int32_t N = MKparams->N;


    cout << "accum: [";
    TorusPolynomial *phase = new_TorusPolynomial(N);
    MKtLwePhase(phase, temp1, MKrlwekey);

    for (int i = 0; i < N-1; ++i)
    {
        cout << (phase->coefsT[i] > 0 ? 1 : 0) << ", " ;
    }
    cout << (phase->coefsT[N-1] > 0 ? 1 : 0) << "]" << endl;  
*/





    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            const int32_t baraij = bara[n*i+j];

            if (baraij == 0) continue; //indeed, this is an easy case!

            MKtfhe_MuxRotateFFT(temp, temp1, bkFFT + (n*i+j), baraij, RLWEparams, MKparams);
            swap(temp, temp1);

/*
            cout << bkFFT[n*i+j].party;
            cout << " accum: party " << i << " key " << j << "  [";
            MKtLwePhase(phase, temp, MKrlwekey);

            for (int i = 0; i < N-1; ++i)
            {
                cout << (phase->coefsT[i] > 0 ? 1 : 0) << ", " ;
            }
            cout << (phase->coefsT[N-1] > 0 ? 1 : 0) << "]" << endl;  
*/



        }
    }

    if (temp1 != accum) {
        MKtLweCopy(accum, temp1, MKparams);
    }


    delete_MKTLweSample(temp1);
    delete_MKTLweSample(temp);
}






















EXPORT void MKtfhe_blindRotateAndExtract(MKLweSample *result,
                                       const TorusPolynomial *v,
                                       const MKTGswExpSample *bk,
                                       const int32_t barb,
                                       const int32_t *bara,
                                       const TLweParams* RLWEparams, 
                                       const MKTFHEParams *MKparams) 
{
    const int32_t N = MKparams->N;
    const int32_t _2N = 2 * N;

    TorusPolynomial *testvectbis = new_TorusPolynomial(N);
    MKTLweSample *acc = new_MKTLweSample(RLWEparams, MKparams);

    if (barb !=0)
    {
        torusPolynomialMulByXai(testvectbis, _2N - barb, v);
    }
    else
    {
        torusPolynomialCopy(testvectbis, v);
    }

    MKtLweNoiselessTrivial(acc, testvectbis, MKparams);
    MKtfhe_blindRotate(acc, bk, bara, RLWEparams, MKparams); 
    MKtLweExtractMKLweSample(result, acc, MKparams);


    delete_MKTLweSample(acc);
    delete_TorusPolynomial(testvectbis);
}



EXPORT void MKtfhe_blindRotateAndExtractFFT(MKLweSample *result,
                                       const TorusPolynomial *v,
                                       const MKTGswExpSampleFFT *bkFFT,
                                       const int32_t barb,
                                       const int32_t *bara,
                                       const TLweParams* RLWEparams, 
                                       const MKTFHEParams *MKparams, 
                                       const MKRLweKey *MKrlwekey) 
{
    const int32_t N = MKparams->N;
    const int32_t _2N = 2 * N;

    TorusPolynomial *testvectbis = new_TorusPolynomial(N);
    MKTLweSample *acc = new_MKTLweSample(RLWEparams, MKparams);

    if (barb !=0)
    {
        torusPolynomialMulByXai(testvectbis, _2N - barb, v);
    }
    else
    {
        torusPolynomialCopy(testvectbis, v);
    }

    MKtLweNoiselessTrivial(acc, testvectbis, MKparams);
    MKtfhe_blindRotateFFT(acc, bkFFT, bara, RLWEparams, MKparams, MKrlwekey);
    MKtLweExtractMKLweSample(result, acc, MKparams);


    delete_MKTLweSample(acc);
    delete_TorusPolynomial(testvectbis);
}
















/**
 * result = MKLWE(mu) iff phase(x)>0, MKLWE(-mu) iff phase(x)<0
 * @param result The resulting MKLweSample
 * @param bk The bootstrapping + keyswitch key
 * @param mu The output message (if phase(x)>0)
 * @param x The input sample
 */
EXPORT void MKtfhe_bootstrap_woKS(MKLweSample *result, const MKLweBootstrappingKey *bk, 
        Torus32 mu, const MKLweSample *x, const TLweParams* RLWEparams, const MKTFHEParams *MKparams) 
{
    const int32_t parties = MKparams->parties;
    const int32_t N = MKparams->N;
    const int32_t Nx2 = 2 * N;
    const int32_t n = MKparams->n;

    TorusPolynomial *testvect = new_TorusPolynomial(N);
    int32_t *bara = new int32_t[parties*n];

    // b*2N
    int32_t barb = modSwitchFromTorus32(x->b, Nx2);
    // a*2N
    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            bara[n*i+j] = modSwitchFromTorus32(x->a[n*i+j], Nx2);
        }
    }
    

    //the initial testvec = [mu,mu,mu,...,mu]
    for (int32_t i = 0; i < N; i++) 
    {
        testvect->coefsT[i] = mu;
    }

    MKtfhe_blindRotateAndExtract(result, testvect, bk->bk, barb, bara, RLWEparams, MKparams);


    delete[] bara;
    delete_TorusPolynomial(testvect);
}



EXPORT void MKtfhe_bootstrap_woKSFFT(MKLweSample *result, const MKLweBootstrappingKeyFFT *bkFFT, 
        Torus32 mu, const MKLweSample *x, const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey) 
{
    const int32_t parties = MKparams->parties;
    const int32_t N = MKparams->N;
    const int32_t Nx2 = 2 * N;
    const int32_t n = MKparams->n;

    TorusPolynomial *testvect = new_TorusPolynomial(N);
    int32_t *bara = new int32_t[parties*n];

    // b*2N
    int32_t barb = modSwitchFromTorus32(x->b, Nx2);
    // a*2N
    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            bara[n*i+j] = modSwitchFromTorus32(x->a[n*i+j], Nx2);
        }
    }
    

    //the initial testvec = [mu,mu,mu,...,mu]
    for (int32_t i = 0; i < N; i++) 
    {
        testvect->coefsT[i] = mu;
    }

    MKtfhe_blindRotateAndExtractFFT(result, testvect, bkFFT->bkFFT, barb, bara, RLWEparams, MKparams, MKrlwekey);


    delete[] bara;
    delete_TorusPolynomial(testvect);
}






















/**
 * result = MKLWE(mu) iff phase(x)>0, MKLWE(-mu) iff phase(x)<0
 * @param result The resulting MKLweSample
 * @param bk The bootstrapping + keyswitch key
 * @param mu The output message (if phase(x)>0)
 * @param x The input sample
 */
EXPORT void MKtfhe_bootstrap(MKLweSample *result, const MKLweBootstrappingKey *bk, Torus32 mu, 
        const MKLweSample *x, const LweParams* LWEparams, const LweParams* extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams) 
{
    MKLweSample *u = new_MKLweSample(extractedLWEparams, MKparams);

    MKtfhe_bootstrap_woKS(u, bk, mu, x, RLWEparams, MKparams);
    // MK Key Switching
    //MKlweKeySwitch(result, bk->ks, u, MKparams);
    MKlweKeySwitch(result, bk->ks, u, LWEparams, MKparams);


    delete_MKLweSample(u);
}




EXPORT void MKtfhe_bootstrapFFT(MKLweSample *result, const MKLweBootstrappingKeyFFT *bkFFT, Torus32 mu, 
        const MKLweSample *x, const LweParams* LWEparams, const LweParams* extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey) 
{
    MKLweSample *u = new_MKLweSample(extractedLWEparams, MKparams);

    MKtfhe_bootstrap_woKSFFT(u, bkFFT, mu, x, RLWEparams, MKparams, MKrlwekey);
    // MK Key Switching
    //MKlweKeySwitch(result, bkFFT->ks, u, MKparams);
    MKlweKeySwitch(result, bkFFT->ks, u, LWEparams, MKparams);


    delete_MKLweSample(u);
}



























// Encrypt and decrypt for gate bootstrap
/** encrypts a boolean */
EXPORT void MKbootsSymEncrypt(MKLweSample *result, int32_t message, const MKLweKey* key) {

    Torus32 _1s8 = modSwitchToTorus32(1, 8);
    Torus32 mu = message ? _1s8 : -_1s8;
    double alpha = key->MKparams->stdevLWE; //TODO: specify noise

    MKlweSymEncrypt(result, mu, alpha, key);
}
/** decrypts a boolean */
EXPORT int32_t MKbootsSymDecrypt(const MKLweSample *sample, const MKLweKey* key) {

    Torus32 mu = MKlwePhase(sample, key);
    return (mu > 0 ? 1 : 0); //we have to do that because of the C binding
}














// MK Bootstrapped NAND (no FFT)
EXPORT void MKbootsNAND(MKLweSample *result, const MKLweSample *ca, const MKLweSample *cb, 
        const MKLweBootstrappingKey *bk, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams) 
{
    static const Torus32 MU = modSwitchToTorus32(1, 8);

    MKLweSample *temp_result = new_MKLweSample(LWEparams, MKparams);

    //compute: (0,1/8) - ca - cb
    static const Torus32 NandConst = modSwitchToTorus32(1, 8);
    MKlweNoiselessTrivial(temp_result, NandConst, MKparams);
    MKlweSubTo(temp_result, ca, MKparams);
    MKlweSubTo(temp_result, cb, MKparams);


    //if the phase is positive, the result is 1/8
    //if the phase is positive, else the result is -1/8
    MKtfhe_bootstrap(result, bk, MU, temp_result, LWEparams, extractedLWEparams, RLWEparams, MKparams);
   

    delete_MKLweSample(temp_result);
}



EXPORT void MKbootsNAND_FFT(MKLweSample *result, const MKLweSample *ca, const MKLweSample *cb, 
        const MKLweBootstrappingKeyFFT *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey) 
{
    static const Torus32 MU = modSwitchToTorus32(1, 8);

    MKLweSample *temp_result = new_MKLweSample(LWEparams, MKparams);

    //compute: (0,1/8) - ca - cb
    static const Torus32 NandConst = modSwitchToTorus32(1, 8);
    MKlweNoiselessTrivial(temp_result, NandConst, MKparams);
    MKlweSubTo(temp_result, ca, MKparams);
    MKlweSubTo(temp_result, cb, MKparams);


    //if the phase is positive, the result is 1/8
    //if the phase is positive, else the result is -1/8
    MKtfhe_bootstrapFFT(result, bkFFT, MU, temp_result, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey); 
       

    delete_MKLweSample(temp_result);
}
























































/* **************************************************************************
***************************** VERSION 2 *************************************
************************************************************************** */




/* *************************
******** MK-RGSW ***********
************************* */


/* Uni-Encrypt */
// Encrypt a integer polynomial as (d, F) = (d, f0, f1)
EXPORT void MKTGswUniEncrypt_v2(MKTGswUESample_v2 *result, IntPolynomial *message, int32_t party, double alpha, const MKRLweKey *key) {
    const int32_t N = key->RLWEparams->N;
    const int32_t dg = key->MKparams->dg;
    const int32_t parties = key->MKparams->parties;

    // generate r, the randomness
    uniform_int_distribution<int32_t> distribution(0, 1);
    IntPolynomial* r = new_IntPolynomial(N);
    for (int j = 0; j < N; ++j){
        r->coefs[j] = distribution(generator);
    }
            

    // d = r*Pkey_parties + m*g + E1 \in T^dg
    for (int i = 0; i < dg; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            // d = E1
            result->d[i].coefsT[j] = gaussian32(0, alpha); // E1
            // d = E1 + m*g[i]
            result->d[i].coefsT[j] += message->coefs[j] * key->MKparams->g[i]; // m*g[i]
        }

        // d = r*Pkey_parties[i] + E1 + m*g[i]
        torusPolynomialAddMulR(&result->d[i], r, &key->Pkey[parties*dg + i]);   
    }


    // F = (f0,f1) \in T^2dg, with f0 = s_party*f1 + e_f + r*g
    for (int i = 0; i < dg; ++i)
    {
        // f1 
        torusPolynomialUniform(&result->f1[i]); 

        // f0 = e_f[i] + r*g[i]
        for (int j = 0; j < N; ++j)
        {
            result->f0[i].coefsT[j] = gaussian32(0, alpha); // e_f
            result->f0[i].coefsT[j] += r->coefs[j] * key->MKparams->g[i]; // r*g[i]
        }

        torusPolynomialAddMulR(&result->f0[i], key->key[party].key, &result->f1[i]);       
    }
    
   
    result->current_variance = alpha * alpha;
    delete_IntPolynomial(r);
}







// Encrypt an integer value as (d, F) = (d, f0, f1)
EXPORT void MKTGswUniEncryptI_v2(MKTGswUESample_v2 *result, int32_t message, int32_t party, double alpha, const MKRLweKey *key) {
    // cout << "alpha for uniencrypt = " << alpha << endl;
    const int32_t N = key->RLWEparams->N;
    const int32_t dg = key->MKparams->dg;
    const int32_t parties = key->MKparams->parties;

    // generate r, the randomness
    uniform_int_distribution<int32_t> distribution(0, 1);
    IntPolynomial* r = new_IntPolynomial(N);

    
    cout << "r = "; // just for verification of f part
    for (int j = 0; j < N; ++j){
        r->coefs[j] = distribution(generator);
        cout << r->coefs[j] << " "; // just for verification of f part 
    }
    cout << endl; // just for verification of f part 
    
   



    // d = r*Pkey_parties + m*g + E1 \in T^dg
    for (int i = 0; i < dg; ++i)
    {
        //cout << "E:";
        for (int j = 0; j < N; ++j)
        {
            // d = E1
            result->d[i].coefsT[j] = gaussian32(0, alpha); // E1
            // cout << result->d[i].coefsT[j] << ", "; 
        }
        //cout << endl;
        // d = E1 + m*g[i]
        result->d[i].coefsT[0] += message * key->MKparams->g[i]; // m*g[i]

        // d1 = r*Pkey_parties[i] + E1 + m*g[i] 
        torusPolynomialAddMulR(&result->d[i], r, &key->Pkey[dg*parties + i]); 
        //cout << "d1 vec = ";
        //for(int j = 0; j < N; j++){
        //    cout << result->d[dg+i].coefsT[j] << ", "; 
        // }
        //cout << endl;  
    }


    // F = (f0,f1) \in T^2dg, with f0 = s_party*f1 + e_f + r*g
    for (int i = 0; i < dg; ++i)
    {
        // f1 
        torusPolynomialUniform(&result->f1[i]); 

        // f0 = e_f[i] + r*g[i]
        for (int j = 0; j < N; ++j)
        {
            result->f0[i].coefsT[j] = gaussian32(0, alpha); // e_f
            result->f0[i].coefsT[j] += r->coefs[j] * key->MKparams->g[i]; // r*g[i]
        }
        // f0 = s_party*f1 + e_f + r*g
        torusPolynomialAddMulR(&result->f0[i], key->key[party].key, &result->f1[i]);       
    }
      

    result->current_variance = alpha * alpha;
    delete_IntPolynomial(r);
}










/**
 * This function computes the decryption (actually the phase) of sample by using key
 * The constant Msize indicates the message space and is used to approximate the phase
 */
// result is an array composed by dg torus polynomials ~r*g[j]
EXPORT void MKtGswSymDecrypt_v2(TorusPolynomial *result, const MKTGswUESample_v2 *sample, const MKRLweKey *key) {
    const int32_t dg = key->MKparams->dg;
    const int32_t party = sample->party;
    

    for (int j = 0; j < dg; ++j)
    {
        // f part 
        torusPolynomialCopy(&result[j], &sample->f0[j]); // phi = f0[j] 
        // phi = f0 - f1*s_party 
        torusPolynomialSubMulR(&result[j], key->key[party].key, &sample->f1[j]);
    }

}










/* EXPAND */
// (d,F) = (d,f0,f1) -> D_i=(x_0, ..., x_{parties-1}, x_parties + d_i, y_0, ..., d_i+y_i, ..., y_parties, d_i)
EXPORT void MKTGswExpand_v2(MKTGswExpSample_v2 *result, const MKTGswUESample_v2 *sample, const MKRLweKey *key, 
    const MKTFHEParams* MKparams) 
{
    const int32_t N = key->RLWEparams->N;
    const int32_t dg = key->MKparams->dg;
    const int32_t party = sample->party;
    const int32_t parties = key->MKparams->parties;


    // INITIALIZE
    // D_i=(0, ..., 0, d_i, 0, ..., d_i, ..., 0, d_i)
    for (int j = 0; j < dg; ++j)
    {
        // initialize x_0, ..., x_{parties-1} as 0
        // initialize x_parties as d_i (x[parties*dg +j] = d[j])
        for (int i = 0; i < parties; ++i)
        {
            torusPolynomialClearN(&result->x[i*dg + j], N);
        }
        torusPolynomialCopyN(&result->x[parties*dg + j], &sample->d[j], N);
        
        // initialize all the y_i as 0
        // initialize y_party as d_i (y[party*dg +j] = d[j])
        for (int i = 0; i <= parties; ++i)
        {
            torusPolynomialClearN(&result->y[i*dg + j], N);
        }
        torusPolynomialCopyN(&result->y[party*dg + j], &sample->d[j], N);

        // d_i = d_i (d[j] = d[j])
        torusPolynomialCopyN(&result->d[j], &sample->d[j], N);
    }


    TorusPolynomial* X = new_TorusPolynomial(N);
    TorusPolynomial* Y = new_TorusPolynomial(N);
    IntPolynomial* u = new_IntPolynomial_array(dg, N);

    // i < parties 
    for (int i = 0; i < parties; ++i) 
    {
        for (int j = 0; j < dg; ++j)
        {
            // g^{-1}(b_i[j]) = [u_0, ...,u_dg-1] intPolynomials
            MKtGswTorus32PolynomialDecompGassembly(u, &key->Pkey[i*dg + j], MKparams);

            // X=0 and Y=0
            torusPolynomialClearN(X, N);
            torusPolynomialClearN(Y, N);
            for (int l = 0; l < dg; ++l)
            {
                // X = x_i[j] = <g^{-1}(b_i[j]), f0>
                torusPolynomialAddMulRFFTN(X, &u[l], &sample->f0[l], N);
                // Y = y_i[j] = <g^{-1}(b_i[j]), f1> 
                torusPolynomialAddMulRFFTN(Y, &u[l], &sample->f1[l], N);          
            }
            
            // x_i
            torusPolynomialAddTo1(&result->x[i*dg + j], X); // N = X->N
            // y_i
            torusPolynomialAddTo1(&result->y[i*dg + j], Y); // N = Y->N
        }   
    }
    
    // i = parties, i.e. b_i = a
    for (int j = 0; j < dg; ++j)
    {
        // g^{-1}(a[j]) = [u_0, ...,u_dg-1] intPolynomials
        MKtGswTorus32PolynomialDecompGassembly(u, &key->Pkey[parties*dg + j], MKparams);

        // X=0 and Y=0
        torusPolynomialClearN(X, N);
        torusPolynomialClearN(Y, N);
        for (int l = 0; l < dg; ++l)
        {
            // X = x_i[j] = <g^{-1}(a[j]), f0>
            torusPolynomialAddMulRFFTN(X, &u[l], &sample->f0[l], N);
            // Y = y_i[j] = <g^{-1}(a[j]), f1> 
            torusPolynomialAddMulRFFTN(Y, &u[l], &sample->f1[l], N);          
        }
        
        // x_i
        torusPolynomialSubTo1(&result->x[parties*dg + j], X); // N = X->N
        // y_i
        torusPolynomialSubTo1(&result->y[parties*dg + j], Y); // N = Y->N
    }   


    result->party = sample->party;
    // TODO: fix this
    result->current_variance = sample->current_variance;
    
    delete_TorusPolynomial(X);
    delete_TorusPolynomial(Y);
    delete_IntPolynomial_array(dg, u);
}






/* EXPAND */
// (d,F) = (d,f0,f1) -> D_i=(x_0, ..., x_{parties-1}, x_parties + d_i, y_0, ..., d_i+y_i, ..., y_perties, d_i)
// sample UE --> resultFFT expand
EXPORT void MKTGswExpandFFT_v2(MKTGswExpSampleFFT_v2 *resultFFT, const MKTGswUESampleFFT_v2* sampleFFT, const MKRLweKey *key, 
        const TLweParams* RLWEparams, const MKTFHEParams* MKparams) 
{
    const int32_t N = key->RLWEparams->N;
    const int32_t dg = key->MKparams->dg;
    const int32_t party = sampleFFT->party;
    const int32_t parties = key->MKparams->parties;


    LagrangeHalfCPolynomial* tempFFT = new_LagrangeHalfCPolynomial(N);

    
    // INITIALIZE
    // D_i=(0, ..., 0, d_i, 0, ..., d_i, ..., 0, d_i)
    for (int j = 0; j < dg; ++j)
    {
        // initialize x_0, ..., x_{parties-1} as 0
        // initialize x_parties as d_i (x[parties*dg +j] = d[j])
        for (int i = 0; i < parties; ++i)
        {
            LagrangeHalfCPolynomialClear(&resultFFT->x[i*dg + j]);
        }
        LagrangeHalfCPolynomialCopy(&resultFFT->x[parties*dg + j], &sampleFFT->d[j]);

        // initialize all the y_i as 0
        // initialize y_party as d_i (y[party*dg +j] = d[j])
        for (int i = 0; i <= parties; ++i)
        {
            LagrangeHalfCPolynomialClear(&resultFFT->y[i*dg + j]);
        }
        LagrangeHalfCPolynomialCopy(&resultFFT->y[party*dg + j], &sampleFFT->d[j]);
        
        // d_i = d_i (d[j] = d[j])
        LagrangeHalfCPolynomialCopy(&resultFFT->d[j], &sampleFFT->d[j]);
    }






    LagrangeHalfCPolynomial* X = new_LagrangeHalfCPolynomial(N);
    LagrangeHalfCPolynomial* Y = new_LagrangeHalfCPolynomial(N);
    IntPolynomial* u = new_IntPolynomial_array(dg, N);
    LagrangeHalfCPolynomial *uFFT = new_LagrangeHalfCPolynomial_array(dg, N); //fft version


    // i < parties 
    for (int i = 0; i < parties; ++i) 
    {
        for (int j = 0; j < dg; ++j)
        {
            // g^{-1}(b_i[j]) = [u_0, ...,u_dg-1] intPolynomials
            MKtGswTorus32PolynomialDecompGassembly(u, &key->Pkey[i*dg + j], MKparams);
            for (int p = 0; p < dg; ++p){
                IntPolynomial_ifft(&uFFT[p], &u[p]); // FFT
            }

            // X=0 and Y=0
            LagrangeHalfCPolynomialClear(X);
            LagrangeHalfCPolynomialClear(Y);
            for (int l = 0; l < dg; ++l)
            {
                // X = xi[j] = <g^{-1}(b_i[j]), f0>
                LagrangeHalfCPolynomialMul(tempFFT, &uFFT[l], &sampleFFT->f0[l]);
                LagrangeHalfCPolynomialAddTo(X, tempFFT);
                // Y = yi[j] = <g^{-1}(b_i[j]), f1> 
                LagrangeHalfCPolynomialMul(tempFFT, &uFFT[l], &sampleFFT->f1[l]);
                LagrangeHalfCPolynomialAddTo(Y, tempFFT);         
            }
            
            // x_i
            LagrangeHalfCPolynomialAddTo(&resultFFT->x[i*dg + j], X);
            // y_i
            LagrangeHalfCPolynomialAddTo(&resultFFT->y[i*dg + j], Y);
        }   
    }

    // i = parties, i.e. b_i = a 
    for (int j = 0; j < dg; ++j)
    {
        // g^{-1}(a[j]) = [u_0, ...,u_dg-1] intPolynomials
        MKtGswTorus32PolynomialDecompGassembly(u, &key->Pkey[parties*dg + j], MKparams);
        for (int p = 0; p < dg; ++p){
            IntPolynomial_ifft(&uFFT[p], &u[p]); // FFT
        }

        // X=0 and Y=0
        LagrangeHalfCPolynomialClear(X);
        LagrangeHalfCPolynomialClear(Y);
        for (int l = 0; l < dg; ++l)
        {
            // X = xi[j] = <g^{-1}(a[j]), f0>
            LagrangeHalfCPolynomialMul(tempFFT, &uFFT[l], &sampleFFT->f0[l]);
            LagrangeHalfCPolynomialAddTo(X, tempFFT);
            // Y = yi[j] = <g^{-1}(a[j]), f1> 
            LagrangeHalfCPolynomialMul(tempFFT, &uFFT[l], &sampleFFT->f1[l]);
            LagrangeHalfCPolynomialAddTo(Y, tempFFT);         
        }
        
        // x_i
        LagrangeHalfCPolynomialSubTo(&resultFFT->x[parties*dg + j], X);
        // y_i
        LagrangeHalfCPolynomialSubTo(&resultFFT->y[parties*dg + j], Y);
    }   
    





    resultFFT->party = sampleFFT->party;
    // TODO: fix this
    resultFFT->current_variance = sampleFFT->current_variance;
    


    // delete 
    delete_LagrangeHalfCPolynomial_array(dg, uFFT);
    delete_IntPolynomial_array(dg, u);
    delete_LagrangeHalfCPolynomial(Y);
    delete_LagrangeHalfCPolynomial(X);
    delete_LagrangeHalfCPolynomial(tempFFT);
   
}











/* ********************************************************************************
*********************** EXTERNAL PRODUCT method 1 *********************************
******************************************************************************** */


// c' = G^{-1}(c)*C, with C = (d, F) = (d, f0, f1) 
EXPORT void MKtGswUEExternMulToMKtLwe_v2m1(MKTLweSample* result, MKTLweSample* sample, 
        const MKTGswUESample_v2* sampleUE, 
        const TLweParams* RLWEparams,
        const MKTFHEParams* MKparams,
        const MKRLweKey *RLWEkey)
{
    const int32_t N = MKparams->N;
    const int32_t dg = MKparams->dg;
    const int32_t party = sampleUE->party;
    const int32_t parties = MKparams->parties;


    // Expand UE sample 
    MKTGswExpSample_v2* sampleExp = new_MKTGswExpSample_v2(RLWEparams, MKparams);
    MKTGswExpand_v2(sampleExp, sampleUE, RLWEkey, MKparams);


    // DECOMPOSE sample
    // u[i] = g^{-1}(a_i), u[parties] = g^{-1}(b), 
    IntPolynomial* u = new_IntPolynomial_array((parties+1)*dg, N);
    for (int i = 0; i <= parties; ++i)
    {
        MKtGswTorus32PolynomialDecompGassembly(&u[i*dg], &sample->a[i], MKparams);
    }




    // c'_i = g^{-1}(a_i)*d, i<parties, i!=party
    for (int i = 0; i < party; ++i)
    {
        torusPolynomialClearN(&result->a[i], N);
        for (int l = 0; l < dg; ++l)
        {
            torusPolynomialAddMulRFFTN(&result->a[i], &u[i*dg+l], &sampleExp->d[l], N);
        }
    }
    for (int i = party+1; i < parties; ++i)
    {
        torusPolynomialClearN(&result->a[i], N);
        for (int l = 0; l < dg; ++l)
        {
            torusPolynomialAddMulRFFTN(&result->a[i], &u[i*dg+l], &sampleExp->d[l], N);
        }
    }

    // c'_party = \sum g^{-1}(a_i)*yi 
    torusPolynomialClearN(&result->a[party], N);
    for (int i = 0; i <= parties; ++i)
    {
        for (int l = 0; l < dg; ++l)
        {
            torusPolynomialAddMulRFFTN(&result->a[party], &u[i*dg + l], &sampleExp->y[dg*i + l], N);
        }
    }

    // c'_parties = \sum g^{-1}(a_i)*xi
    torusPolynomialClearN(&result->a[parties], N);
    for (int i = 0; i <= parties; ++i)
    {
        for (int l = 0; l < dg; ++l)
        {
            torusPolynomialAddMulRFFTN(&result->a[parties], &u[i*dg + l], &sampleExp->x[dg*i + l], N);
        }
    }



    // TODO current_variance

    delete_MKTGswExpSample_v2(sampleExp);
    delete_IntPolynomial_array((parties+1)*dg, u);
}










// c' = G^{-1}(c)*C, with C = (d, F) = (d, f0, f1) 
// result is not in FFT
EXPORT void MKtGswUEExternMulToMKtLwe_FFT_v2m1(MKTLweSample* result, MKTLweSample* sample, 
        const MKTGswUESampleFFT_v2* sampleUEFFT, 
        const TLweParams* RLWEparams,
        const MKTFHEParams* MKparams,
        const MKRLweKey *RLWEkey)
{
    const int32_t N = MKparams->N;
    const int32_t dg = MKparams->dg;
    const int32_t party = sampleUEFFT->party;
    const int32_t parties = MKparams->parties;
    const int parties1dg = (parties+1)*dg;

    


    // TEST
     // Expand UE sample 
    MKTGswExpSample_v2* sampleExp = new_MKTGswExpSample_v2(RLWEparams, MKparams);
    MKTGswUESample_v2* sampleUE = new_MKTGswUESample_v2(RLWEparams, MKparams);

    for (int i = 0; i < 3*dg; ++i)
    {
        TorusPolynomial_fft(&sampleUE->d[i], &sampleUEFFT->d[i]);
    } 

    MKTGswExpand_v2(sampleExp, sampleUE, RLWEkey, MKparams);

    MKTGswExpSampleFFT_v2* sampleExpFFT = new_MKTGswExpSampleFFT_v2(RLWEparams,  MKparams, sampleUEFFT->current_variance);
    for (int i = 0; i < (2*(parties+1)+1)*dg; ++i)
    {
        TorusPolynomial_ifft(&sampleExpFFT->x[i], &sampleExp->x[i]);
    }
    // TEST end 

/*
    // Expand UE sample in FFT 
    MKTGswExpSampleFFT_v2* sampleExpFFT = new_MKTGswExpSampleFFT_v2(RLWEparams,  MKparams, sampleUEFFT->current_variance);
    MKTGswExpandFFT_v2(sampleExpFFT, sampleUEFFT, RLWEkey, RLWEparams, MKparams);
*/
    






    // DECOMPOSE sample and convert it to FFT
    // u[i*dg] = g^{-1}(a_i), u[parties*dg] = g^{-1}(b), 
    IntPolynomial* uDec = new_IntPolynomial_array(parties1dg, N);
    LagrangeHalfCPolynomial *uDecFFT = new_LagrangeHalfCPolynomial_array(parties1dg, N); //fft version
    for (int i = 0; i <= parties; ++i){
        MKtGswTorus32PolynomialDecompGassembly(&uDec[i*dg], &sample->a[i], MKparams);
    }
    for (int p = 0; p < parties1dg; ++p){
        IntPolynomial_ifft(&uDecFFT[p], &uDec[p]); // FFT
    }


    // c'_i = g^{-1}(a_i)*d, i<parties, i!=party
    for (int i = 0; i < party; ++i)
    {
        torusPolynomialClearN(&result->a[i], N);
        for (int l = 0; l < dg; ++l)
        {
            MulFFTAndAddTo(&result->a[i], &uDecFFT[i*dg+l], &sampleExpFFT->d[l], N);
        }
    }
    for (int i = party+1; i < parties; ++i)
    {
        torusPolynomialClearN(&result->a[i], N);
        for (int l = 0; l < dg; ++l)
        {
            MulFFTAndAddTo(&result->a[i], &uDecFFT[i*dg+l], &sampleExpFFT->d[l], N);
        }
    }

    
    // c'_party = \sum g^{-1}(a_i)*yi
    torusPolynomialClearN(&result->a[party], N);
    for (int i = 0; i <= parties; ++i)
    {
        for (int l = 0; l < dg; ++l)
        {
            MulFFTAndAddTo(&result->a[party], &uDecFFT[i*dg+l], &sampleExpFFT->y[dg*i + l], N);
        }
    }


    // c'_parties = \sum g^{-1}(a_i)*xi
    torusPolynomialClearN(&result->a[parties], N);
    for (int i = 0; i <= parties; ++i)
    {
        for (int l = 0; l < dg; ++l)
        {
            MulFFTAndAddTo(&result->a[parties], &uDecFFT[i*dg+l], &sampleExpFFT->x[dg*i + l], N);
        }
    }

    
    // TODO current_variance
    delete_LagrangeHalfCPolynomial_array(parties1dg, uDecFFT);
    delete_IntPolynomial_array(parties1dg, uDec);
    delete_MKTGswExpSampleFFT_v2(sampleExpFFT);
}

























/* ********************************************************************************
*********************** EXTERNAL PRODUCT method 2 *********************************
******************************************************************************** */


// c' = G^{-1}(c)*C, with C = (d, F) = (d, f0, f1) 
EXPORT void MKtGswUEExternMulToMKtLwe_v2m2(MKTLweSample* result, MKTLweSample* sample, 
        const MKTGswUESample_v2* sampleUE, 
        const TLweParams* RLWEparams,
        const MKTFHEParams* MKparams,
        const MKRLweKey *RLWEkey)
{
    const int32_t N = MKparams->N;
    const int32_t dg = MKparams->dg;
    const int32_t party = sampleUE->party;
    const int32_t parties = MKparams->parties;


    // DECOMPOSE sample
    // uDec[i] = g^{-1}(a_i), uDec[parties] = g^{-1}(b), 
    IntPolynomial* uDec = new_IntPolynomial_array((parties+1)*dg, N);
    for (int i = 0; i <= parties; ++i)
    {
        MKtGswTorus32PolynomialDecompGassembly(&uDec[i*dg], &sample->a[i], MKparams);
    }



    // u[i] = uDec[i] * d 
    TorusPolynomial* u = new_TorusPolynomial_array(parties+1, N);
    for (int i = 0; i <= parties; ++i)
    {
        torusPolynomialClearN(&u[i], N);
        for (int j = 0; j < dg; ++j)
        {
            torusPolynomialAddMulRFFTN(&u[i], &uDec[i*dg+j], &sampleUE->d[j], N);
        }
    }

    // v[i] = uDec[i] * b_i, for i < parties  
    TorusPolynomial* v = new_TorusPolynomial_array(parties+1, N);
    for (int i = 0; i < parties; ++i)
    {
        torusPolynomialClearN(&v[i], N);
        for (int j = 0; j < dg; ++j)
        {
            torusPolynomialAddMulRFFTN(&v[i], &uDec[i*dg+j], &RLWEkey->Pkey[i*dg + j], N);
        }
    }
    // v[parties] = - uDec[parties] * a
    torusPolynomialClearN(&v[parties], N);
    for (int j = 0; j < dg; ++j)
    {
        torusPolynomialSubMulRFFTN(&v[parties], &uDec[parties*dg+j], &RLWEkey->Pkey[parties*dg + j], N);
    }
    


    // Decompose v
    // vDec[i] = g^{-1}(v[i]) 
    IntPolynomial* vDec = new_IntPolynomial_array((parties+1)*dg, N);
    for (int i = 0; i <= parties; ++i)
    {
        MKtGswTorus32PolynomialDecompGassembly(&vDec[i*dg], &v[i], MKparams);
    }


    // w0[i] = vDec[i] * f0 
    TorusPolynomial* w0 = new_TorusPolynomial_array(parties+1, N);
    for (int i = 0; i <= parties; ++i)
    {
        torusPolynomialClearN(&w0[i], N);
        for (int j = 0; j < dg; ++j)
        {
            torusPolynomialAddMulRFFTN(&w0[i], &vDec[i*dg+j], &sampleUE->f0[j], N);
        }
    }
    // w1[i] = vDec[i] * f1 
    TorusPolynomial* w1 = new_TorusPolynomial_array(parties+1, N);
    for (int i = 0; i <= parties; ++i)
    {
        torusPolynomialClearN(&w1[i], N);
        for (int j = 0; j < dg; ++j)
        {
            torusPolynomialAddMulRFFTN(&w1[i], &vDec[i*dg+j], &sampleUE->f1[j], N);
        }
    }





    // c'_i = u[i], i<parties, i!=party
    for (int i = 0; i < party; ++i)
    {
        torusPolynomialCopyN(&result->a[i], &u[i], N);
        
    }
    for (int i = party+1; i < parties; ++i)
    {
        torusPolynomialCopyN(&result->a[i], &u[i], N);
    }


    // c'_party = u[party] + \sum w1[i] 
    torusPolynomialCopyN(&result->a[party], &u[party], N);
    for (int i = 0; i <= parties; ++i)
    {
        torusPolynomialAddTo1(&result->a[party], &w1[i]);
    }


    // c'_parties = u[parties] + \sum w0[i] 
    torusPolynomialCopyN(&result->a[parties], &u[parties], N);
    for (int i = 0; i <= parties; ++i)
    {
        torusPolynomialAddTo1(&result->a[parties], &w0[i]);
    }


    // TODO current_variance
    delete_TorusPolynomial_array(parties+1, w1);
    delete_TorusPolynomial_array(parties+1, w0);
    delete_IntPolynomial_array((parties+1)*dg, vDec);
    delete_TorusPolynomial_array(parties+1, v);
    delete_TorusPolynomial_array(parties+1, u);
    delete_IntPolynomial_array((parties+1)*dg, uDec);
}









// c' = G^{-1}(c)*C, with C = (d, F) = (d, f0, f1) 
// result is not in FFT
EXPORT void MKtGswUEExternMulToMKtLwe_FFT_v2m2(MKTLweSample* result, MKTLweSample* sample, 
        const MKTGswUESampleFFT_v2* sampleUEFFT, 
        const TLweParams* RLWEparams,
        const MKTFHEParams* MKparams,
        const MKRLweKey *RLWEkey)
{
    const int32_t N = MKparams->N;
    const int32_t dg = MKparams->dg;
    const int32_t party = sampleUEFFT->party;
    const int32_t parties = MKparams->parties;
    const int parties1dg = (parties+1)*dg;


    // DECOMPOSE sample and convert it to FFT
    // uDec[i*dg] = g^{-1}(a_i), uDec[parties*dg] = g^{-1}(b), 
    IntPolynomial* uDec = new_IntPolynomial_array(parties1dg, N);
    LagrangeHalfCPolynomial *uDecFFT = new_LagrangeHalfCPolynomial_array(parties1dg, N); //fft version
    for (int i = 0; i <= parties; ++i){
        MKtGswTorus32PolynomialDecompGassembly(&uDec[i*dg], &sample->a[i], MKparams);
    }
    for (int p = 0; p < parties1dg; ++p){
        IntPolynomial_ifft(&uDecFFT[p], &uDec[p]); // FFT
    }


    

    

    // uFFT[i] = uDecFFT[i] * dFFT
    // LagrangeHalfCPolynomial *uFFT = new_LagrangeHalfCPolynomial_array(parties+1, N); 
    TorusPolynomial *u = new_TorusPolynomial_array(parties+1, N);
    for (int i = 0; i <= parties; ++i)
    {
        // LagrangeHalfCPolynomialClear(&uFFT[i]);
        torusPolynomialClearN(&u[i], N);
        for (int j = 0; j < dg; ++j)

        {
            MulFFTAndAddTo(&u[i], &uDecFFT[i*dg+j], &sampleUEFFT->d[j], N);
            // LagrangeHalfCPolynomialClear(tempFFT);
            // LagrangeHalfCPolynomialMul(tempFFT, &uDecFFT[i*dg+j], &sampleUEFFT->d[j]);
            // LagrangeHalfCPolynomialAddTo(&uFFT[i], tempFFT);
        }
    }


    // computed non in FFT because it needs to be decomposed
    // v[i] = uDec[i] * b_i, for i < parties  
    TorusPolynomial* v = new_TorusPolynomial_array(parties+1, N);
    for (int i = 0; i < parties; ++i)
    {
        torusPolynomialClearN(&v[i], N);
        for (int j = 0; j < dg; ++j)
        {
            torusPolynomialAddMulRFFTN(&v[i], &uDec[i*dg+j], &RLWEkey->Pkey[i*dg + j], N);
        }
    }
    // v[parties] = - uDec[parties] * a
    torusPolynomialClearN(&v[parties], N);
    for (int j = 0; j < dg; ++j)
    {
        torusPolynomialSubMulRFFTN(&v[parties], &uDec[parties*dg+j], &RLWEkey->Pkey[parties*dg + j], N);
    }
    // Decompose v and convert it in FFT
    // vDec[i] = g^{-1}(v[i]) 
    IntPolynomial* vDec = new_IntPolynomial_array(parties1dg, N);
    LagrangeHalfCPolynomial *vDecFFT = new_LagrangeHalfCPolynomial_array(parties1dg, N); //fft version
    for (int i = 0; i <= parties; ++i)
    {
        MKtGswTorus32PolynomialDecompGassembly(&vDec[i*dg], &v[i], MKparams);
    }
    for (int p = 0; p < parties1dg; ++p){
        IntPolynomial_ifft(&vDecFFT[p], &vDec[p]); // FFT
    } 


    // w0FFT[i] = vDecFFT[i] * f0FFT
    //LagrangeHalfCPolynomial *w0FFT = new_LagrangeHalfCPolynomial_array(parties+1, N); 
    TorusPolynomial *w0 = new_TorusPolynomial_array(parties+1, N);
    for (int i = 0; i <= parties; ++i)
    {
        // LagrangeHalfCPolynomialClear(&w0FFT[i]);
        torusPolynomialClearN(&w0[i], N);
        for (int j = 0; j < dg; ++j)
        {
            MulFFTAndAddTo(&w0[i], &vDecFFT[i*dg+j], &sampleUEFFT->f0[j], N);
            // LagrangeHalfCPolynomialClear(tempFFT);
            // LagrangeHalfCPolynomialMul(tempFFT, &vDecFFT[i*dg+j], &sampleUEFFT->f0[j]);
            // LagrangeHalfCPolynomialAddTo(&w0FFT[i], tempFFT);
        }
    }
    // w1FFT[i] = vDecFFT[i] * f1FFT
    // LagrangeHalfCPolynomial *w1FFT = new_LagrangeHalfCPolynomial_array(parties+1, N); 
    TorusPolynomial *w1 = new_TorusPolynomial_array(parties+1, N);
    for (int i = 0; i <= parties; ++i)
    {
        // LagrangeHalfCPolynomialClear(&w1FFT[i]);
        torusPolynomialClearN(&w1[i], N);
        for (int j = 0; j < dg; ++j)
        {
            MulFFTAndAddTo(&w1[i], &vDecFFT[i*dg+j], &sampleUEFFT->f1[j], N);
            // LagrangeHalfCPolynomialClear(tempFFT);
            // LagrangeHalfCPolynomialMul(tempFFT, &vDecFFT[i*dg+j], &sampleUEFFT->f1[j]);
            // LagrangeHalfCPolynomialAddTo(&w1FFT[i], tempFFT);
        }
    }







    // the result is not in FFT
    // c'_i = invFFT(uFFT[i]), i<parties, i!=party
    for (int i = 0; i < party; ++i)
    {
        torusPolynomialCopyN(&result->a[i], &u[i], N);
        //TorusPolynomial_fft(&result->a[i], &uFFT[i]); // invFFT        
    }
    for (int i = party+1; i < parties; ++i)
    {
        torusPolynomialCopyN(&result->a[i], &u[i], N);
        //TorusPolynomial_fft(&result->a[i], &uFFT[i]); // invFFT 
    }



    // c'_party = invFFT( uFFT[party] + \sum w1FFT[i] ) 
    // LagrangeHalfCPolynomialClear(tempFFT);
    // LagrangeHalfCPolynomialCopy(tempFFT, &uFFT[party]);
    torusPolynomialCopyN(&result->a[party], &u[party], N);
    for (int i = 0; i <= parties; ++i)
    {
        torusPolynomialAddTo1(&result->a[party], &w1[i]);
        // LagrangeHalfCPolynomialAddTo(tempFFT, &w1FFT[i]);
    }
    // TorusPolynomial_fft(&result->a[party], tempFFT); // invFFT  



    // c'_parties = invFFT( uFFT[parties] + \sum w0FFT[i] ) 
    // LagrangeHalfCPolynomialClear(tempFFT);
    // LagrangeHalfCPolynomialCopy(tempFFT, &uFFT[parties]);
    torusPolynomialCopyN(&result->a[parties], &u[parties], N);
    for (int i = 0; i <= parties; ++i)
    {
        torusPolynomialAddTo1(&result->a[parties], &w0[i]);
        // LagrangeHalfCPolynomialAddTo(tempFFT, &w0FFT[i]);
    }
    //TorusPolynomial_fft(&result->a[parties], tempFFT); // invFFT  



    // TODO current_variance   
    // delete_LagrangeHalfCPolynomial_array(parties+1, w1FFT); 
    // delete_LagrangeHalfCPolynomial_array(parties+1, w0FFT); 
    delete_TorusPolynomial_array(parties+1, w1); 
    delete_TorusPolynomial_array(parties+1, w0); 
    delete_LagrangeHalfCPolynomial_array(parties1dg, vDecFFT); 
    delete_IntPolynomial_array(parties1dg, vDec);
    delete_TorusPolynomial_array(parties+1, v);
    delete_TorusPolynomial_array(parties+1, u);
    // delete_LagrangeHalfCPolynomial(tempFFT);
    // delete_LagrangeHalfCPolynomial_array(parties1dg, uDecFFT);
    delete_LagrangeHalfCPolynomial_array(parties1dg, uDecFFT);
    delete_IntPolynomial_array(parties1dg, uDec);
}


    

    


   

























































/* ********************************************************************************
************************** Bootstrapping method 1 *********************************
******************************************************************************** */


// MUX -> rotate
// Only the PK part of RLWEkey is used 
void MKtfhe_MuxRotate_v2m1(MKTLweSample *result, MKTLweSample *accum, const MKTGswUESample_v2* bki, 
    const int32_t barai, const TLweParams* RLWEparams, const MKTFHEParams* MKparams, const MKRLweKey *RLWEkey) 
{
    MKTLweSample *temp_result = new_MKTLweSample(RLWEparams, MKparams);

    // ACC = BKi*[(X^barai-1)*ACC]+ACC
    // temp = (X^barai-1)*ACC
    MKtLweMulByXaiMinusOne(temp_result, barai, accum, MKparams);
    // temp *= BKi
    MKtGswUEExternMulToMKtLwe_v2m1(result, temp_result, bki, RLWEparams, MKparams, RLWEkey);
    // ACC += temp
    MKtLweAddTo(result, accum, MKparams);

    delete_MKTLweSample(temp_result);
}









// MUX -> rotate
// Only the PK part of RLWEkey is used 
void MKtfhe_MuxRotateFFT_v2m1(MKTLweSample *result, MKTLweSample *accum, const MKTGswUESampleFFT_v2 *bkiFFT, 
    const int32_t barai, const TLweParams* RLWEparams, const MKTFHEParams* MKparams, const MKRLweKey *RLWEkey) 
{
    MKTLweSample *temp_result = new_MKTLweSample(RLWEparams, MKparams);

    // ACC = BKi*[(X^barai-1)*ACC]+ACC
    // temp = (X^barai-1)*ACC
    MKtLweMulByXaiMinusOne(temp_result, barai, accum, MKparams);
    // temp *= BKi
    MKtGswUEExternMulToMKtLwe_FFT_v2m1(result, temp_result,bkiFFT, RLWEparams, MKparams, RLWEkey);
    // ACC += temp
    MKtLweAddTo(result, accum, MKparams);

    delete_MKTLweSample(temp_result);
}










// MK Blind rotate
// Only the PK part of RLWEkey is used 
EXPORT void MKtfhe_blindRotate_v2m1(MKTLweSample *accum, const MKTGswUESample_v2 *bk, const int32_t *bara, 
    const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *RLWEkey) 
{
    const int32_t parties = MKparams->parties;
    const int32_t n = MKparams->n;

    MKTLweSample *temp = new_MKTLweSample(RLWEparams, MKparams);

    // MKTLweSample *temp1 = accum;
    MKTLweSample *temp1 = new_MKTLweSample(RLWEparams, MKparams);
    MKtLweCopy(temp1, accum, MKparams);    


    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            const int32_t baraij = bara[n*i+j];

            if (baraij == 0) continue; //indeed, this is an easy case!

            MKtfhe_MuxRotate_v2m1(temp, temp1, bk + (n*i+j), baraij, RLWEparams, MKparams, RLWEkey);

            swap(temp, temp1);
        }
    }

    if (temp1 != accum) {
        MKtLweCopy(accum, temp1, MKparams);
    }


    delete_MKTLweSample(temp1);
    delete_MKTLweSample(temp);
}





// MK Blind rotate
// Only the PK part of RLWEkey is used 
EXPORT void MKtfhe_blindRotateFFT_v2m1(MKTLweSample *accum, const MKTGswUESampleFFT_v2 *bkFFT, 
    const int32_t *bara, const TLweParams* RLWEparams, const MKTFHEParams *MKparams, 
    const MKRLweKey *MKrlwekey) 
{
    const int32_t parties = MKparams->parties;
    const int32_t n = MKparams->n;

    MKTLweSample *temp = new_MKTLweSample(RLWEparams, MKparams);

    // MKTLweSample *temp1 = accum;
    MKTLweSample *temp1 = new_MKTLweSample(RLWEparams, MKparams);
    MKtLweCopy(temp1, accum, MKparams); 



    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            const int32_t baraij = bara[n*i+j];

            if (baraij == 0) continue; //indeed, this is an easy case!

            MKtfhe_MuxRotateFFT_v2m1(temp, temp1, bkFFT + (n*i+j), baraij, RLWEparams, MKparams, MKrlwekey); 
            swap(temp, temp1);

        }
    }

    if (temp1 != accum) {
        MKtLweCopy(accum, temp1, MKparams);
    }


    delete_MKTLweSample(temp1);
    delete_MKTLweSample(temp);
}





















// MK Blind rotate and extract 
// Only the PK part of RLWEkey is used 
EXPORT void MKtfhe_blindRotateAndExtract_v2m1(MKLweSample *result,
                                       const TorusPolynomial *v,
                                       const MKTGswUESample_v2 *bk,
                                       const int32_t barb,
                                       const int32_t *bara,
                                       const TLweParams* RLWEparams, 
                                       const MKTFHEParams *MKparams,
                                       const MKRLweKey *RLWEkey) 
{
    const int32_t N = MKparams->N;
    const int32_t _2N = 2 * N;

    TorusPolynomial *testvectbis = new_TorusPolynomial(N);
    MKTLweSample *acc = new_MKTLweSample(RLWEparams, MKparams);

    if (barb !=0)
    {
        torusPolynomialMulByXai(testvectbis, _2N - barb, v);
    }
    else
    {
        torusPolynomialCopy(testvectbis, v);
    }

    MKtLweNoiselessTrivial(acc, testvectbis, MKparams);
    MKtfhe_blindRotate_v2m1(acc, bk, bara, RLWEparams, MKparams, RLWEkey);
    MKtLweExtractMKLweSample(result, acc, MKparams);


    delete_MKTLweSample(acc);
    delete_TorusPolynomial(testvectbis);
}



// MK Blind rotate and extract 
// Only the PK part of RLWEkey is used 
EXPORT void MKtfhe_blindRotateAndExtractFFT_v2m1(MKLweSample *result,
                                       const TorusPolynomial *v,
                                       const MKTGswUESampleFFT_v2 *bkFFT,
                                       const int32_t barb,
                                       const int32_t *bara,
                                       const TLweParams* RLWEparams, 
                                       const MKTFHEParams *MKparams, 
                                       const MKRLweKey *MKrlwekey) 
{
    const int32_t N = MKparams->N;
    const int32_t _2N = 2 * N;

    TorusPolynomial *testvectbis = new_TorusPolynomial(N);
    MKTLweSample *acc = new_MKTLweSample(RLWEparams, MKparams);

    if (barb !=0)
    {
        torusPolynomialMulByXai(testvectbis, _2N - barb, v);
    }
    else
    {
        torusPolynomialCopy(testvectbis, v);
    }

    MKtLweNoiselessTrivial(acc, testvectbis, MKparams);
    MKtfhe_blindRotateFFT_v2m1(acc, bkFFT, bara, RLWEparams, MKparams, MKrlwekey);
    MKtLweExtractMKLweSample(result, acc, MKparams);


    delete_MKTLweSample(acc);
    delete_TorusPolynomial(testvectbis);
}


















// MK Bootstrap without key switching 
// Only the PK part of RLWEkey is used 
EXPORT void MKtfhe_bootstrap_woKS_v2m1(MKLweSample *result, const MKLweBootstrappingKey_v2 *bk, 
        Torus32 mu, const MKLweSample *x, const TLweParams* RLWEparams, const MKTFHEParams *MKparams,
        const MKRLweKey *RLWEkey) 
{
    const int32_t parties = MKparams->parties;
    const int32_t N = MKparams->N;
    const int32_t Nx2 = 2 * N;
    const int32_t n = MKparams->n;

    TorusPolynomial *testvect = new_TorusPolynomial(N);
    int32_t *bara = new int32_t[parties*n];

    // b*2N
    int32_t barb = modSwitchFromTorus32(x->b, Nx2);
    // a*2N
    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            bara[n*i+j] = modSwitchFromTorus32(x->a[n*i+j], Nx2);
        }
    }
    

    //the initial testvec = [mu,mu,mu,...,mu]
    for (int32_t i = 0; i < N; i++) 
    {
        testvect->coefsT[i] = mu;
    }

    MKtfhe_blindRotateAndExtract_v2m1(result, testvect, bk->bk, barb, bara, RLWEparams, MKparams, RLWEkey); 


    delete[] bara;
    delete_TorusPolynomial(testvect);
}






// MK Bootstrap without key switching 
// Only the PK part of RLWEkey is used 
EXPORT void MKtfhe_bootstrap_woKSFFT_v2m1(MKLweSample *result, const MKLweBootstrappingKeyFFT_v2 *bkFFT, 
        Torus32 mu, const MKLweSample *x, const TLweParams* RLWEparams, const MKTFHEParams *MKparams, 
        const MKRLweKey *MKrlwekey) 
{
    const int32_t parties = MKparams->parties;
    const int32_t N = MKparams->N;
    const int32_t Nx2 = 2 * N;
    const int32_t n = MKparams->n;

    TorusPolynomial *testvect = new_TorusPolynomial(N);
    int32_t *bara = new int32_t[parties*n];

    // b*2N
    int32_t barb = modSwitchFromTorus32(x->b, Nx2);
    // a*2N
    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            bara[n*i+j] = modSwitchFromTorus32(x->a[n*i+j], Nx2);
        }
    }
    

    //the initial testvec = [mu,mu,mu,...,mu]
    for (int32_t i = 0; i < N; i++) 
    {
        testvect->coefsT[i] = mu;
    }

    MKtfhe_blindRotateAndExtractFFT_v2m1(result, testvect, bkFFT->bkFFT, barb, bara, RLWEparams, MKparams, MKrlwekey);

    delete[] bara;
    delete_TorusPolynomial(testvect);
}
















// MK Bootstrap
// Only the PK part of RLWEkey is used 
EXPORT void MKtfhe_bootstrap_v2m1(MKLweSample *result, const MKLweBootstrappingKey_v2 *bk, Torus32 mu, 
        const MKLweSample *x, const LweParams* LWEparams, const LweParams* extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *RLWEkey) 
{
    MKLweSample *u = new_MKLweSample(extractedLWEparams, MKparams);

    MKtfhe_bootstrap_woKS_v2m1(u, bk, mu, x, RLWEparams, MKparams, RLWEkey);
    
    // MK Key Switching
    //MKlweKeySwitch(result, bk->ks, u, MKparams);
    MKlweKeySwitch(result, bk->ks, u, LWEparams, MKparams);


    delete_MKLweSample(u);
}




// MK Bootstrap
// Only the PK part of RLWEkey is used 
EXPORT void MKtfhe_bootstrapFFT_v2m1(MKLweSample *result, const MKLweBootstrappingKeyFFT_v2 *bkFFT, Torus32 mu, 
        const MKLweSample *x, const LweParams* LWEparams, const LweParams* extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey) 
{
    MKLweSample *u = new_MKLweSample(extractedLWEparams, MKparams);

    MKtfhe_bootstrap_woKSFFT_v2m1(u, bkFFT, mu, x, RLWEparams, MKparams, MKrlwekey);
    // MK Key Switching
    //MKlweKeySwitch(result, bkFFT->ks, u, MKparams);
    MKlweKeySwitch(result, bkFFT->ks, u, LWEparams, MKparams);


    delete_MKLweSample(u);
}











// MK Bootstrapped NAND 
// Only the PK part of RLWEkey is used 
EXPORT void MKbootsNAND_v2m1(MKLweSample *result, const MKLweSample *ca, const MKLweSample *cb, 
        const MKLweBootstrappingKey_v2 *bk, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *RLWEkey) 
{
    static const Torus32 MU = modSwitchToTorus32(1, 8);

    MKLweSample *temp_result = new_MKLweSample(LWEparams, MKparams);

    //compute: (0,1/8) - ca - cb
    static const Torus32 NandConst = modSwitchToTorus32(1, 8);
    MKlweNoiselessTrivial(temp_result, NandConst, MKparams);
    MKlweSubTo(temp_result, ca, MKparams);
    MKlweSubTo(temp_result, cb, MKparams);


    //if the phase is positive, the result is 1/8
    //if the phase is positive, else the result is -1/8

    MKtfhe_bootstrap_v2m1(result, bk, MU, temp_result, LWEparams, extractedLWEparams, RLWEparams, MKparams, RLWEkey);

    delete_MKLweSample(temp_result);
}




// MK Bootstrapped NAND 
// Only the PK part of RLWEkey is used 
EXPORT void MKbootsNAND_FFT_v2m1(MKLweSample *result, const MKLweSample *ca, const MKLweSample *cb, 
        const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey) 
{
    static const Torus32 MU = modSwitchToTorus32(1, 8);

    MKLweSample *temp_result = new_MKLweSample(LWEparams, MKparams);

    //compute: (0,1/8) - ca - cb
    static const Torus32 NandConst = modSwitchToTorus32(1, 8);
    MKlweNoiselessTrivial(temp_result, NandConst, MKparams);
    MKlweSubTo(temp_result, ca, MKparams);
    MKlweSubTo(temp_result, cb, MKparams);


    //if the phase is positive, the result is 1/8
    //if the phase is positive, else the result is -1/8
    MKtfhe_bootstrapFFT_v2m1(result, bkFFT, MU, temp_result, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);   

    delete_MKLweSample(temp_result);
}


























































/* ********************************************************************************
************************** Bootstrapping method 2 *********************************
******************************************************************************** */


// MUX -> rotate
// Only the PK part of RLWEkey is used 
void MKtfhe_MuxRotate_v2m2(MKTLweSample *result, MKTLweSample *accum, const MKTGswUESample_v2* bki, 
    const int32_t barai, const TLweParams* RLWEparams, const MKTFHEParams* MKparams, const MKRLweKey *RLWEkey) 
{
    MKTLweSample *temp_result = new_MKTLweSample(RLWEparams, MKparams);

    // ACC = BKi*[(X^barai-1)*ACC]+ACC
    // temp = (X^barai-1)*ACC
    MKtLweMulByXaiMinusOne(temp_result, barai, accum, MKparams);
    // temp *= BKi
    MKtGswUEExternMulToMKtLwe_v2m2(result, temp_result, bki, RLWEparams, MKparams, RLWEkey);
    // ACC += temp
    MKtLweAddTo(result, accum, MKparams);

    delete_MKTLweSample(temp_result);
}




// MUX -> rotate
// Only the PK part of RLWEkey is used 
void MKtfhe_MuxRotateFFT_v2m2(MKTLweSample *result, MKTLweSample *accum, const MKTGswUESampleFFT_v2 *bkiFFT, 
    const int32_t barai, const TLweParams* RLWEparams, const MKTFHEParams* MKparams, const MKRLweKey *RLWEkey) 
{
    MKTLweSample *temp_result = new_MKTLweSample(RLWEparams, MKparams);

    // ACC = BKi*[(X^barai-1)*ACC]+ACC
    // temp = (X^barai-1)*ACC
    MKtLweMulByXaiMinusOne(temp_result, barai, accum, MKparams);
    // temp *= BKi
    MKtGswUEExternMulToMKtLwe_FFT_v2m2(result, temp_result,bkiFFT, RLWEparams, MKparams, RLWEkey);
    // ACC += temp
    MKtLweAddTo(result, accum, MKparams);

    delete_MKTLweSample(temp_result);
}










// MK Blind rotate
// Only the PK part of RLWEkey is used 
EXPORT void MKtfhe_blindRotate_v2m2(MKTLweSample *accum, const MKTGswUESample_v2 *bk, const int32_t *bara, 
    const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *RLWEkey) 
{
    const int32_t parties = MKparams->parties;
    const int32_t n = MKparams->n;

    MKTLweSample *temp = new_MKTLweSample(RLWEparams, MKparams);

    // MKTLweSample *temp1 = accum;
    MKTLweSample *temp1 = new_MKTLweSample(RLWEparams, MKparams);
    MKtLweCopy(temp1, accum, MKparams);    


    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            const int32_t baraij = bara[n*i+j];

            if (baraij == 0) continue; //indeed, this is an easy case!

            MKtfhe_MuxRotate_v2m2(temp, temp1, bk + (n*i+j), baraij, RLWEparams, MKparams, RLWEkey);

            swap(temp, temp1);
        }
    }

    if (temp1 != accum) {
        MKtLweCopy(accum, temp1, MKparams);
    }


    delete_MKTLweSample(temp1);
    delete_MKTLweSample(temp);
}






// MK Blind rotate
// Only the PK part of RLWEkey is used 
EXPORT void MKtfhe_blindRotateFFT_v2m2(MKTLweSample *accum, const MKTGswUESampleFFT_v2 *bkFFT, 
    const int32_t *bara, const TLweParams* RLWEparams, const MKTFHEParams *MKparams, 
    const MKRLweKey *MKrlwekey) 
{
    const int32_t parties = MKparams->parties;
    const int32_t n = MKparams->n;

    MKTLweSample *temp = new_MKTLweSample(RLWEparams, MKparams);

    // MKTLweSample *temp1 = accum;
    MKTLweSample *temp1 = new_MKTLweSample(RLWEparams, MKparams);
    MKtLweCopy(temp1, accum, MKparams); 



    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            const int32_t baraij = bara[n*i+j];

            if (baraij == 0) continue; //indeed, this is an easy case!

            MKtfhe_MuxRotateFFT_v2m2(temp, temp1, bkFFT + (n*i+j), baraij, RLWEparams, MKparams, MKrlwekey); 
            swap(temp, temp1);

        }
    }

    if (temp1 != accum) {
        MKtLweCopy(accum, temp1, MKparams);
    }


    delete_MKTLweSample(temp1);
    delete_MKTLweSample(temp);
}


























// MK Blind rotate and extract 
// Only the PK part of RLWEkey is used 
EXPORT void MKtfhe_blindRotateAndExtract_v2m2(MKLweSample *result,
                                       const TorusPolynomial *v,
                                       const MKTGswUESample_v2 *bk,
                                       const int32_t barb,
                                       const int32_t *bara,
                                       const TLweParams* RLWEparams, 
                                       const MKTFHEParams *MKparams,
                                       const MKRLweKey *RLWEkey) 
{
    const int32_t N = MKparams->N;
    const int32_t _2N = 2 * N;

    TorusPolynomial *testvectbis = new_TorusPolynomial(N);
    MKTLweSample *acc = new_MKTLweSample(RLWEparams, MKparams);

    if (barb !=0)
    {
        torusPolynomialMulByXai(testvectbis, _2N - barb, v);
    }
    else
    {
        torusPolynomialCopy(testvectbis, v);
    }

    MKtLweNoiselessTrivial(acc, testvectbis, MKparams);
    MKtfhe_blindRotate_v2m2(acc, bk, bara, RLWEparams, MKparams, RLWEkey);
    MKtLweExtractMKLweSample(result, acc, MKparams);


    delete_MKTLweSample(acc);
    delete_TorusPolynomial(testvectbis);
}




// MK Blind rotate and extract 
// Only the PK part of RLWEkey is used 
EXPORT void MKtfhe_blindRotateAndExtractFFT_v2m2(MKLweSample *result,
                                       const TorusPolynomial *v,
                                       const MKTGswUESampleFFT_v2 *bkFFT,
                                       const int32_t barb,
                                       const int32_t *bara,
                                       const TLweParams* RLWEparams, 
                                       const MKTFHEParams *MKparams, 
                                       const MKRLweKey *MKrlwekey) 
{
    const int32_t N = MKparams->N;
    const int32_t _2N = 2 * N;

    TorusPolynomial *testvectbis = new_TorusPolynomial(N);
    MKTLweSample *acc = new_MKTLweSample(RLWEparams, MKparams);

    if (barb !=0)
    {
        torusPolynomialMulByXai(testvectbis, _2N - barb, v);
    }
    else
    {
        torusPolynomialCopy(testvectbis, v);
    }

    MKtLweNoiselessTrivial(acc, testvectbis, MKparams);
    MKtfhe_blindRotateFFT_v2m2(acc, bkFFT, bara, RLWEparams, MKparams, MKrlwekey);
    MKtLweExtractMKLweSample(result, acc, MKparams);


    delete_MKTLweSample(acc);
    delete_TorusPolynomial(testvectbis);
}



















// MK Bootstrap without key switching 
// Only the PK part of RLWEkey is used 
EXPORT void MKtfhe_bootstrap_woKS_v2m2(MKLweSample *result, const MKLweBootstrappingKey_v2 *bk, 
        Torus32 mu, const MKLweSample *x, const TLweParams* RLWEparams, const MKTFHEParams *MKparams,
        const MKRLweKey *RLWEkey) 
{
    const int32_t parties = MKparams->parties;
    const int32_t N = MKparams->N;
    const int32_t Nx2 = 2 * N;
    const int32_t n = MKparams->n;

    TorusPolynomial *testvect = new_TorusPolynomial(N);
    int32_t *bara = new int32_t[parties*n];

    // b*2N
    int32_t barb = modSwitchFromTorus32(x->b, Nx2);
    // a*2N
    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            bara[n*i+j] = modSwitchFromTorus32(x->a[n*i+j], Nx2);
        }
    }
    

    //the initial testvec = [mu,mu,mu,...,mu]
    for (int32_t i = 0; i < N; i++) 
    {
        testvect->coefsT[i] = mu;
    }

    MKtfhe_blindRotateAndExtract_v2m2(result, testvect, bk->bk, barb, bara, RLWEparams, MKparams, RLWEkey); 


    delete[] bara;
    delete_TorusPolynomial(testvect);
}





// MK Bootstrap without key switching 
// Only the PK part of RLWEkey is used 
EXPORT void MKtfhe_bootstrap_woKSFFT_v2m2(MKLweSample *result, const MKLweBootstrappingKeyFFT_v2 *bkFFT, 
        Torus32 mu, const MKLweSample *x, const TLweParams* RLWEparams, const MKTFHEParams *MKparams, 
        const MKRLweKey *MKrlwekey) 
{
    const int32_t parties = MKparams->parties;
    const int32_t N = MKparams->N;
    const int32_t Nx2 = 2 * N;
    const int32_t n = MKparams->n;

    TorusPolynomial *testvect = new_TorusPolynomial(N);
    int32_t *bara = new int32_t[parties*n];

    // b*2N
    int32_t barb = modSwitchFromTorus32(x->b, Nx2);
    // a*2N
    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            bara[n*i+j] = modSwitchFromTorus32(x->a[n*i+j], Nx2);
        }
    }
    

    //the initial testvec = [mu,mu,mu,...,mu]
    for (int32_t i = 0; i < N; i++) 
    {
        testvect->coefsT[i] = mu;
    }

    MKtfhe_blindRotateAndExtractFFT_v2m2(result, testvect, bkFFT->bkFFT, barb, bara, RLWEparams, MKparams, MKrlwekey);

    delete[] bara;
    delete_TorusPolynomial(testvect);
}













// MK Bootstrap
// Only the PK part of RLWEkey is used 
EXPORT void MKtfhe_bootstrap_v2m2(MKLweSample *result, const MKLweBootstrappingKey_v2 *bk, Torus32 mu, 
        const MKLweSample *x, const LweParams* LWEparams, const LweParams* extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *RLWEkey) 
{
    MKLweSample *u = new_MKLweSample(extractedLWEparams, MKparams);

    MKtfhe_bootstrap_woKS_v2m2(u, bk, mu, x, RLWEparams, MKparams, RLWEkey);
    
    // MK Key Switching
    //MKlweKeySwitch(result, bk->ks, u, MKparams);
    MKlweKeySwitch(result, bk->ks, u, LWEparams, MKparams);


    delete_MKLweSample(u);
}




// MK Bootstrap
// Only the PK part of RLWEkey is used 
EXPORT void MKtfhe_bootstrapFFT_v2m2(MKLweSample *result, const MKLweBootstrappingKeyFFT_v2 *bkFFT, Torus32 mu, 
        const MKLweSample *x, const LweParams* LWEparams, const LweParams* extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey) 
{
    MKLweSample *u = new_MKLweSample(extractedLWEparams, MKparams);

    MKtfhe_bootstrap_woKSFFT_v2m2(u, bkFFT, mu, x, RLWEparams, MKparams, MKrlwekey);
    // MK Key Switching
    //MKlweKeySwitch(result, bkFFT->ks, u, MKparams);
    MKlweKeySwitch(result, bkFFT->ks, u, LWEparams, MKparams);


    delete_MKLweSample(u);
}













// MK Bootstrapped NAND 
// Only the PK part of RLWEkey is used 
EXPORT void MKbootsNAND_v2m2(MKLweSample *result, const MKLweSample *ca, const MKLweSample *cb, 
        const MKLweBootstrappingKey_v2 *bk, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *RLWEkey) 
{
    static const Torus32 MU = modSwitchToTorus32(1, 8);

    MKLweSample *temp_result = new_MKLweSample(LWEparams, MKparams);

    //compute: (0,1/8) - ca - cb
    static const Torus32 NandConst = modSwitchToTorus32(1, 8);
    MKlweNoiselessTrivial(temp_result, NandConst, MKparams);
    MKlweSubTo(temp_result, ca, MKparams);
    MKlweSubTo(temp_result, cb, MKparams);


    //if the phase is positive, the result is 1/8
    //if the phase is positive, else the result is -1/8

    MKtfhe_bootstrap_v2m2(result, bk, MU, temp_result, LWEparams, extractedLWEparams, RLWEparams, MKparams, RLWEkey);

    delete_MKLweSample(temp_result);
}





// MK Bootstrapped NAND 
// Only the PK part of RLWEkey is used 
EXPORT void MKbootsNAND_FFT_v2m2(MKLweSample *result, const MKLweSample *ca, const MKLweSample *cb, 
        const MKLweBootstrappingKeyFFT_v2 *bkFFT, const LweParams* LWEparams, const LweParams *extractedLWEparams, 
        const TLweParams* RLWEparams, const MKTFHEParams *MKparams, const MKRLweKey *MKrlwekey) 
{
    static const Torus32 MU = modSwitchToTorus32(1, 8);

    MKLweSample *temp_result = new_MKLweSample(LWEparams, MKparams);

    //compute: (0,1/8) - ca - cb
    static const Torus32 NandConst = modSwitchToTorus32(1, 8);
    MKlweNoiselessTrivial(temp_result, NandConst, MKparams);
    MKlweSubTo(temp_result, ca, MKparams);
    MKlweSubTo(temp_result, cb, MKparams);


    //if the phase is positive, the result is 1/8
    //if the phase is positive, else the result is -1/8
    MKtfhe_bootstrapFFT_v2m2(result, bkFFT, MU, temp_result, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);   

    delete_MKLweSample(temp_result);
}



#include <cstdlib>
#include <iostream>
#include <random>
#include <cassert>
#include "tlwe_functions.h"
#include "numeric_functions.h"
#include "polynomials_arithmetic.h"
#include "lagrangehalfc_arithmetic.h"
#include "tfhe_generic_templates.h"
#include "lwe-functions.h"
#include "lweparams.h"
#include "lwekeyswitch.h"

#include "mkTFHEparams.h"
#include "mkTFHEsamples.h"
#include "mkTFHEkeys.h"
#include "mkTFHEfunctions.h"

using namespace std;


// MKLwe 
// key generation for every party 
EXPORT void MKLweKeyGen(MKLweKey* result) {

    const int32_t parties = result->MKparams->parties;

    for (int i = 0; i < parties; ++i)
    {
        lweKeyGen(&result->key[i]);
    }
}





// MKRLwe
// key generation for every party
// secret and public keys
EXPORT void MKRLweKeyGen(MKRLweKey *result) {

    const int32_t parties = result->MKparams->parties;
    const int32_t dg = result->MKparams->dg;
    const int32_t N = result->MKparams->N;
    const double stdevRLWEkey = result->MKparams->stdevRLWEkey; 
    // secret keys
    for (int i = 0; i < parties; ++i)
    {
        tLweKeyGen(&result->key[i]); 
    }

    // public keys
    // a
    for (int j = 0; j < dg; ++j)
    {
        torusPolynomialUniform(&result->Pkey[parties*dg + j]);
    }
    // b_i = +a*s_i + e_i
    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < dg; ++j)
        {
            // b_i = e_i 
            for (int l = 0; l < N; ++l)
            {
                result->Pkey[i*dg + j].coefsT[l] = gaussian32(0, stdevRLWEkey);
            }      
            // b_i = e_i + a*s_i
            torusPolynomialAddMulRFFT1(&result->Pkey[i*dg + j], result->key[i].key, &result->Pkey[parties*dg + j]); 
        }
    }
    
}



//extractions Ring Lwe -> Lwe (extracted)
EXPORT void MKtLweExtractKey(MKLweKey* LWEkey, const MKRLweKey* RLWEkey) 
{
    const int32_t parties = RLWEkey->MKparams->parties;
    const int32_t N = RLWEkey->MKparams->N;
    
    assert(result->params->n_extract == N);

    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            LWEkey->key[i].key[j]=RLWEkey->key[i].key->coefs[j];
        }
    }
}
































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
        const MKTFHEParams* MKparams)
{
    const int32_t n_in = result->n_in;
    const int32_t parties = result->parties;
    const int32_t Bksbit = result->Bksbit;
    const int32_t Bks = result->Bks;
    const int32_t dks = result->dks;
    const double stdevKS = MKparams->stdevKS;
    const int32_t sizeks = parties*n_in*dks*(Bks-1);
    

    double err = 0;

    // chose a random vector of gaussian noises
    double* noise = new double[sizeks];
    for (int32_t i = 0; i < sizeks; ++i){
        normal_distribution<double> distribution(0.0,stdevKS); 
        noise[i] = distribution(generator);
        err += noise[i];
    }
    // recenter the noises
    err = err/sizeks;
    for (int32_t i = 0; i < sizeks; ++i) noise[i] -= err;




    // generate the ks
    int32_t index = 0; 
    for (int p = 0; p < parties; ++p)
    {
        for (int32_t i = 0; i < n_in; ++i) 
        {
            for (int32_t j = 0; j < dks; ++j) 
            {
                // term h=0 as trivial encryption of 0 (it will not be used in the KeySwitching)
                MKlweNoiselessTrivial(&result->ks[p][i][j][0], 0, MKparams);

                for (int32_t h = 1; h < Bks; ++h) { // not the 0 term
                    
                    Torus32 mess = (in_key->key[p].key[i]*h)*(1<<(32-(j+1)*Bksbit));

                    MKlweSymEncryptWithExternalNoise(&result->ks[p][i][j][h], mess, noise[index], stdevKS, out_key);
                    index += 1;
                }
            }
        }
    }

    delete[] noise; 
}
























/* *******************************************************
*************** Bootstrapping Key v2 *********************
******************************************************* */


EXPORT void init_MKLweBootstrappingKey_v2(MKLweBootstrappingKey_v2 *obj,
        const LweParams* LWEparams, const TLweParams* RLWEparams, const MKTFHEParams* MKparams) 
{
    const int32_t n = MKparams->n;
    const int32_t parties = MKparams->parties;
    const int32_t n_extract = MKparams->n_extract;
    const int32_t dks = MKparams->dks;
    const int32_t Bksbit = MKparams->Bksbit;

    MKTGswUESample_v2* bk = new_MKTGswUESample_v2_array(n*parties, RLWEparams, MKparams);    
    //MKLweKeySwitchKey *ks = new_MKLweKeySwitchKey(n_in, LWEparams, MKparams);
    LweKeySwitchKey *ks = new_LweKeySwitchKey_array(parties, n_extract, dks, Bksbit, LWEparams);

    new(obj) MKLweBootstrappingKey_v2(MKparams, bk, ks);
}

EXPORT void destroy_MKLweBootstrappingKey_v2(MKLweBootstrappingKey_v2 *obj) {
    delete_LweKeySwitchKey_array(obj->MKparams->parties, obj->ks);
    //delete_MKLweKeySwitchKey(obj->ks);
    delete_MKTGswUESample_v2_array(obj->MKparams->parties*obj->MKparams->n, obj->bk);
    obj->~MKLweBootstrappingKey_v2();
}







EXPORT void MKlweCreateBootstrappingKey_v2(MKLweBootstrappingKey_v2* result, const MKLweKey* LWEkey, 
        const MKRLweKey* RLWEkey, const MKLweKey* extractedLWEkey, const LweParams *extractedLWEparams,
        const LweParams *LWEparams, const TLweParams *RLWEparams, const MKTFHEParams* MKparams)
{
    const int32_t parties = MKparams->parties;
    const int32_t n = MKparams->n;

    
    // bootstrapping key generation -> expansion. 
    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            // fix the party 
            result->bk[i*n+j].party = i; 
            // LWE key si encrypted with RLWE key Si
            MKTGswUniEncryptI_v2(&result->bk[i*n+j], LWEkey->key[i].key[j], i, MKparams->stdevBK, RLWEkey); // party = i
        }
    }
    

    // key switching key
    //MKlweCreateKeySwitchKey(result->ks, extractedLWEkey, LWEkey, MKparams);
    for (int i = 0; i < parties; ++i)
    {
        // every party generates his KS key independently 
        lweCreateKeySwitchKey(&result->ks[i], &extractedLWEkey->key[i], &LWEkey->key[i]);
    }
    

    result->MKparams = MKparams;
}










// FFT
EXPORT void init_MKLweBootstrappingKeyFFT_v2(MKLweBootstrappingKeyFFT_v2 *obj, 
    const MKLweBootstrappingKey_v2 *bk, const LweParams* LWEparams, const TLweParams* RLWEparams, const MKTFHEParams* MKparams) 
{
    const int32_t n = MKparams->n;
    const int32_t dks = MKparams->dks;
    const int32_t dg = MKparams->dg;
    const int32_t Bksbit = MKparams->Bksbit;
    const int32_t Bks = 1 << Bksbit;
    const int32_t n_extract = MKparams->n_extract;
    const int32_t parties = MKparams->parties;
    

    //MKLweKeySwitchKey *ks = new_MKLweKeySwitchKey(n_extract, LWEparams, MKparams);
    LweKeySwitchKey *ks = new_LweKeySwitchKey_array(parties, n_extract, dks, Bksbit, LWEparams);
    // Copy the KeySwitching key
    for (int p = 0; p < parties; ++p)
    {
        for (int32_t i = 0; i < n_extract; i++) 
        {
            for (int32_t j = 0; j < dks; j++) 
            {
                for (int32_t l = 0; l < Bks; l++) {
                    lweCopy(&ks[p].ks[i][j][l], &bk->ks[p].ks[i][j][l], LWEparams);
                }
            }
        }
    }

    
    // Bootstrapping Key FFT 
    int32_t nb_polys = 3*dg;
    MKTGswUESampleFFT_v2 *bkFFT = new_MKTGswUESampleFFT_v2_array(n*parties, RLWEparams, MKparams, 0, 0.0);
    // convert bk to bkFFT
    clock_t begin = clock();
    for (int p = 0; p < parties; ++p)
    {
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < nb_polys; ++j)
            {
                TorusPolynomial_ifft(&bkFFT[p*n+i].d[j], &bk->bk[p*n+i].d[j]);
            }
            bkFFT[p*n+i].party = bk->bk[p*n+i].party; 
        }
    }
    clock_t end = clock();
    double time = ((double) end - begin)/CLOCKS_PER_SEC;
    cout << "Time BK FFT conversion: " << time << " seconds" << endl;

    
    new(obj) MKLweBootstrappingKeyFFT_v2(MKparams, bkFFT, ks);
}



//destroys the MKLweBootstrappingKeyFFT_v2 structure
EXPORT void destroy_MKLweBootstrappingKeyFFT_v2(MKLweBootstrappingKeyFFT_v2 *obj) {
    delete_LweKeySwitchKey_array(obj->MKparams->parties, obj->ks);
    //delete_MKLweKeySwitchKey((MKLweKeySwitchKey *) obj->ks);
    delete_MKTGswUESampleFFT_v2_array(obj->MKparams->n*obj->MKparams->parties, obj->bkFFT);
    obj->~MKLweBootstrappingKeyFFT_v2();
}




#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include "tfhe.h"
#include "polynomials.h"
#include "lwesamples.h"
#include "lwekey.h"
#include "lweparams.h"
#include "tlwe.h"
#include "tgsw.h"



#include "mkTFHEparams.h"
#include "mkTFHEkeys.h"
#include "mkTFHEkeygen.h"
#include "mkTFHEsamples.h"
#include "mkTFHEfunctions.h"





 

using namespace std;



// **********************************************************************************
// ********************************* MAIN *******************************************
// **********************************************************************************


void dieDramatically(string message) {
    cerr << message << endl;
    abort();
}


void print_error(TorusPolynomial &poly, vector<int32_t> &plain, int32_t Bgscale){
    double maxerr = 0;
    int N = plain.size();  
    for (int h = 0; h < N; ++h)
    {
        double a = (double) poly.coefsT[h] / pow(2,32);
        double b = (double) plain[h] / (double) Bgscale;
        maxerr = max(maxerr, abs(a-b)); 
    }
    cout << maxerr << endl; 
}




//EXPORT void tLweExtractKey(LweKey* result, const TLweKey* key); //TODO: change the name and put in a .h
//EXPORT void tfhe_createLweBootstrappingKeyFFT(LweBootstrappingKeyFFT* bk, const LweKey* key_in, const TGswKey* rgsw_key);
//EXPORT void tfhe_bootstrapFFT(LweSample* result, const LweBootstrappingKeyFFT* bk, Torus32 mu1, Torus32 mu0, const LweSample* x);




int32_t main(int32_t argc, char **argv) {

    // Test trials
    const int32_t nb_trials = 1;

    // generate params 
    static const int32_t k = 1;
    //static const int32_t bk_l = 2;
    //static const int32_t bk_Bgbit = 10;
    //static const int32_t ks_basebit = 2;
    //static const int32_t ks_length = 8;
    static const double ks_stdev = 2.4e-5; // 2.44e-5; //standard deviation
    static const double bk_stdev = 7.2e-9; // 7.18e-9; //standard deviation
    static const double max_stdev = 0.012467; //max standard deviation for a 1/4 msg space
    // new params
    static const int32_t n = 500;            // LWE modulus
    static const int32_t n_extract = 1024;    // LWE extract modulus (used in bootstrapping)
    static const int32_t hLWE = 0;         // HW secret key LWE --> not used
    static const double stdevLWE = 2.4e-5; // 0.012467;      // LWE ciphertexts standard deviation
    //static const double stdevLWE = 0; // 0.012467;      // LWE ciphertexts standard deviation
    static const int32_t Bksbit = 2;       // Base bit key switching
    static const int32_t dks = 6;          // dimension key switching
    static const double stdevKS = 2.4e-5; // 2.44e-5;       // KS key standard deviation
    //static const double stdevKS = 0; // 2.44e-5;       // KS key standard deviation
    static const int32_t N = 1024;            // RLWE,RGSW modulus
    static const int32_t hRLWE = 0;        // HW secret key RLWE,RGSW --> not used
    static const double stdevRLWEkey = 7.2e-9; // 7.18e-11; //7.18e-10; // 0.012467;  // RLWE key standard deviation
    static const double stdevRLWE = 7.2e-9; // 0.012467;     // RLWE ciphertexts standard deviation
    static const double stdevRGSW = 7.2e-9; // 7.18e-11; //7.18e-10;     // RGSW ciphertexts standard deviation 
    //static const double stdevRLWEkey = 0; // 7.18e-11; //7.18e-10; // 0.012467;  // RLWE key standard deviation
    //static const double stdevRLWE = 0; // 0.012467;     // RLWE ciphertexts standard deviation
    //static const double stdevRGSW = 0; //
    static const int32_t Bgbit = 5; // 10;        // Base bit gadget
    static const int32_t dg = 5; // 2;           // dimension gadget
    //static const double stdevBK = 7.2e-9; // 7.18e-9;       // BK standard deviation
    static const double stdevBK = 7.18e-9;       // BK standard deviation
    static const int32_t parties = 2;      // number of parties


    // params
    LweParams *extractedLWEparams = new_LweParams(n_extract, ks_stdev, max_stdev);
    LweParams *LWEparams = new_LweParams(n, ks_stdev, max_stdev);
    TLweParams *RLWEparams = new_TLweParams(N, k, bk_stdev, max_stdev);
    // TGswParams *TGSWparams = new_TGswParams(dg, Bgbit, RLWEparams);
    MKTFHEParams *MKparams = new_MKTFHEParams(n, n_extract, hLWE, stdevLWE, Bksbit, dks, stdevKS, N, 
                            hRLWE, stdevRLWEkey, stdevRLWE, stdevRGSW, Bgbit, dg, stdevBK, parties);


    cout << "Params: DONE!" << endl;






    // Key generation 
    // LWE key
    cout << "Starting KEY GENERATION" << endl;
    clock_t begin_KG = clock();
    // RLWE key 
    MKRLweKey* MKrlwekey = new_MKRLweKey(RLWEparams, MKparams);
    MKRLweKeyGen(MKrlwekey);

    // verify MKrlwekey
    
    // for (int i = 0; i <= parties; ++i)
    // {
    //     cout << "key party: " << i << " - ";
    //     for (int j = 0; j < dg; ++j)
    //     {
    //         for (int h = 0; h < N; ++h)
    //         {
    //             cout << MKrlwekey->Pkey[dg*i + j].coefsT[h] << ", ";
    //         }
    //         cout << " - ";
    //     }
    //     cout << endl;        
    // }
    
    TorusPolynomial* phase_rlwekey = new_TorusPolynomial(N);
    vector<int32_t> zeroplain(N);
    for(int i = 0; i < N; i++){
        zeroplain[i] = 0;
    }
    for (int i = 0; i < parties; ++i)
    {
        for (int j = 0; j < dg; ++j)
        {
            torusPolynomialCopy(phase_rlwekey, &MKrlwekey->Pkey[i*dg + j]); // phi = b
            torusPolynomialSubMulR(phase_rlwekey, MKrlwekey->key[i].key, &MKrlwekey->Pkey[parties*dg + j]); // phi = b - a*s_party
            cout << "phase of pub key"; 
            print_error(*phase_rlwekey, zeroplain, 1); 
            // for (int h = 0; h < N; ++h)
            // {
            //     cout << phase_rlwekey->coefsT[h] << ", ";
            //     // cout << modSwitchFromTorus32(phase_rlwekey->coefsT[h], 2) << ", ";
            // }
            // cout << endl;
        }
        cout << endl;
    }
    delete_TorusPolynomial(phase_rlwekey);
    

    cout << "KeyGen MKrlwekey: DONE!" << endl;
    clock_t end_KG = clock();
    double time_KG = ((double) end_KG - begin_KG)/CLOCKS_PER_SEC;
    cout << "Finished KEY GENERATION" << endl;
    










    int32_t Msize = 2; // 7;
    static uniform_int_distribution<int32_t> distributionRLWE(0, Msize - 1);
    static uniform_int_distribution<int32_t> distributionRGSW(0,1);
  











    for (int trial = 0; trial < nb_trials; ++trial)
    {
        cout << "****************" << endl;
        cout << "Trial: " << trial << endl;
        cout << "****************" << endl;


        // sample RLWE
        TorusPolynomial *sample_clear = new_TorusPolynomial(N);
        for (int32_t i = 0; i < N; ++i) {
            int32_t temp = distributionRLWE(generator);
            sample_clear->coefsT[i] = modSwitchToTorus32(temp, Msize);
        }
        MKTLweSample* sample_enc = new_MKTLweSample(RLWEparams, MKparams);
        MKtLweSymEncrypt(sample_enc, sample_clear, stdevRLWE, MKrlwekey);

        // sample TGSW_UE -> just one integer
        int32_t sampleUE_clear = distributionRGSW(generator);
        MKTGswUESample_v2* sampleUE_enc = new_MKTGswUESample_v2(RLWEparams, MKparams);
        MKTGswUniEncryptI_v2(sampleUE_enc, sampleUE_clear, 0, stdevRGSW, MKrlwekey); // party = 0

        // Convert the UE sample to UE_FFT
        LagrangeHalfCPolynomial *arrUE = new_LagrangeHalfCPolynomial_array(3*dg, N);
        for (int i = 0; i < 3*dg; ++i)
        {
            TorusPolynomial_ifft(&arrUE[i], &sampleUE_enc->d[i]);
        }
        MKTGswUESampleFFT_v2* sampleUE_FFT_enc = new_MKTGswUESampleFFT_v2(RLWEparams, MKparams, arrUE, sampleUE_enc->party, sampleUE_enc->current_variance);
        sampleUE_FFT_enc->party = sampleUE_enc->party;

        // result RLWE
        TorusPolynomial *result_clear = new_TorusPolynomial(N);
        


        




        // verify Encrypt
        // sample RLWE
        TorusPolynomial *sample_decrypt = new_TorusPolynomial(N);
        MKtLweSymDecrypt(sample_decrypt, sample_enc, MKrlwekey, Msize);
        cout << "sample decrypt, errors?" << endl;
        for (int i = 0; i < N; ++i)
        {
            if (modSwitchFromTorus32(sample_decrypt->coefsT[i], Msize) != modSwitchFromTorus32(sample_clear->coefsT[i], Msize))
            {
                cout << modSwitchFromTorus32(sample_decrypt->coefsT[i], Msize) << " -- " << modSwitchFromTorus32(sample_clear->coefsT[i], Msize) << endl;
            }
        }



        // sample RGSW UE
        TorusPolynomial *sampleUE_decryptF = new_TorusPolynomial_array(dg, N);
        IntPolynomial* r = new_IntPolynomial(N);
        MKtGswSymDecrypt_v2(sampleUE_decryptF, sampleUE_enc, MKrlwekey); // ~r*g[j]

        // VERIFY F
        cout << "r = ";
        for (int h = 0; h < N; ++h)
        {
            Torus32 temp = sampleUE_decryptF[0].coefsT[h];
            double tempd = (double) temp / pow(2,32 - MKparams->Bgbit);
            if (tempd > 0.5)
            {
                r->coefs[h] = 1;
                cout << 1 << " ";
            }
            else{
                r->coefs[h] = 0;
                cout << 0 << " ";
            }
        }        
        cout << endl;    
        for (int j = 1; j < dg; ++j)
        {
            cout << "r = ";
            for (int h = 0; h < N; ++h)
            {
                Torus32 temp = sampleUE_decryptF[j].coefsT[h];
                double tempd = (double) temp / pow(2,32 - (j + 1) * MKparams->Bgbit);
                if (tempd > 0.5)
                {
                    cout << 1 << " ";
                }
                else{
                    cout << 0 << " ";
                }
            }        
            cout << endl;       
        }

        // VERIFY d
        TorusPolynomial *sampleUE_decryptD = new_TorusPolynomial_array(dg, N);
        for (int j = 0; j < dg; ++j)
        {
            torusPolynomialCopy(&sampleUE_decryptD[j], &sampleUE_enc->d[j]); // phi = d[j] 
            // phi = d - r*a
            torusPolynomialSubMulR(&sampleUE_decryptD[j], r, &MKrlwekey->Pkey[dg*parties + j]);
        }

        for (int j = 0; j < dg; ++j)
        {
            cout << "clear = " << sampleUE_clear << endl;
            cout << "mu = ";
            for (int h = 0; h < N; ++h)
            {
                Torus32 temp = sampleUE_decryptD[j].coefsT[h];
                double tempd = (double) temp / pow(2,32 - (j + 1) * MKparams->Bgbit);
                if (tempd > 0.5)
                {
                    cout << 1 << " ";
                }
                else{
                    cout << 0 << " ";
                }
            }        
            cout << endl;       
        }


        delete_IntPolynomial(r);
        delete_TorusPolynomial_array(dg, sampleUE_decryptD);
        delete_TorusPolynomial_array(dg, sampleUE_decryptF);


















        // Multiplication in clear
        for (int i = 0; i < N; ++i)
        {
            result_clear->coefsT[i] = sampleUE_clear*sample_clear->coefsT[i];
        }







        /* ********************************************************************************
        *********************** EXTERNAL PRODUCT method 1 *********************************
        ******************************************************************************** */



        MKTLweSample* result_enc_m1 = new_MKTLweSample(RLWEparams, MKparams);
        // External product 
        clock_t begin_EP_m1 = clock();
        MKtGswUEExternMulToMKtLwe_v2m1(result_enc_m1, sample_enc, sampleUE_enc, RLWEparams, MKparams, MKrlwekey);
        clock_t end_EP_m1 = clock();
        cout << "Time per MK_extProd_v2_m1 (microsecs)... " << (end_EP_m1 - begin_EP_m1) << endl;



        // Verify 
        TorusPolynomial *result_decrypt_m1 = new_TorusPolynomial(N);
        MKtLweSymDecrypt(result_decrypt_m1, result_enc_m1, MKrlwekey, Msize);

        cout << "EXTERNAL PRODUCT v2m1" << endl;
        cout << "integer: " << sampleUE_clear << endl;
        cout << "result errors?" << endl;

        for (int i = 0; i < N; ++i)
        {
            if (modSwitchFromTorus32(result_decrypt_m1->coefsT[i], Msize) != modSwitchFromTorus32(result_clear->coefsT[i], Msize))
            {
                cout << modSwitchFromTorus32(result_decrypt_m1->coefsT[i], Msize) << " -- " << modSwitchFromTorus32(result_clear->coefsT[i], Msize) << endl;
            }
        }
        cout << endl; 







        MKTLweSample* resultFFT_enc_m1 = new_MKTLweSample(RLWEparams, MKparams);
        // External product FFT (inputs and outputs are of the same type but the function should be faster)
        clock_t begin_EP_m1_FFT = clock();
        MKtGswUEExternMulToMKtLwe_FFT_v2m1(resultFFT_enc_m1, sample_enc, sampleUE_FFT_enc, RLWEparams, MKparams, MKrlwekey);
        clock_t end_EP_m1_FFT = clock();
        cout << "Time per MK_extProd_FFT_v2_m1 (microsecs)... " << (end_EP_m1_FFT - begin_EP_m1_FFT) << endl;



        // Verify 
        TorusPolynomial *resultFFT_decrypt_m1 = new_TorusPolynomial(N);
        MKtLweSymDecrypt(resultFFT_decrypt_m1, resultFFT_enc_m1, MKrlwekey, Msize);

        cout << "EXTERNAL PRODUCT v2m1 FFT" << endl;
        cout << "integer: " << sampleUE_clear << endl;
        cout << "result errors?" << endl;

        for (int i = 0; i < N; ++i)
        {
            if (modSwitchFromTorus32(resultFFT_decrypt_m1->coefsT[i], Msize) != modSwitchFromTorus32(result_clear->coefsT[i], Msize))
            {
                cout << modSwitchFromTorus32(resultFFT_decrypt_m1->coefsT[i], Msize) << " -- " << modSwitchFromTorus32(result_clear->coefsT[i], Msize) << endl;
            }
        }
        cout << endl; 








        /* ********************************************************************************
        *********************** EXTERNAL PRODUCT method 2 *********************************
        ******************************************************************************** */


        MKTLweSample* result_enc_m2 = new_MKTLweSample(RLWEparams, MKparams);
        // External product 
        clock_t begin_EP_m2 = clock();
        MKtGswUEExternMulToMKtLwe_v2m2(result_enc_m2, sample_enc, sampleUE_enc, RLWEparams, MKparams, MKrlwekey);
        clock_t end_EP_m2 = clock();
        cout << "Time per MK_extProd_v2_m2 (microsecs)... " << (end_EP_m2 - begin_EP_m2) << endl;



        // Verify 
        TorusPolynomial *result_decrypt_m2 = new_TorusPolynomial(N);
        MKtLweSymDecrypt(result_decrypt_m2, result_enc_m2, MKrlwekey, Msize);

        cout << "EXTERNAL PRODUCT v2m2" << endl;
        cout << "integer: " << sampleUE_clear << endl;
        cout << "result errors?" << endl;

        for (int i = 0; i < N; ++i)
        {
            if (modSwitchFromTorus32(result_decrypt_m2->coefsT[i], Msize) != modSwitchFromTorus32(result_clear->coefsT[i], Msize))
            {
                cout << modSwitchFromTorus32(result_decrypt_m2->coefsT[i], Msize) << " -- " << modSwitchFromTorus32(result_clear->coefsT[i], Msize) << endl;
            }
        }
        cout << endl; 






        
        MKTLweSample* resultFFT_enc_m2 = new_MKTLweSample(RLWEparams, MKparams);
        // External product FFT (inputs and outputs are of the same type but the function should be faster)
        clock_t begin_EP_m2_FFT = clock();
        MKtGswUEExternMulToMKtLwe_FFT_v2m2(resultFFT_enc_m2, sample_enc, sampleUE_FFT_enc, RLWEparams, MKparams, MKrlwekey);
        clock_t end_EP_m2_FFT = clock();
        cout << "Time per MK_extProd_FFT_v2_m2 (microsecs)... " << (end_EP_m2_FFT - begin_EP_m2_FFT) << endl;



        // Verify 
        TorusPolynomial *resultFFT_decrypt_m2 = new_TorusPolynomial(N);
        MKtLweSymDecrypt(resultFFT_decrypt_m2, resultFFT_enc_m2, MKrlwekey, Msize);

        cout << "EXTERNAL PRODUCT v2m2 FFT" << endl;
        cout << "integer: " << sampleUE_clear << endl;
        cout << "result errors?" << endl;

        for (int i = 0; i < N; ++i)
        {
            if (modSwitchFromTorus32(resultFFT_decrypt_m2->coefsT[i], Msize) != modSwitchFromTorus32(result_clear->coefsT[i], Msize))
            {
                cout << modSwitchFromTorus32(resultFFT_decrypt_m2->coefsT[i], Msize) << " -- " << modSwitchFromTorus32(result_clear->coefsT[i], Msize) << endl;
            }
        }
        cout << endl; 
        















        // delete samples
        //delete_TorusPolynomial(resultFFT_decrypt_m2);
        //delete_MKTLweSample(resultFFT_enc_m2);
        delete_TorusPolynomial(result_decrypt_m2);
        delete_MKTLweSample(result_enc_m2);
        delete_TorusPolynomial(resultFFT_decrypt_m1);
        delete_TorusPolynomial(result_decrypt_m1);
        delete_TorusPolynomial(sample_decrypt);
        delete_MKTLweSample(resultFFT_enc_m1);
        delete_MKTLweSample(result_enc_m1);
        delete_TorusPolynomial(result_clear);
        delete_LagrangeHalfCPolynomial_array(3*dg, arrUE);
        delete_MKTGswUESampleFFT_v2(sampleUE_FFT_enc);
        delete_MKTGswUESample_v2(sampleUE_enc);
        delete_MKTLweSample(sample_enc);
        delete_TorusPolynomial(sample_clear);
    }

    
    cout << "Time per KEY GENERATION (seconds)... " << time_KG << endl;



    // delete keys
    delete_MKRLweKey(MKrlwekey);
    // delete params
    delete_MKTFHEParams(MKparams);
    delete_TLweParams(RLWEparams);
    delete_LweParams(LWEparams);
    delete_LweParams(extractedLWEparams);

    return 0;
}

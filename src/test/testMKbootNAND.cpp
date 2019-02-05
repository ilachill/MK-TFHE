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


//EXPORT void tLweExtractKey(LweKey* result, const TLweKey* key); //TODO: change the name and put in a .h
//EXPORT void tfhe_createLweBootstrappingKeyFFT(LweBootstrappingKeyFFT* bk, const LweKey* key_in, const TGswKey* rgsw_key);
//EXPORT void tfhe_bootstrapFFT(LweSample* result, const LweBootstrappingKeyFFT* bk, Torus32 mu1, Torus32 mu0, const LweSample* x);




int32_t main(int32_t argc, char **argv) {

    // Test trials
    const int32_t nb_trials = 10;

    // generate params 
    static const int32_t k = 1;
    //static const int32_t bk_l = 2;
    //static const int32_t bk_Bgbit = 10;
    //static const int32_t ks_basebit = 2;
    //static const int32_t ks_length = 8;
    static const double ks_stdev = 2.44e-5; //standard deviation
    static const double bk_stdev = 3.29e-10; //standard deviation
    static const double max_stdev = 0.012467; //max standard deviation for a 1/4 msg space
    // new params
    static const int32_t n = 500;            // LWE modulus
    static const int32_t n_extract = 1024;    // LWE extract modulus (used in bootstrapping)
    static const int32_t hLWE = 0;         // HW secret key LWE --> not used
    static const double stdevLWE = 0.012467;      // LWE ciphertexts standard deviation
    static const int32_t Bksbit = 2;       // Base bit key switching
    static const int32_t dks = 8;          // dimension key switching
    static const double stdevKS = 2.44e-5;       // KS key standard deviation
    static const int32_t N = 1024;            // RLWE,RGSW modulus
    static const int32_t hRLWE = 0;        // HW secret key RLWE,RGSW --> not used
    static const double stdevRLWEkey = 3.29e-10; // 0; // 0.012467;  // RLWE key standard deviation
    static const double stdevRLWE = 3.29e-10; // 0; // 0.012467;     // RLWE ciphertexts standard deviation
    static const double stdevRGSW = 3.29e-10;     // RGSW ciphertexts standard deviation 
    static const int32_t Bgbit = 7;        // Base bit gadget
    static const int32_t dg = 4;           // dimension gadget
    static const double stdevBK = 3.29e-10;       // BK standard deviation
    static const int32_t parties = 2;      // number of parties

    // 2 parties, B=2^7, d=4
    // 4 parties, B=2^6, d=5
    // 8 parties, B=2^4, d=8
    // noise 3.29e-10, security level 152 by following TFHE analysis


    // params
    LweParams *extractedLWEparams = new_LweParams(n_extract, ks_stdev, max_stdev);
    LweParams *LWEparams = new_LweParams(n, ks_stdev, max_stdev);
    TLweParams *RLWEparams = new_TLweParams(N, k, bk_stdev, max_stdev);
    MKTFHEParams *MKparams = new_MKTFHEParams(n, n_extract, hLWE, stdevLWE, Bksbit, dks, stdevKS, N, 
                            hRLWE, stdevRLWEkey, stdevRLWE, stdevRGSW, Bgbit, dg, stdevBK, parties);


    cout << "Params: DONE!" << endl;






    // Key generation 
    // LWE key
    MKLweKey* MKlwekey = new_MKLweKey(LWEparams, MKparams);
    MKLweKeyGen(MKlwekey);
    cout << "KeyGen MKlwekey: DONE!" << endl;
    // RLWE key 
    MKRLweKey* MKrlwekey = new_MKRLweKey(RLWEparams, MKparams);
    MKRLweKeyGen(MKrlwekey);
    cout << "KeyGen MKrlwekey: DONE!" << endl;
    // LWE key extracted 
    MKLweKey* MKextractedlwekey = new_MKLweKey(extractedLWEparams, MKparams);
    MKtLweExtractKey(MKextractedlwekey, MKrlwekey);
    cout << "KeyGen MKextractedlwekey: DONE!" << endl;
    // bootstrapping + key switching keys
    MKLweBootstrappingKey* MKlweBK = new_MKLweBootstrappingKey(LWEparams, RLWEparams, MKparams);
    MKlweCreateBootstrappingKey(MKlweBK, MKlwekey, MKrlwekey, MKextractedlwekey, 
                                extractedLWEparams, LWEparams, RLWEparams, MKparams);
    cout << "KeyGen MKlweBK: DONE!" << endl;
    cout << "KeyGen: DONE!" << endl;





    for (int trial = 0; trial < nb_trials; ++trial)
    {
        cout << "****************" << endl;
        cout << "Trial: " << trial << endl;
        cout << "****************" << endl;


        int32_t mess1 = rand() % 2;
        int32_t mess2 = rand() % 2;
        int32_t out = 1 - (mess1 * mess2);
        // generate 2 samples in input
        MKLweSample *test_in1 = new_MKLweSample(LWEparams, MKparams);
        MKLweSample *test_in2 = new_MKLweSample(LWEparams, MKparams);
        MKbootsSymEncrypt(test_in1, mess1, MKlwekey);
        MKbootsSymEncrypt(test_in2, mess2, MKlwekey);
        // generate output sample
        MKLweSample *test_out = new_MKLweSample(LWEparams, MKparams);

        cout << "Encryption: DONE!" << endl;




        // verify encrypt 
        cout << "Message 1: clear = " << mess1 << ", decrypted = " << MKbootsSymDecrypt(test_in1, MKlwekey) << endl;
        cout << "Message 2: clear = " << mess2 << ", decrypted = " << MKbootsSymDecrypt(test_in2, MKlwekey) << endl;






        // evaluate MK bootstrapped NAND 
        cout << "Starting MK bootstrapped NAND: trial " << trial << endl;
        clock_t begin = clock();
        MKbootsNAND(test_out, test_in1, test_in2, MKlweBK, LWEparams, extractedLWEparams, RLWEparams, MKparams);
        clock_t end = clock();
        cout << "Finished MK bootstrapped" << endl;
        cout << "Time per MKbootNAND gate (microsecs)... " << (end - begin) << endl;






        // verify NAND
        int32_t outNAND = MKbootsSymDecrypt(test_out, MKlwekey);
        cout << "NAND: clear = " << out << ", decrypted = " << outNAND << endl;
        if (outNAND != out) {
            cout << "ERROR!!! " << trial << "," << trial << " - ";
            cout << t32tod(MKlwePhase(test_in1, MKlwekey)) << " - ";
            cout << t32tod(MKlwePhase(test_in2, MKlwekey)) << " - ";
            cout << t32tod(MKlwePhase(test_out, MKlwekey)) << endl;
        }






        // delete samples
        delete_MKLweSample(test_out);
        delete_MKLweSample(test_in2);
        delete_MKLweSample(test_in1);
    }


   

    // delete keys
    delete_MKLweBootstrappingKey(MKlweBK);
    delete_MKLweKey(MKextractedlwekey);
    delete_MKRLweKey(MKrlwekey);
    delete_MKLweKey(MKlwekey);
    // delete params
    delete_MKTFHEParams(MKparams);
    delete_TLweParams(RLWEparams);
    delete_LweParams(LWEparams);
    delete_LweParams(extractedLWEparams);

    return 0;
}

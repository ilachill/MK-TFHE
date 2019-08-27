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


        



int32_t main(int32_t argc, char **argv) {

    // Test trials
    const int32_t nb_trials = 10;


    // generate params 
    static const int32_t k = 1;
    static const double ks_stdev = 3.05e-5;// 2.44e-5; //standard deviation
    static const double bk_stdev = 3.72e-9; // 3.29e-10; //standard deviation
    static const double max_stdev = 0.012467; //max standard deviation for a 1/4 msg space
    static const int32_t n = 560; //500;            // LWE modulus
    static const int32_t n_extract = 1024;    // LWE extract modulus (used in bootstrapping)
    static const int32_t hLWE = 0;         // HW secret key LWE --> not used
    static const double stdevLWE = 0.012467;      // LWE ciphertexts standard deviation
    static const int32_t Bksbit = 2;       // Base bit key switching
    static const int32_t dks = 8;          // dimension key switching
    static const double stdevKS = ks_stdev; // 2.44e-5;       // KS key standard deviation
    static const int32_t N = 1024;            // RLWE,RGSW modulus
    static const int32_t hRLWE = 0;        // HW secret key RLWE,RGSW --> not used
    static const double stdevRLWEkey = bk_stdev; // 3.29e-10; // 0; // 0.012467;  // RLWE key standard deviation
    static const double stdevRLWE = bk_stdev; // 3.29e-10; // 0; // 0.012467;     // RLWE ciphertexts standard deviation
    static const double stdevRGSW = bk_stdev; // 3.29e-10;     // RGSW ciphertexts standard deviation 
    static const int32_t Bgbit = 9;        // Base bit gadget
    static const int32_t dg = 3;           // dimension gadget
    static const double stdevBK = bk_stdev; // 3.29e-10;       // BK standard deviation
    static const int32_t parties = 2;      // number of parties

    // new parameters 
    // 2 parties, B=2^9, d=3 -> works
    // 4 parties, B=2^8, d=4 -> works
    // 8 parties, B=2^6, d=5 -> works 
    

    // params
    LweParams *extractedLWEparams = new_LweParams(n_extract, ks_stdev, max_stdev);
    LweParams *LWEparams = new_LweParams(n, ks_stdev, max_stdev);
    TLweParams *RLWEparams = new_TLweParams(N, k, bk_stdev, max_stdev);
    MKTFHEParams *MKparams = new_MKTFHEParams(n, n_extract, hLWE, stdevLWE, Bksbit, dks, stdevKS, N, 
                            hRLWE, stdevRLWEkey, stdevRLWE, stdevRGSW, Bgbit, dg, stdevBK, parties);


    cout << "Params: DONE!" << endl;






   
    // Key generation 
    cout << "Starting KEY GENERATION" << endl;
    clock_t begin_KG = clock();

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
    MKLweBootstrappingKey_v2* MKlweBK = new_MKLweBootstrappingKey_v2(LWEparams, RLWEparams, MKparams);
    MKlweCreateBootstrappingKey_v2(MKlweBK, MKlwekey, MKrlwekey, MKextractedlwekey, 
                                extractedLWEparams, LWEparams, RLWEparams, MKparams);
    cout << "KeyGen MKlweBK: DONE!" << endl;

    // bootstrapping FFT + key switching keys
    MKLweBootstrappingKeyFFT_v2* MKlweBK_FFT = new_MKLweBootstrappingKeyFFT_v2(MKlweBK, LWEparams, RLWEparams, MKparams);
    cout << "KeyGen MKlweBK_FFT: DONE!" << endl;   

    clock_t end_KG = clock();
    double time_KG = ((double) end_KG - begin_KG)/CLOCKS_PER_SEC;
    cout << "Finished KEY GENERATION" << endl;





    



    int32_t error_count_EncDec = 0;
    
    int32_t error_count_v2m2 = 0;
    double argv_time_NAND_v2m2 = 0.0;




    for (int trial = 0; trial < nb_trials; ++trial)
    {
        cout << "****************" << endl;
        cout << "Trial: " << trial << endl;
        cout << "****************" << endl;

        // use current time as seed for the random generator
        srand(time(0));

        int32_t mess1 = rand() % 2;
        int32_t mess2 = rand() % 2;
        int32_t out = 1 - (mess1 * mess2);
        // generate 2 samples in input
        MKLweSample *test_in1 = new_MKLweSample(LWEparams, MKparams);
        MKLweSample *test_in2 = new_MKLweSample(LWEparams, MKparams);
        MKbootsSymEncrypt(test_in1, mess1, MKlwekey);
        MKbootsSymEncrypt(test_in2, mess2, MKlwekey);
        // generate output sample
        MKLweSample *test_out_v2m2 = new_MKLweSample(LWEparams, MKparams);

        cout << "Encryption: DONE!" << endl;




        // verify encrypt 
        int32_t mess1_dec = MKbootsSymDecrypt(test_in1, MKlwekey);
        int32_t mess2_dec = MKbootsSymDecrypt(test_in2, MKlwekey);
        cout << "Message 1: clear = " << mess1 << ", decrypted = " << mess1_dec << endl;
        cout << "Message 2: clear = " << mess2 << ", decrypted = " << mess2_dec << endl;

        // count encrypt/decrypt errors
        if (mess1 != mess1_dec)
        {
            error_count_EncDec += 1;
        }
        if (mess2 != mess2_dec)
        {
            error_count_EncDec += 1;
        }








        // evaluate MK bootstrapped NAND 
        cout << "Starting MK bootstrapped NAND FFT version 2 method 2: trial " << trial << endl;
        clock_t begin_NAND_v2m2 = clock();
        MKbootsNAND_FFT_v2m2(test_out_v2m2, test_in1, test_in2, MKlweBK_FFT, LWEparams, extractedLWEparams, RLWEparams, MKparams, MKrlwekey);
        clock_t end_NAND_v2m2 = clock();
        double time_NAND_v2m2 = ((double) end_NAND_v2m2 - begin_NAND_v2m2)/CLOCKS_PER_SEC;
        cout << "Finished MK bootstrapped NAND FFT v2m2" << endl;
        cout << "Time per MKbootNAND_FFT gate v2m2 (seconds)... " << time_NAND_v2m2 << endl;

        argv_time_NAND_v2m2 += time_NAND_v2m2;

        // verify NAND
        int32_t outNAND_v2m2 = MKbootsSymDecrypt(test_out_v2m2, MKlwekey);
        cout << "NAND: clear = " << out << ", decrypted = " << outNAND_v2m2 << endl;
        if (outNAND_v2m2 != out) {
            error_count_v2m2 +=1;
            cout << "ERROR!!! " << trial << "," << trial << " - ";
            cout << t32tod(MKlwePhase(test_in1, MKlwekey)) << " - ";
            cout << t32tod(MKlwePhase(test_in2, MKlwekey)) << " - ";
            cout << t32tod(MKlwePhase(test_out_v2m2, MKlwekey)) << endl;
        }








        // delete samples
        delete_MKLweSample(test_out_v2m2);
        delete_MKLweSample(test_in2);
        delete_MKLweSample(test_in1);
    }

    cout << endl;
    cout << "Time per KEY GENERATION (seconds)... " << time_KG << endl;
    
    cout << "ERRORS v2m2: " << error_count_v2m2 << " over " << nb_trials << " tests!" << endl;
    cout << "Average time per bootNAND_FFT_v2m2: " << argv_time_NAND_v2m2/nb_trials << " seconds" << endl;

    cout << endl << "ERRORS Encrypt/Decrypt: " << error_count_EncDec << " over " << nb_trials << " tests!" << endl;
    

   

    // delete keys
    delete_MKLweBootstrappingKeyFFT_v2(MKlweBK_FFT);
    delete_MKLweBootstrappingKey_v2(MKlweBK);
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

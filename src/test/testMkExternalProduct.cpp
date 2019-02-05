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
    const int32_t nb_trials = 10;

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
        MKTGswUESample* sampleUE_enc = new_MKTGswUESample(RLWEparams, MKparams);
        MKTGswUniEncryptI(sampleUE_enc, sampleUE_clear, 0, stdevRGSW, MKrlwekey); // party = 0

        // sample TRGSW_Exp
        MKTGswExpSample* sampleExp_enc = new_MKTGswExpSample(RLWEparams, MKparams);
        MKTGswExpand(sampleExp_enc, sampleUE_enc, MKrlwekey, MKparams);
        /*
        MKTGswExpSampleFFT* sampleExpFFT_enc = new_MKTGswExpSampleFFT(RLWEparams, MKparams, 0.0);
        MKTGswExpandFFT(sampleExpFFT_enc, sampleUE_enc, MKrlwekey, RLWEparams, MKparams);
        */

        // result RLWE
        TorusPolynomial *result_clear = new_TorusPolynomial(N);
        MKTLweSample* result_enc = new_MKTLweSample(RLWEparams, MKparams);

        




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
        TorusPolynomial *sampleUE_decrypt = new_TorusPolynomial_array(3*dg, N);
        MKtGswSymDecrypt(sampleUE_decrypt, sampleUE_enc, MKrlwekey);
        cout << "sampleUE decrypt: clear = " << sampleUE_clear << endl;
        cout << "sampleUE decrypt: errors? " << endl;  
        vector<int32_t> plain(N); 
        plain[0] = sampleUE_clear;
        for (int j = 1; j < N; j++){
            plain[j] = 0;
        }   
        vector<int32_t> plain_times_secret(N); 
        for (int j = 0; j < N; j++){
            plain_times_secret[j] = sampleUE_clear * MKrlwekey->key[sampleUE_enc->party].key->coefs[j]; 
        }   

        for (int i = 0; i < dg; ++i)
        {
            int32_t MsizeBg = 1 << (Bgbit*(i+1)); // Bg^{i+1}


            // c part
            cout << "c part" << endl;
            print_error(sampleUE_decrypt[i], plain, MsizeBg); 
            // first coeff = message 
            // if (modSwitchFromTorus32(sampleUE_decrypt[i].coefsT[0], MsizeBg) != sampleUE_clear)
            // {
            //     cout << "Error: i=" << i << ", j=" << 0 << " - " << modSwitchFromTorus32(sampleUE_decrypt[i].coefsT[0], MsizeBg) << ", ";
            // }
            // for (int j = 1; j < N; ++j)
            // {
            //     // everithing equal to 0
            //     if (modSwitchFromTorus32(sampleUE_decrypt[i].coefsT[j], MsizeBg) != 0)
            //     {
            //         cout << "Error: i=" << i << ", j=" << j << " - " << modSwitchFromTorus32(sampleUE_decrypt[i].coefsT[j], MsizeBg) << ", ";
            //     }                    
            // }
            // cout << endl;

            // d part
            cout << "d part" << endl;
            // -d0 + d1 * s[party] ~= mu*g * s[party]
            
            print_error(sampleUE_decrypt[dg + i], plain_times_secret, MsizeBg); 

            // for (int j = 0; j < N; ++j)
            // {   
            //     // equal to the secret key 

            //     if (modSwitchFromTorus32(sampleUE_decrypt[dg + i].coefsT[j], MsizeBg) != sampleUE_clear * MKrlwekey->key[sampleUE_enc->party].key->coefs[j])
            //     {
            //         cout << "Error: i=" << i << ", j=" << j << " - " << modSwitchFromTorus32(sampleUE_decrypt[dg + i].coefsT[j], MsizeBg) << " - " << MKrlwekey->key[sampleUE_enc->party].key->coefs[j] << ", ";
            //     }                   
            // }
            // cout << endl;

            /*
            // f part : just print it and compare it with the printed r before 
            cout << "f part (r)" << endl;
            for (int j = 0; j < N; ++j)
            {
                cout << modSwitchFromTorus32(sampleUE_decrypt[2*dg + i].coefsT[j], 2) << ", ";
            }
            cout << endl;
            */
        }             







        // sample RGSW Expand
        TorusPolynomial *sampleExp_decrypt = new_TorusPolynomial_array((parties+1)*dg, N);
        MKtGswEXPSymDecrypt(sampleExp_decrypt, sampleExp_enc, MKrlwekey);
        cout << "sampleExp decrypt: clear = " << sampleUE_clear << endl;
        cout << "sampleExp decrypt: errors? " << endl;     
        // verify that: 
        // for i= 0,...,parties-1, phi_i ~ mu*s_i
        for (int i = 0; i < parties; ++i)
        {
            for (int j = 0; j < dg; ++j)
            {
                int32_t MsizeBg = 1 << (Bgbit*(j+1)); // Bg^{j+1}

                double maxerr = 0; 
                for (int h = 0; h < N; ++h)
                {
                    //if (modSwitchFromTorus32(sampleExp_decrypt[i*dg+j].coefsT[h], MsizeBg) != sampleUE_clear * MKrlwekey->key[i].key->coefs[h])
                    //{
                        //cout << "Error " << i << "," << j << "," << h << " - ";
                        // cout << modSwitchFromTorus32(sampleExp_decrypt[i*dg+j].coefsT[h], MsizeBg) << ";"; 
                        double a = (double) sampleExp_decrypt[i*dg+j].coefsT[h] / pow(2,32);
                        double b = (double) sampleUE_clear * MKrlwekey->key[i].key->coefs[h] / (double) MsizeBg;
                        //cout << abs(a-b) << endl;
                        maxerr = max(maxerr, abs(a-b)); 
                    //}
                }
                cout << "Error " << i << "," << j << " : " << maxerr << endl;

            }
        }
        cout << endl;
        // verify that: 
        // phi_parties = phi_b ~ mu
        for (int j = 0; j < dg; ++j)
        {
            int32_t MsizeBg = 1 << (Bgbit*(j+1)); // Bg^{j+1}


            double maxerr = 0; 
            //if (modSwitchFromTorus32(sampleExp_decrypt[parties*dg+j].coefsT[0], MsizeBg) != sampleUE_clear)
            //{
            double a = (double) sampleExp_decrypt[parties*dg+j].coefsT[0] / pow(2,32);
            double b = (double) sampleUE_clear / (double) MsizeBg;
            maxerr = max(maxerr, abs(a-b)); 
            for (int h = 1; h < N; ++h)
            {
                double a = (double) sampleExp_decrypt[parties*dg+j].coefsT[h] / pow(2,32);
                double b = 0; // msg = 0 for higher coefficients
                maxerr = max(maxerr, abs(a-b)); 
                // everithing equal to 0
                //if (modSwitchFromTorus32(sampleExp_decrypt[parties*dg+j].coefsT[h], MsizeBg) != 0)
                //{
                //    cout << "Error " << parties << "," << j << "," << h << " - ";
                //}                    
            }
            cout << "Error " << parties << "," << j << ":" << maxerr << endl;

        }
        cout << endl;

        





























        // External FFT product 
        MKtGswExpExternMulToMKTLwe(result_enc, sample_enc, sampleExp_enc, MKparams);
        /*
        MKtGswExpExternMulToMKTLweFFT(result_enc, sample_enc, sampleExpFFT_enc, RLWEparams, MKparams);
        */

        // Multiplication in clear
        for (int i = 0; i < N; ++i)
        {
            result_clear->coefsT[i] = sampleUE_clear*sample_clear->coefsT[i];
        }



        // Verify 
        TorusPolynomial *result_decrypt = new_TorusPolynomial(N);
        MKtLweSymDecrypt(result_decrypt, result_enc, MKrlwekey, Msize);

        cout << "integer: " << sampleUE_clear << endl;
        cout << "result errors?" << endl;
            
        for (int i = 0; i < N; ++i)
        {
            int32_t res_dec = modSwitchFromTorus32(result_decrypt->coefsT[i], Msize);
            int32_t res_cle = modSwitchFromTorus32(result_clear->coefsT[i], Msize);
            if (res_dec != res_cle)
            {
                cout << "error " << i << ", ";
                //cout << res_dec << " -- " << res_cle << ", ";
            }
        }
        cout << endl; 
    











        
        // test decomp
        cout << endl << "TEST DECOMP" << endl;

        IntPolynomial* u = new_IntPolynomial_array(dg, N);
        TorusPolynomial *Tpoly = new_TorusPolynomial(N);
        torusPolynomialUniform(Tpoly);

        // decomp_g
        //MKtGswTorus32PolynomialDecompG(u, Tpoly, MKparams);
        MKtGswTorus32PolynomialDecompGassembly(u, Tpoly, MKparams); 
        //tGswTorus32PolynomialDecompH(u, Tpoly, TGSWparams);

        
        // g*decomp_g(Tpoly) ?= Tpoly
        int32_t eps = 1 << (32 - Bgbit*dg); // 1/(Bg^dg)
        int32_t beta = 1 << (Bgbit - 1); // Bg/2

        cout << "Checking MK params ..." << endl;
        for (int i = 0; i < dg; i++){
            cout << log2(MKparams->g[i]) << ",";
        }
        cout << endl;

        for (int j = 0; j < N; ++j)
        {
            Torus32 Tpoly_out = 0;

            for (int l = 0; l < dg; ++l)
            {
                if (abs(u[l].coefs[j]) > beta)
                {
                    cout << "Error decomp: " << u[l].coefs[j] << ", ";
                }

                Tpoly_out += u[l].coefs[j] * MKparams->g[l];
            }   

            if (abs(Tpoly->coefsT[j] - Tpoly_out) > eps)
            {
                cout << "j=" << j << ": " << abs(Tpoly_out - Tpoly->coefsT[j]) << endl;
            }      
        }
        cout << endl;
        
        delete_TorusPolynomial(Tpoly);
        delete_IntPolynomial_array(dg, u);

        








        // delete samples
        delete_TorusPolynomial(result_decrypt);
        delete_TorusPolynomial_array((parties+1)*dg, sampleExp_decrypt);
        delete_TorusPolynomial_array(3*dg, sampleUE_decrypt);
        delete_TorusPolynomial(sample_decrypt);
        delete_MKTLweSample(result_enc);
        delete_TorusPolynomial(result_clear);
        
        delete_MKTGswExpSample(sampleExp_enc);
        /*
        delete_MKTGswExpSampleFFT(sampleExpFFT_enc);
        */
        delete_MKTGswUESample(sampleUE_enc);
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

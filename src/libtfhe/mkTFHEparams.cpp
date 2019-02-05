#include <cmath>
#include "mkTFHEparams.h"
#include <new>
#include <stdlib.h>


using namespace std;

MKTFHEParams::MKTFHEParams(int32_t n, int32_t n_extract, int32_t hLWE, double stdevLWE, 
							int32_t Bksbit, int32_t dks, double stdevKS,
							int32_t N, int32_t hRLWE, double stdevRLWEkey, double stdevRLWE, double stdevRGSW, 
							int32_t Bgbit, int32_t dg, double stdevBK,
							int32_t parties):

		n(n),						// LWE modulus
		n_extract(n_extract),		// LWE extract modulus (used in bootstrapping)
		hLWE(hLWE),					// HW secret key LWE
		stdevLWE(stdevLWE),			// LWE ciphertexts standard deviation
		Bksbit(Bksbit),				// Base bit key switching
		dks(dks),					// dimension key switching
		stdevKS(stdevKS),			// KS key standard deviation
		N(N), 						// RLWE,RGSW modulus
		hRLWE(hRLWE), 				// HW secret key RLWE,RGSW
		stdevRLWEkey(stdevRLWEkey),	// RLWE key standard deviation
		stdevRLWE(stdevRLWE), 		// RLWE ciphertexts standard deviation
		stdevRGSW(stdevRGSW),		// RGSW ciphertexts standard deviation 
		Bgbit(Bgbit),				// Base bit gadget
		dg(dg),						// dimension gadget
		stdevBK(stdevBK),			// BK standard deviation
		parties(parties),			// number of parties

		maskMod((1 << Bgbit) - 1),	// Bg - 1
    	halfBg((1 << Bgbit) / 2)	// Bg / 2
{
	g = new Torus32[dg]; 			// gadget
    for (int i = 0; i < dg; ++i) {
        int32_t kk = (32 - (i + 1) * Bgbit);
        g[i] = 1 << kk; // 1/(Bg^(i+1)) as a Torus32
    }

    // offset = Bg/2 * (2^(32-Bgbit) + 2^(32-2*Bgbit) + ... + 2^(32-l*Bgbit))
    uint32_t temp1 = 0;
    for (int i = 0; i < dg; ++i) {
        uint32_t temp0 = 1 << (32 - (i + 1) * Bgbit);
        temp1 += temp0;
    }
    offset = temp1 * halfBg;
}

MKTFHEParams::~MKTFHEParams() {
	delete[] g;
}




// alloc
EXPORT MKTFHEParams* alloc_MKTFHEParams() {
    return (MKTFHEParams*) malloc(sizeof(MKTFHEParams));
}
EXPORT MKTFHEParams* alloc_MKTFHEParams_array(int32_t nbelts) {
    return (MKTFHEParams*) malloc(nbelts*sizeof(MKTFHEParams));
}

// free memory space 
EXPORT void free_MKTFHEParams(MKTFHEParams* ptr) {
    free(ptr);
}
EXPORT void free_MKTFHEParams_array(int32_t nbelts, MKTFHEParams* ptr) {
    free(ptr);
}

// initialize the key structure
EXPORT void init_MKTFHEParams(MKTFHEParams* obj, int32_t n, int32_t n_extract, int32_t hLWE, double stdevLWE, 
		int32_t Bksbit, int32_t dks, double stdevKS, int32_t N, int32_t hRLWE, double stdevRLWEkey, 
		double stdevRLWE, double stdevRGSW, int32_t Bgbit, int32_t dg, double stdevBK, int32_t parties) 
{
    new(obj) MKTFHEParams(n, n_extract, hLWE, stdevLWE, Bksbit, dks, stdevKS, N, hRLWE, stdevRLWEkey, 
		stdevRLWE, stdevRGSW, Bgbit, dg, stdevBK, parties);
}
EXPORT void init_MKTFHEParams_array(int32_t nbelts, MKTFHEParams* obj, int32_t n, int32_t n_extract, int32_t hLWE, 
		double stdevLWE, int32_t Bksbit, int32_t dks, double stdevKS, int32_t N, int32_t hRLWE, double stdevRLWEkey, 
		double stdevRLWE, double stdevRGSW, int32_t Bgbit, int32_t dg, double stdevBK, int32_t parties) 
{
    for (int i = 0; i < nbelts; i++) {
    	new(obj+i) MKTFHEParams(n, n_extract, hLWE, stdevLWE, Bksbit, dks, stdevKS, N, hRLWE, stdevRLWEkey, 
		stdevRLWE, stdevRGSW, Bgbit, dg, stdevBK, parties);
    }
}

// destroys the structure
EXPORT void destroy_MKTFHEParams(MKTFHEParams* obj) {
    obj->~MKTFHEParams();
}
EXPORT void destroy_MKTFHEParams_array(int32_t nbelts, MKTFHEParams* obj) {
    for (int i = 0; i < nbelts; i++) 
    {
		(obj+i)->~MKTFHEParams();
    }
}
 
// new = alloc + init
EXPORT MKTFHEParams* new_MKTFHEParams(int32_t n, int32_t n_extract, int32_t hLWE, double stdevLWE, 
		int32_t Bksbit, int32_t dks, double stdevKS, int32_t N, int32_t hRLWE, double stdevRLWEkey, 
		double stdevRLWE, double stdevRGSW, int32_t Bgbit, int32_t dg, double stdevBK, int32_t parties) 
{
    return new MKTFHEParams(n, n_extract, hLWE, stdevLWE, Bksbit, dks, stdevKS, N, hRLWE, stdevRLWEkey, 
		stdevRLWE, stdevRGSW, Bgbit, dg, stdevBK, parties);
}
EXPORT MKTFHEParams* new_MKTFHEParams_array(int32_t nbelts, int32_t n, int32_t n_extract, int32_t hLWE, double stdevLWE, 
		int32_t Bksbit, int32_t dks, double stdevKS, int32_t N, int32_t hRLWE, double stdevRLWEkey, 
		double stdevRLWE, double stdevRGSW, int32_t Bgbit, int32_t dg, double stdevBK, int32_t parties) 
{
    MKTFHEParams* obj = alloc_MKTFHEParams_array(nbelts);
    init_MKTFHEParams_array(nbelts, obj, n, n_extract, hLWE, stdevLWE, Bksbit, dks, stdevKS, N, hRLWE, stdevRLWEkey, 
		stdevRLWE, stdevRGSW, Bgbit, dg, stdevBK, parties);
    return obj;
}

// delete = destroy + free
EXPORT void delete_MKTFHEParams(MKTFHEParams* obj) {
    delete obj;
}
EXPORT void delete_MKTFHEParams_array(int32_t nbelts, MKTFHEParams* obj) {
    destroy_MKTFHEParams_array(nbelts,obj);
    free_MKTFHEParams_array(nbelts,obj);
}





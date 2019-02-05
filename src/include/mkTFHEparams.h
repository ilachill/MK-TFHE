#ifndef MKTFHEPARAMS_H
#define MKTFHEPARAMS_H


#include "tfhe_core.h"

//this structure contains MK parameters
//this structure is constant (cannot be modified once initialized): 
//the pointer to the param can be passed directly
//to all the keys that use these params.
struct MKTFHEParams{
	const int32_t n;			// LWE modulus
	const int32_t n_extract;	// LWE extract modulus (used in bootstrapping)
	const int32_t hLWE;			// HW secret key LWE
	const double stdevLWE;		// LWE ciphertexts standard deviation
	const int32_t Bksbit;		// Base bit key switching
	const int32_t dks;			// dimension key switching
	const double stdevKS;		// KS key standard deviation
	const int32_t N;			// RLWE,RGSW modulus
	const int32_t hRLWE;		// HW secret key RLWE,RGSW
	const double stdevRLWEkey;	// RLWE key standard deviation
	const double stdevRLWE;		// RLWE ciphertexts standard deviation
	const double stdevRGSW;		// RGSW ciphertexts standard deviation 
	const int32_t Bgbit;		// Base bit gadget
	const int32_t dg;			// dimension gadget
	const double stdevBK;		// BK standard deviation
	const int32_t parties;		// number of parties
	const uint32_t maskMod;		// Bg - 1
    const int32_t halfBg; 		// Bg/2
	Torus32 *g; 				// gadget: powers of Bgbit
    uint32_t offset; 			// offset = Bg/2 * (2^(32-Bgbit) + 2^(32-2*Bgbit) + ... + 2^(32-l*Bgbit))

//since all members are declared constant, a constructor is 
//required in the structure.
#ifdef __cplusplus
	MKTFHEParams(int32_t n, int32_t n_extract, int32_t hLWE, double stdevLWE, 
				int32_t Bksbit, int32_t dks, double stdevKS,
				int32_t N, int32_t hRLWE, double stdevRLWEkey, double stdevRLWE, double stdevRGSW, 
				int32_t Bgbit, int32_t dg, double stdevBK,
				int32_t parties);
	~MKTFHEParams();
	MKTFHEParams(const MKTFHEParams&) = delete; //forbidden
	MKTFHEParams& operator=(const MKTFHEParams& ) = delete; //forbidden
#endif
};




// alloc
EXPORT MKTFHEParams* alloc_MKTFHEParams();
EXPORT MKTFHEParams* alloc_MKTFHEParams_array(int32_t nbelts);
// free memory space 
EXPORT void free_MKTFHEParams(MKTFHEParams* ptr);
EXPORT void free_MKTFHEParams_array(int32_t nbelts, MKTFHEParams* ptr);
// initialize the key structure
EXPORT void init_MKTFHEParams(MKTFHEParams* obj, int32_t n, int32_t n_extract, int32_t hLWE, double stdevLWE, 
		int32_t Bksbit, int32_t dks, double stdevKS, int32_t N, int32_t hRLWE, double stdevRLWEkey, 
		double stdevRLWE, double stdevRGSW, int32_t Bgbit, int32_t dg, double stdevBK, int32_t parties) ;
EXPORT void init_MKTFHEParams_array(int32_t nbelts, MKTFHEParams* obj, int32_t n, int32_t n_extract, int32_t hLWE, 
		double stdevLWE, int32_t Bksbit, int32_t dks, double stdevKS, int32_t N, int32_t hRLWE, double stdevRLWEkey, 
		double stdevRLWE, double stdevRGSW, int32_t Bgbit, int32_t dg, double stdevBK, int32_t parties);
// destroys the structure
EXPORT void destroy_MKTFHEParams(MKTFHEParams* obj);
EXPORT void destroy_MKTFHEParams_array(int32_t nbelts, MKTFHEParams* obj);
// new = alloc + init
EXPORT MKTFHEParams* new_MKTFHEParams(int32_t n, int32_t n_extract, int32_t hLWE, double stdevLWE, 
		int32_t Bksbit, int32_t dks, double stdevKS, int32_t N, int32_t hRLWE, double stdevRLWEkey, 
		double stdevRLWE, double stdevRGSW, int32_t Bgbit, int32_t dg, double stdevBK, int32_t parties);
EXPORT MKTFHEParams* new_MKTFHEParams_array(int32_t nbelts, int32_t n, int32_t n_extract, int32_t hLWE, double stdevLWE, 
		int32_t Bksbit, int32_t dks, double stdevKS, int32_t N, int32_t hRLWE, double stdevRLWEkey, 
		double stdevRLWE, double stdevRGSW, int32_t Bgbit, int32_t dg, double stdevBK, int32_t parties);
// delete = destroy + free
EXPORT void delete_MKTFHEParams(MKTFHEParams* obj);
EXPORT void delete_MKTFHEParams_array(int32_t nbelts, MKTFHEParams* obj);

#endif //MKTFHEPARAMS_H

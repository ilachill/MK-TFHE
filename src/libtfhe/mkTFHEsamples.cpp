#include "tfhe_core.h"
#include "lwesamples.h"
#include "lweparams.h"
#include "tlwe.h"
#include "polynomials.h"
#include "tgsw.h"
#include "polynomials_arithmetic.h"
#include "lagrangehalfc_arithmetic.h"
#include <new>
#include <stdlib.h>


#include "mkTFHEparams.h"
#include "mkTFHEkeys.h"
#include "mkTFHEsamples.h"


// MK LWE sample (a_1, ..., a_k, b)
MKLweSample::MKLweSample(const LweParams* LWEparams, const MKTFHEParams* MKparams) : 
		parties(MKparams->parties), n(LWEparams->n)
{
	this->a = new Torus32[parties*n];
    this->b = 0;
    this->current_variance = 0.0;
}

MKLweSample::~MKLweSample() {
    delete[] a;
}



// alloc 
EXPORT MKLweSample* alloc_MKLweSample(){
    return (MKLweSample*) malloc(sizeof(MKLweSample));
}
EXPORT MKLweSample* alloc_MKLweSample_array(int32_t nbelts) {
    return (MKLweSample*) malloc(nbelts*sizeof(MKLweSample));
}

//free memory space 
EXPORT void free_MKLweSample(MKLweSample* ptr) {
    free(ptr);
}
EXPORT void free_MKLweSample_array(int32_t nbelts, MKLweSample* ptr){
    free(ptr);
}

// init
EXPORT void init_MKLweSample(MKLweSample* obj, const LweParams* LWEparams, const MKTFHEParams* MKparams) {
    new(obj) MKLweSample(LWEparams, MKparams);
}
EXPORT void init_MKLweSample_array(int32_t nbelts, MKLweSample* obj, const LweParams* LWEparams, 
        const MKTFHEParams* MKparams)
{
    for (int i = 0; i < nbelts; i++) {
        new(obj+i) MKLweSample(LWEparams, MKparams);
    }
}

// destroys the structure
EXPORT void destroy_MKLweSample(MKLweSample* obj) {
    obj->~MKLweSample();
}
EXPORT void destroy_MKLweSample_array(int32_t nbelts, MKLweSample* obj) {
    for (int i = 0; i < nbelts; i++) {
        (obj+i)->~MKLweSample();
    }
}

// new = alloc + init
EXPORT MKLweSample* new_MKLweSample(const LweParams* LWEparams, const MKTFHEParams* MKparams) {
    return new MKLweSample(LWEparams, MKparams);
}
EXPORT MKLweSample* new_MKLweSample_array(int32_t nbelts, const LweParams* LWEparams, const MKTFHEParams* MKparams){
    MKLweSample* obj = alloc_MKLweSample_array(nbelts);
    init_MKLweSample_array(nbelts,obj,LWEparams,MKparams);
    return obj;
}

// delete = destroy + free
EXPORT void delete_MKLweSample(MKLweSample* obj) {
    delete obj;
}
EXPORT void delete_MKLweSample_array(int32_t nbelts, MKLweSample* obj) {
    destroy_MKLweSample_array(nbelts,obj);
    free_MKLweSample_array(nbelts,obj);
}

































// MK RLWE sample (a_1, ..., a_k, b)
MKTLweSample::MKTLweSample(const TLweParams* RLWEparams, const MKTFHEParams* MKparams) :
		parties(MKparams->parties), N(RLWEparams->N)
{
	//a is a table of parties+1 polynomials, b is an alias for &a[parties]
    a = new_TorusPolynomial_array(parties+1, N);
    b = a + parties;
    current_variance = 0.0;
}

MKTLweSample::~MKTLweSample() {
    delete_TorusPolynomial_array(parties+1, a);
}


// alloc 
EXPORT MKTLweSample* alloc_MKTLweSample(){
    return (MKTLweSample*) malloc(sizeof(MKTLweSample));
}
EXPORT MKTLweSample* alloc_MKTLweSample_array(int32_t nbelts) {
    return (MKTLweSample*) malloc(nbelts*sizeof(MKTLweSample));
}

//free memory space 
EXPORT void free_MKTLweSample(MKTLweSample* ptr) {
    free(ptr);
}
EXPORT void free_MKTLweSample_array(int32_t nbelts, MKTLweSample* ptr){
    free(ptr);
}

// init
EXPORT void init_MKTLweSample(MKTLweSample* obj, const TLweParams* RLWEparams, const MKTFHEParams* MKparams) {
    new(obj) MKTLweSample(RLWEparams, MKparams);
}
EXPORT void init_MKTLweSample_array(int32_t nbelts, MKTLweSample* obj, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams)
{
    for (int i = 0; i < nbelts; i++) {
        new(obj+i) MKTLweSample(RLWEparams, MKparams);
    }
}

// destroys the structure
EXPORT void destroy_MKTLweSample(MKTLweSample* obj) {
    obj->~MKTLweSample();
}
EXPORT void destroy_MKTLweSample_array(int32_t nbelts, MKTLweSample* obj) {
    for (int i = 0; i < nbelts; i++) {
        (obj+i)->~MKTLweSample();
    }
}

// new = alloc + init
EXPORT MKTLweSample* new_MKTLweSample(const TLweParams* RLWEparams, const MKTFHEParams* MKparams) {
    MKTLweSample* obj = alloc_MKTLweSample();
    init_MKTLweSample(obj,RLWEparams,MKparams);
    return obj;
    //return new MKTLweSample(RLWEparams, MKparams);
}
EXPORT MKTLweSample* new_MKTLweSample_array(int32_t nbelts, const TLweParams* RLWEparams, const MKTFHEParams* MKparams){
    MKTLweSample* obj = alloc_MKTLweSample_array(nbelts);
    init_MKTLweSample_array(nbelts,obj,RLWEparams,MKparams);
    return obj;
}

// delete = destroy + free
EXPORT void delete_MKTLweSample(MKTLweSample* obj) {
    destroy_MKTLweSample(obj);
    free_MKTLweSample(obj);
    //delete obj;
}
EXPORT void delete_MKTLweSample_array(int32_t nbelts, MKTLweSample* obj) {
    destroy_MKTLweSample_array(nbelts,obj);
    free_MKTLweSample_array(nbelts,obj);
}































// MK RLWE sample FFT (a_1, ..., a_k, b)
MKTLweSampleFFT::MKTLweSampleFFT(const TLweParams *params, const MKTFHEParams* MKparams, LagrangeHalfCPolynomial *arr, 
		double current_variance) : parties(MKparams->parties)
{
    //a is a table of parties+1 polynomials, b is an alias for &a[parties]
    a = arr;
    b = a + parties;
    current_variance = 0.0;
}

MKTLweSampleFFT::~MKTLweSampleFFT() {
}



// alloc 
EXPORT MKTLweSampleFFT* alloc_MKTLweSampleFFT(){
    return (MKTLweSampleFFT*) malloc(sizeof(MKTLweSampleFFT));
}
EXPORT MKTLweSampleFFT* alloc_MKTLweSampleFFT_array(int32_t nbelts) {
    return (MKTLweSampleFFT*) malloc(nbelts*sizeof(MKTLweSampleFFT));
}

//free memory space 
EXPORT void free_MKTLweSampleFFT(MKTLweSampleFFT* ptr) {
    free(ptr);
}
EXPORT void free_MKTLweSampleFFT_array(int32_t nbelts, MKTLweSampleFFT* ptr){
    free(ptr);
}

// init
EXPORT void init_MKTLweSampleFFT(MKTLweSampleFFT* obj, const TLweParams* RLWEparams, const MKTFHEParams* MKparams, 
        LagrangeHalfCPolynomial *arr, double current_variance) 
{
    new(obj) MKTLweSampleFFT(RLWEparams, MKparams, arr, current_variance);
}
EXPORT void init_MKTLweSampleFFT_array(int32_t nbelts, MKTLweSampleFFT* obj, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams, LagrangeHalfCPolynomial *arr, double current_variance)
{
    for (int i = 0; i < nbelts; i++) {
        new(obj+i) MKTLweSampleFFT(RLWEparams, MKparams, arr, current_variance);
    }
}

// destroys the structure
EXPORT void destroy_MKTLweSampleFFT(MKTLweSampleFFT* obj) {
    obj->~MKTLweSampleFFT();
}
EXPORT void destroy_MKTLweSampleFFT_array(int32_t nbelts, MKTLweSampleFFT* obj) {
    for (int i = 0; i < nbelts; i++) {
        (obj+i)->~MKTLweSampleFFT();
    }
}

// new = alloc + init
EXPORT MKTLweSampleFFT* new_MKTLweSampleFFT(const TLweParams* RLWEparams, const MKTFHEParams* MKparams, 
        LagrangeHalfCPolynomial *arr, double current_variance) {
    return new MKTLweSampleFFT(RLWEparams, MKparams, arr, current_variance);
}
EXPORT MKTLweSampleFFT* new_MKTLweSampleFFT_array(int32_t nbelts, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams, LagrangeHalfCPolynomial *arr, double current_variance)
{
    MKTLweSampleFFT* obj = alloc_MKTLweSampleFFT_array(nbelts);
    init_MKTLweSampleFFT_array(nbelts,obj,RLWEparams,MKparams,arr,current_variance);
    return obj;
}

// delete = destroy + free
EXPORT void delete_MKTLweSampleFFT(MKTLweSampleFFT* obj) {
    delete obj;
}
EXPORT void delete_MKTLweSampleFFT_array(int32_t nbelts, MKTLweSampleFFT* obj) {
    destroy_MKTLweSampleFFT_array(nbelts,obj);
    free_MKTLweSampleFFT_array(nbelts,obj);
}







































/* **************************************************************************
***************************** VERSION 2 *************************************
************************************************************************** */



// MK RGSW UniEnc sample (d,F)=(d,f0,f1)
MKTGswUESample_v2::MKTGswUESample_v2(const TLweParams* RLWEparams, const MKTFHEParams* MKparams) :
        dg(MKparams->dg), N(RLWEparams->N)
{
    //a is a table of parties+1 polynomials, b is an alias for &a[parties]
    d = new_TorusPolynomial_array(3*dg, N);
    f0 = d + dg;
    f1 = d + 2*dg;
    party = 0; // party (from 0 to parties-1)
    current_variance = 0.0;
}

MKTGswUESample_v2::~MKTGswUESample_v2() {
    delete_TorusPolynomial_array(3*dg, d);
}


// alloc
EXPORT MKTGswUESample_v2* alloc_MKTGswUESample_v2() {
    return (MKTGswUESample_v2*) malloc(sizeof(MKTGswUESample_v2));
}
EXPORT MKTGswUESample_v2* alloc_MKTGswUESample_v2_array(int32_t nbelts) {
    return (MKTGswUESample_v2*) malloc(nbelts*sizeof(MKTGswUESample_v2));
}

//free memory space 
EXPORT void free_MKTGswUESample_v2(MKTGswUESample_v2* ptr) {
    free(ptr);
}
EXPORT void free_MKTGswUESample_v2_array(int32_t nbelts, MKTGswUESample_v2* ptr) {
    free(ptr);
}

// initialize the structure
EXPORT void init_MKTGswUESample_v2(MKTGswUESample_v2* obj, const TLweParams* RLWEparams, const MKTFHEParams* MKparams) {
    new(obj) MKTGswUESample_v2(RLWEparams, MKparams);
}
EXPORT void init_MKTGswUESample_v2_array(int32_t nbelts, MKTGswUESample_v2* obj, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams) 
{
    for (int i = 0; i < nbelts; i++) {
        new(obj+i) MKTGswUESample_v2(RLWEparams, MKparams);
    }
}

//destroys the structure
EXPORT void destroy_MKTGswUESample_v2(MKTGswUESample_v2* obj) {
    obj->~MKTGswUESample_v2();
}
EXPORT void destroy_MKTGswUESample_v2_array(int32_t nbelts, MKTGswUESample_v2* obj) {
    for (int i = 0; i < nbelts; i++) {
        (obj+i)->~MKTGswUESample_v2();
    }
}
 
// new = alloc + init
EXPORT MKTGswUESample_v2* new_MKTGswUESample_v2(const TLweParams* RLWEparams, const MKTFHEParams* MKparams) {
    return new MKTGswUESample_v2(RLWEparams, MKparams);
}
EXPORT MKTGswUESample_v2* new_MKTGswUESample_v2_array(int32_t nbelts, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams) 
{
    MKTGswUESample_v2* obj = alloc_MKTGswUESample_v2_array(nbelts);
    init_MKTGswUESample_v2_array(nbelts,obj,RLWEparams,MKparams);
    return obj;
}

// delete = destroy + free
EXPORT void delete_MKTGswUESample_v2(MKTGswUESample_v2* obj) {
    delete obj;
}
EXPORT void delete_MKTGswUESample_v2_array(int32_t nbelts, MKTGswUESample_v2* obj) {
    destroy_MKTGswUESample_v2_array(nbelts,obj);
    free_MKTGswUESample_v2_array(nbelts,obj);
}
























// MK RGSW UniEnc sample FFT (d,F)=(d,f0,f1)
MKTGswUESampleFFT_v2::MKTGswUESampleFFT_v2(const TLweParams* RLWEparams, const MKTFHEParams* MKparams, 
        LagrangeHalfCPolynomial *arr, int32_t party, double current_variance) : 
        dg(MKparams->dg)
{

    //a is a table of parties+1 polynomials, b is an alias for &a[parties]
    d = arr;
    f0 = d + dg;
    f1 = d + 2*dg;
    party = 0; // party (from 0 to parties-1)
    current_variance = 0.0;
}

MKTGswUESampleFFT_v2::~MKTGswUESampleFFT_v2() {
}


// alloc
EXPORT MKTGswUESampleFFT_v2* alloc_MKTGswUESampleFFT_v2() {
    return (MKTGswUESampleFFT_v2*) malloc(sizeof(MKTGswUESampleFFT_v2));
}
EXPORT MKTGswUESampleFFT_v2* alloc_MKTGswUESampleFFT_v2_array(int32_t nbelts) {
    return (MKTGswUESampleFFT_v2*) malloc(nbelts*sizeof(MKTGswUESampleFFT_v2));
}

//free memory space 
EXPORT void free_MKTGswUESampleFFT_v2(MKTGswUESampleFFT_v2* ptr) {
    free(ptr);
}
EXPORT void free_MKTGswUESampleFFT_v2_array(int32_t nbelts, MKTGswUESampleFFT_v2* ptr) {
    free(ptr);
}

// initialize the structure
EXPORT void init_MKTGswUESampleFFT_v2(MKTGswUESampleFFT_v2* obj, const TLweParams* RLWEparams, const MKTFHEParams* MKparams, 
        int32_t party, double current_variance) {
    LagrangeHalfCPolynomial *arr = new_LagrangeHalfCPolynomial_array(3*MKparams->dg, MKparams->N);
    new(obj) MKTGswUESampleFFT_v2(RLWEparams, MKparams, arr, party, current_variance);
}
EXPORT void init_MKTGswUESampleFFT_v2_array(int32_t nbelts, MKTGswUESampleFFT_v2* obj, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams, int32_t party, double current_variance) 
{
    LagrangeHalfCPolynomial *arr = new_LagrangeHalfCPolynomial_array(nbelts*3*MKparams->dg, MKparams->N);
    for (int i = 0; i < nbelts; i++) {
        new(obj+i) MKTGswUESampleFFT_v2(RLWEparams, MKparams, arr + i*3*MKparams->dg, party, current_variance);
    }
}

//destroys the structure
EXPORT void destroy_MKTGswUESampleFFT_v2(MKTGswUESampleFFT_v2* obj) {
    delete_LagrangeHalfCPolynomial_array(3*obj->dg, obj->d);
    obj->~MKTGswUESampleFFT_v2();
}
EXPORT void destroy_MKTGswUESampleFFT_v2_array(int32_t nbelts, MKTGswUESampleFFT_v2* obj) {
    delete_LagrangeHalfCPolynomial_array(nbelts*3*obj->dg, obj->d);
    for (int i = 0; i < nbelts; i++) {
        (obj+i)->~MKTGswUESampleFFT_v2();
    }
}
 
// new = alloc + init
EXPORT MKTGswUESampleFFT_v2* new_MKTGswUESampleFFT_v2(const TLweParams* RLWEparams, const MKTFHEParams* MKparams, 
        int32_t party, double current_variance) 
{
    MKTGswUESampleFFT_v2* obj = alloc_MKTGswUESampleFFT_v2();
    init_MKTGswUESampleFFT_v2(obj, RLWEparams, MKparams, party, current_variance);
    return obj;
    // return new MKTGswUESampleFFT_v2(RLWEparams, MKparams, arr, party, current_variance);
}
EXPORT MKTGswUESampleFFT_v2* new_MKTGswUESampleFFT_v2_array(int32_t nbelts, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams, int32_t party, double current_variance) 
{
    MKTGswUESampleFFT_v2* obj = alloc_MKTGswUESampleFFT_v2_array(nbelts);
    init_MKTGswUESampleFFT_v2_array(nbelts,obj,RLWEparams,MKparams, party, current_variance);
    return obj;
}

// delete = destroy + free
EXPORT void delete_MKTGswUESampleFFT_v2(MKTGswUESampleFFT_v2* obj) {
    destroy_MKTGswUESampleFFT_v2(obj);
    free_MKTGswUESampleFFT_v2(obj);
    //delete obj;
}
EXPORT void delete_MKTGswUESampleFFT_v2_array(int32_t nbelts, MKTGswUESampleFFT_v2* obj) {
    destroy_MKTGswUESampleFFT_v2_array(nbelts,obj);
    free_MKTGswUESampleFFT_v2_array(nbelts,obj);
}































// MK RGSW Expanded sample: party i D_i=(x_0, ..., x_{parties-1}, x_parties + d_i, y_0, ..., d_i+y_i, ..., y_perties, d_i)
MKTGswExpSample_v2::MKTGswExpSample_v2(const TLweParams* RLWEparams, const MKTFHEParams* MKparams) :
        parties(MKparams->parties), dg(MKparams->dg), N(RLWEparams->N)
{
    //a is a table of parties+1 polynomials, b is an alias for &a[parties]
    x = new_TorusPolynomial_array((2*(parties+1) + 1)*dg, N);
    y = x + (parties+1)*dg;
    d = y + (parties+1)*dg;
    party = 0; // party (from 0 to parties-1)
    current_variance = 0.0;
}

MKTGswExpSample_v2::~MKTGswExpSample_v2() {
    delete_TorusPolynomial_array((2*(parties+1)+1)*dg, x);
}



// alloc
EXPORT MKTGswExpSample_v2* alloc_MKTGswExpSample_v2() {
    return (MKTGswExpSample_v2*) malloc(sizeof(MKTGswExpSample_v2));
}
EXPORT MKTGswExpSample_v2* alloc_MKTGswExpSample_v2_array(int32_t nbelts) {
    return (MKTGswExpSample_v2*) malloc(nbelts*sizeof(MKTGswExpSample_v2));
}

//free memory space 
EXPORT void free_MKTGswExpSample_v2(MKTGswExpSample_v2* ptr) {
    free(ptr);
}
EXPORT void free_MKTGswExpSample_v2_array(int32_t nbelts, MKTGswExpSample_v2* ptr) {
    free(ptr);
}

// initialize the structure
EXPORT void init_MKTGswExpSample_v2(MKTGswExpSample_v2* obj, const TLweParams* RLWEparams, const MKTFHEParams* MKparams) {
    new(obj) MKTGswExpSample_v2(RLWEparams, MKparams);
}
EXPORT void init_MKTGswExpSample_v2_array(int32_t nbelts, MKTGswExpSample_v2* obj, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams) 
{
    for (int i = 0; i < nbelts; i++) {
        new(obj+i) MKTGswExpSample_v2(RLWEparams, MKparams);
    }
}

//destroys the structure
EXPORT void destroy_MKTGswExpSample_v2(MKTGswExpSample_v2* obj) {
    obj->~MKTGswExpSample_v2();
}
EXPORT void destroy_MKTGswExpSample_v2_array(int32_t nbelts, MKTGswExpSample_v2* obj) {
    for (int i = 0; i < nbelts; i++) {
        (obj+i)->~MKTGswExpSample_v2();
    }
}
 
// new = alloc + init
EXPORT MKTGswExpSample_v2* new_MKTGswExpSample_v2(const TLweParams* RLWEparams, const MKTFHEParams* MKparams) {
    return new MKTGswExpSample_v2(RLWEparams, MKparams);
}
EXPORT MKTGswExpSample_v2* new_MKTGswExpSample_v2_array(int32_t nbelts, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams) 
{
    MKTGswExpSample_v2* obj = alloc_MKTGswExpSample_v2_array(nbelts);
    init_MKTGswExpSample_v2_array(nbelts,obj,RLWEparams,MKparams);
    return obj;
}

// delete = destroy + free
EXPORT void delete_MKTGswExpSample_v2(MKTGswExpSample_v2* obj) {
    delete obj;
}
EXPORT void delete_MKTGswExpSample_v2_array(int32_t nbelts, MKTGswExpSample_v2* obj) {
    destroy_MKTGswExpSample_v2_array(nbelts,obj);
    free_MKTGswExpSample_v2_array(nbelts,obj);
}























// MK RGSW Expanded sample FFT: party i D_i=(x_0, ..., x_{parties-1}, x_parties + d_i, y_0, ..., d_i+y_i, ..., y_perties, d_i)
MKTGswExpSampleFFT_v2::MKTGswExpSampleFFT_v2(const TLweParams* RLWEparams, const MKTFHEParams* MKparams, 
        LagrangeHalfCPolynomial *arr, double current_variance) :
        parties(MKparams->parties), dg(MKparams->dg)
{

    //a is a table of parties+1 polynomials, b is an alias for &a[parties]
    x = arr;
    y = x + (parties+1)*dg;
    d = y + (parties+1)*dg;
    party = 0; // party (from 0 to parties-1)
    current_variance = 0.0;
}

MKTGswExpSampleFFT_v2::~MKTGswExpSampleFFT_v2() {
}



// alloc
EXPORT MKTGswExpSampleFFT_v2* alloc_MKTGswExpSampleFFT_v2() {
    return (MKTGswExpSampleFFT_v2*) malloc(sizeof(MKTGswExpSampleFFT_v2));
}
EXPORT MKTGswExpSampleFFT_v2* alloc_MKTGswExpSampleFFT_v2_array(int32_t nbelts) {
    return (MKTGswExpSampleFFT_v2*) malloc(nbelts*sizeof(MKTGswExpSampleFFT_v2));
}

//free memory space 
EXPORT void free_MKTGswExpSampleFFT_v2(MKTGswExpSampleFFT_v2* ptr) {
    free(ptr);
}
EXPORT void free_MKTGswExpSampleFFT_v2_array(int32_t nbelts, MKTGswExpSampleFFT_v2* ptr) {
    free(ptr);
}

// initialize the structure
EXPORT void init_MKTGswExpSampleFFT_v2(MKTGswExpSampleFFT_v2* obj, const TLweParams* RLWEparams, const MKTFHEParams* MKparams, 
        double current_variance) 
{
    LagrangeHalfCPolynomial *arr = new_LagrangeHalfCPolynomial_array((2*(1 + MKparams->parties) + 1)*MKparams->dg, MKparams->N);
    new(obj) MKTGswExpSampleFFT_v2(RLWEparams, MKparams, arr, current_variance);
}
EXPORT void init_MKTGswExpSampleFFT_v2_array(int32_t nbelts, MKTGswExpSampleFFT_v2* obj, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams, double current_variance) 
{
    LagrangeHalfCPolynomial *arr = new_LagrangeHalfCPolynomial_array(nbelts*(2*(1 + MKparams->parties) + 1)*MKparams->dg, MKparams->N);
    for (int i = 0; i < nbelts; i++) {
        new(obj+i) MKTGswExpSampleFFT_v2(RLWEparams, MKparams, arr + i*((2*(1 + MKparams->parties) + 1)*MKparams->dg), current_variance);
    }
}

//destroys the structure
EXPORT void destroy_MKTGswExpSampleFFT_v2(MKTGswExpSampleFFT_v2* obj) {
    delete_LagrangeHalfCPolynomial_array((2*(1 + obj->parties) + 1)*obj->dg, obj->x);
    obj->~MKTGswExpSampleFFT_v2();
}
EXPORT void destroy_MKTGswExpSampleFFT_v2_array(int32_t nbelts, MKTGswExpSampleFFT_v2* obj) {
    delete_LagrangeHalfCPolynomial_array(nbelts*(2*(1 + obj->parties) + 1)*obj->dg, obj->x);
    for (int i = 0; i < nbelts; i++) {
        (obj+i)->~MKTGswExpSampleFFT_v2();
    }
}
 
// new = alloc + init
EXPORT MKTGswExpSampleFFT_v2* new_MKTGswExpSampleFFT_v2(const TLweParams* RLWEparams, const MKTFHEParams* MKparams, 
        double current_variance) 
{
    MKTGswExpSampleFFT_v2* obj = alloc_MKTGswExpSampleFFT_v2();
    init_MKTGswExpSampleFFT_v2(obj, RLWEparams, MKparams, current_variance);
    return obj;
    // return new MKTGswExpSampleFFT_v2(RLWEparams, MKparams, current_variance);
}
EXPORT MKTGswExpSampleFFT_v2* new_MKTGswExpSampleFFT_v2_array(int32_t nbelts, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams, double current_variance) 
{
    MKTGswExpSampleFFT_v2* obj = alloc_MKTGswExpSampleFFT_v2_array(nbelts);
    init_MKTGswExpSampleFFT_v2_array(nbelts,obj,RLWEparams,MKparams, current_variance);
    return obj;
}

// delete = destroy + free
EXPORT void delete_MKTGswExpSampleFFT_v2(MKTGswExpSampleFFT_v2* obj) {
    destroy_MKTGswExpSampleFFT_v2(obj);
    free_MKTGswExpSampleFFT_v2(obj);
    //delete obj;
}
EXPORT void delete_MKTGswExpSampleFFT_v2_array(int32_t nbelts, MKTGswExpSampleFFT_v2* obj) {
    destroy_MKTGswExpSampleFFT_v2_array(nbelts,obj);
    free_MKTGswExpSampleFFT_v2_array(nbelts,obj);
}








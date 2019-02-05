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



























// MK RGSW UniEnc sample (C,D,F)=(c0,c1,d0,d1,f0,f1)
MKTGswUESample::MKTGswUESample(const TLweParams* RLWEparams, const MKTFHEParams* MKparams) :
		dg(MKparams->dg), N(RLWEparams->N)
{
	//a is a table of parties+1 polynomials, b is an alias for &a[parties]
	c = new_TorusPolynomial_array(6*dg, N);
    d = c + 2*dg;
    f = c + 4*dg;
    party = 0; // party (from 0 to parties-1)
    current_variance = 0.0;
}

MKTGswUESample::~MKTGswUESample() {
    delete_TorusPolynomial_array(6*dg, c);
}


// alloc
EXPORT MKTGswUESample* alloc_MKTGswUESample() {
    return (MKTGswUESample*) malloc(sizeof(MKTGswUESample));
}
EXPORT MKTGswUESample* alloc_MKTGswUESample_array(int32_t nbelts) {
    return (MKTGswUESample*) malloc(nbelts*sizeof(MKTGswUESample));
}

//free memory space 
EXPORT void free_MKTGswUESample(MKTGswUESample* ptr) {
    free(ptr);
}
EXPORT void free_MKTGswUESample_array(int32_t nbelts, MKTGswUESample* ptr) {
    free(ptr);
}

// initialize the structure
EXPORT void init_MKTGswUESample(MKTGswUESample* obj, const TLweParams* RLWEparams, const MKTFHEParams* MKparams) {
    new(obj) MKTGswUESample(RLWEparams, MKparams);
}
EXPORT void init_MKTGswUESample_array(int32_t nbelts, MKTGswUESample* obj, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams) 
{
    for (int i = 0; i < nbelts; i++) {
        new(obj+i) MKTGswUESample(RLWEparams, MKparams);
    }
}

//destroys the structure
EXPORT void destroy_MKTGswUESample(MKTGswUESample* obj) {
    obj->~MKTGswUESample();
}
EXPORT void destroy_MKTGswUESample_array(int32_t nbelts, MKTGswUESample* obj) {
    for (int i = 0; i < nbelts; i++) {
        (obj+i)->~MKTGswUESample();
    }
}
 
// new = alloc + init
EXPORT MKTGswUESample* new_MKTGswUESample(const TLweParams* RLWEparams, const MKTFHEParams* MKparams) {
    return new MKTGswUESample(RLWEparams, MKparams);
}
EXPORT MKTGswUESample* new_MKTGswUESample_array(int32_t nbelts, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams) 
{
    MKTGswUESample* obj = alloc_MKTGswUESample_array(nbelts);
    init_MKTGswUESample_array(nbelts,obj,RLWEparams,MKparams);
    return obj;
}

// delete = destroy + free
EXPORT void delete_MKTGswUESample(MKTGswUESample* obj) {
    delete obj;
}
EXPORT void delete_MKTGswUESample_array(int32_t nbelts, MKTGswUESample* obj) {
    destroy_MKTGswUESample_array(nbelts,obj);
    free_MKTGswUESample_array(nbelts,obj);
}
























// MK RGSW UniEnc sample FFT (C,D,F)=(c0,c1,d0,d1,f0,f1)
MKTGswUESampleFFT::MKTGswUESampleFFT(const TLweParams* RLWEparams, const MKTFHEParams* MKparams, 
		LagrangeHalfCPolynomial *arr, double current_variance) : 
		dg(MKparams->dg)
{

	//a is a table of parties+1 polynomials, b is an alias for &a[parties]
    c = arr;
    d = c + 2*dg;
    f = c + 4*dg;
    party = 0; // party (from 0 to parties-1)
    current_variance = 0.0;
}

MKTGswUESampleFFT::~MKTGswUESampleFFT() {
}


// alloc
EXPORT MKTGswUESampleFFT* alloc_MKTGswUESampleFFT() {
    return (MKTGswUESampleFFT*) malloc(sizeof(MKTGswUESampleFFT));
}
EXPORT MKTGswUESampleFFT* alloc_MKTGswUESampleFFT_array(int32_t nbelts) {
    return (MKTGswUESampleFFT*) malloc(nbelts*sizeof(MKTGswUESampleFFT));
}

//free memory space 
EXPORT void free_MKTGswUESampleFFT(MKTGswUESampleFFT* ptr) {
    free(ptr);
}
EXPORT void free_MKTGswUESampleFFT_array(int32_t nbelts, MKTGswUESampleFFT* ptr) {
    free(ptr);
}

// initialize the structure
EXPORT void init_MKTGswUESampleFFT(MKTGswUESampleFFT* obj, const TLweParams* RLWEparams, const MKTFHEParams* MKparams, 
        LagrangeHalfCPolynomial *arr, double current_variance) {
    new(obj) MKTGswUESampleFFT(RLWEparams, MKparams, arr, current_variance);
}
EXPORT void init_MKTGswUESampleFFT_array(int32_t nbelts, MKTGswUESampleFFT* obj, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams, LagrangeHalfCPolynomial *arr, double current_variance) 
{
    for (int i = 0; i < nbelts; i++) {
        new(obj+i) MKTGswUESampleFFT(RLWEparams, MKparams, arr, current_variance);
    }
}

//destroys the structure
EXPORT void destroy_MKTGswUESampleFFT(MKTGswUESampleFFT* obj) {
    obj->~MKTGswUESampleFFT();
}
EXPORT void destroy_MKTGswUESampleFFT_array(int32_t nbelts, MKTGswUESampleFFT* obj) {
    for (int i = 0; i < nbelts; i++) {
        (obj+i)->~MKTGswUESampleFFT();
    }
}
 
// new = alloc + init
EXPORT MKTGswUESampleFFT* new_MKTGswUESampleFFT(const TLweParams* RLWEparams, const MKTFHEParams* MKparams, 
        LagrangeHalfCPolynomial *arr, double current_variance) {
    return new MKTGswUESampleFFT(RLWEparams, MKparams, arr, current_variance);
}
EXPORT MKTGswUESampleFFT* new_MKTGswUESampleFFT_array(int32_t nbelts, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams, LagrangeHalfCPolynomial *arr, double current_variance) 
{
    MKTGswUESampleFFT* obj = alloc_MKTGswUESampleFFT_array(nbelts);
    init_MKTGswUESampleFFT_array(nbelts,obj,RLWEparams,MKparams, arr, current_variance);
    return obj;
}

// delete = destroy + free
EXPORT void delete_MKTGswUESampleFFT(MKTGswUESampleFFT* obj) {
    delete obj;
}
EXPORT void delete_MKTGswUESampleFFT_array(int32_t nbelts, MKTGswUESampleFFT* obj) {
    destroy_MKTGswUESampleFFT_array(nbelts,obj);
    free_MKTGswUESampleFFT_array(nbelts,obj);
}



































// MK RGSW Expanded sample: party i C=(y1, ..., d1, ..., yk, c1, d0+x1, ..., d0, ..., d0+xk, c0)
MKTGswExpSample::MKTGswExpSample(const TLweParams* RLWEparams, const MKTFHEParams* MKparams) :
		parties(MKparams->parties), dg(MKparams->dg), N(RLWEparams->N)
{
	//a is a table of parties+1 polynomials, b is an alias for &a[parties]
	c = new_TorusPolynomial_array(2*(parties+1)*dg, N);
    y = c;
    c1 = y + parties*dg;
    x = c + (parties+1)*dg;
    c0 = x + parties*dg;
    party = 0; // party (from 0 to parties-1)
    current_variance = 0.0;
}

MKTGswExpSample::~MKTGswExpSample() {
    delete_TorusPolynomial_array(2*(parties+1)*dg, c);
}



// alloc
EXPORT MKTGswExpSample* alloc_MKTGswExpSample() {
    return (MKTGswExpSample*) malloc(sizeof(MKTGswExpSample));
}
EXPORT MKTGswExpSample* alloc_MKTGswExpSample_array(int32_t nbelts) {
    return (MKTGswExpSample*) malloc(nbelts*sizeof(MKTGswExpSample));
}

//free memory space 
EXPORT void free_MKTGswExpSample(MKTGswExpSample* ptr) {
    free(ptr);
}
EXPORT void free_MKTGswExpSample_array(int32_t nbelts, MKTGswExpSample* ptr) {
    free(ptr);
}

// initialize the structure
EXPORT void init_MKTGswExpSample(MKTGswExpSample* obj, const TLweParams* RLWEparams, const MKTFHEParams* MKparams) {
    new(obj) MKTGswExpSample(RLWEparams, MKparams);
}
EXPORT void init_MKTGswExpSample_array(int32_t nbelts, MKTGswExpSample* obj, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams) 
{
    for (int i = 0; i < nbelts; i++) {
        new(obj+i) MKTGswExpSample(RLWEparams, MKparams);
    }
}

//destroys the structure
EXPORT void destroy_MKTGswExpSample(MKTGswExpSample* obj) {
    obj->~MKTGswExpSample();
}
EXPORT void destroy_MKTGswExpSample_array(int32_t nbelts, MKTGswExpSample* obj) {
    for (int i = 0; i < nbelts; i++) {
        (obj+i)->~MKTGswExpSample();
    }
}
 
// new = alloc + init
EXPORT MKTGswExpSample* new_MKTGswExpSample(const TLweParams* RLWEparams, const MKTFHEParams* MKparams) {
    return new MKTGswExpSample(RLWEparams, MKparams);
}
EXPORT MKTGswExpSample* new_MKTGswExpSample_array(int32_t nbelts, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams) 
{
    MKTGswExpSample* obj = alloc_MKTGswExpSample_array(nbelts);
    init_MKTGswExpSample_array(nbelts,obj,RLWEparams,MKparams);
    return obj;
}

// delete = destroy + free
EXPORT void delete_MKTGswExpSample(MKTGswExpSample* obj) {
    delete obj;
}
EXPORT void delete_MKTGswExpSample_array(int32_t nbelts, MKTGswExpSample* obj) {
    destroy_MKTGswExpSample_array(nbelts,obj);
    free_MKTGswExpSample_array(nbelts,obj);
}























// MK RGSW Expanded sample FFT: party i C=(y1, ..., d1, ..., yk, c1, d0+x1, ..., d0, ..., d0+xk, c0)
MKTGswExpSampleFFT::MKTGswExpSampleFFT(const TLweParams* RLWEparams, const MKTFHEParams* MKparams, 
		LagrangeHalfCPolynomial *arr, double current_variance) :
        parties(MKparams->parties), dg(MKparams->dg)
{

	//a is a table of parties+1 polynomials, b is an alias for &a[parties]
    c = arr;
    y = c;
    c1 = y + parties*dg;
    x = c + (parties+1)*dg;
    c0 = x + parties*dg;
    party = 0; // party (from 0 to parties-1)
    current_variance = 0.0;
}

MKTGswExpSampleFFT::~MKTGswExpSampleFFT() {
}



// alloc
EXPORT MKTGswExpSampleFFT* alloc_MKTGswExpSampleFFT() {
    return (MKTGswExpSampleFFT*) malloc(sizeof(MKTGswExpSampleFFT));
}
EXPORT MKTGswExpSampleFFT* alloc_MKTGswExpSampleFFT_array(int32_t nbelts) {
    return (MKTGswExpSampleFFT*) malloc(nbelts*sizeof(MKTGswExpSampleFFT));
}

//free memory space 
EXPORT void free_MKTGswExpSampleFFT(MKTGswExpSampleFFT* ptr) {
    free(ptr);
}
EXPORT void free_MKTGswExpSampleFFT_array(int32_t nbelts, MKTGswExpSampleFFT* ptr) {
    free(ptr);
}

// initialize the structure
EXPORT void init_MKTGswExpSampleFFT(MKTGswExpSampleFFT* obj, const TLweParams* RLWEparams, const MKTFHEParams* MKparams, 
        double current_variance) 
{
    LagrangeHalfCPolynomial *arr = new_LagrangeHalfCPolynomial_array(2*(1 + MKparams->parties)*MKparams->dg, MKparams->N);
    new(obj) MKTGswExpSampleFFT(RLWEparams, MKparams, arr, current_variance);
}
EXPORT void init_MKTGswExpSampleFFT_array(int32_t nbelts, MKTGswExpSampleFFT* obj, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams, double current_variance) 
{
    LagrangeHalfCPolynomial *arr = new_LagrangeHalfCPolynomial_array(nbelts*2*(1 + MKparams->parties)*MKparams->dg, MKparams->N);
    for (int i = 0; i < nbelts; i++) {
        new(obj+i) MKTGswExpSampleFFT(RLWEparams, MKparams, arr + i*(2*(1 + MKparams->parties)*MKparams->dg), current_variance);
    }
}

//destroys the structure
EXPORT void destroy_MKTGswExpSampleFFT(MKTGswExpSampleFFT* obj) {
    delete_LagrangeHalfCPolynomial_array(2*(1 + obj->parties)*obj->dg, obj->c);
    obj->~MKTGswExpSampleFFT();
}
EXPORT void destroy_MKTGswExpSampleFFT_array(int32_t nbelts, MKTGswExpSampleFFT* obj) {
    delete_LagrangeHalfCPolynomial_array(nbelts*2*(1 + obj->parties)*obj->dg, obj->c);
    for (int i = 0; i < nbelts; i++) {
        (obj+i)->~MKTGswExpSampleFFT();
    }
}
 
// new = alloc + init
EXPORT MKTGswExpSampleFFT* new_MKTGswExpSampleFFT(const TLweParams* RLWEparams, const MKTFHEParams* MKparams, 
        double current_variance) 
{
    MKTGswExpSampleFFT* obj = alloc_MKTGswExpSampleFFT();
    init_MKTGswExpSampleFFT(obj, RLWEparams, MKparams, current_variance);
    return obj;
    // return new MKTGswExpSampleFFT(RLWEparams, MKparams, current_variance);
}
EXPORT MKTGswExpSampleFFT* new_MKTGswExpSampleFFT_array(int32_t nbelts, const TLweParams* RLWEparams, 
        const MKTFHEParams* MKparams, double current_variance) 
{
    MKTGswExpSampleFFT* obj = alloc_MKTGswExpSampleFFT_array(nbelts);
    init_MKTGswExpSampleFFT_array(nbelts,obj,RLWEparams,MKparams, current_variance);
    return obj;
}

// delete = destroy + free
EXPORT void delete_MKTGswExpSampleFFT(MKTGswExpSampleFFT* obj) {
    destroy_MKTGswExpSampleFFT(obj);
    free_MKTGswExpSampleFFT(obj);
    //delete obj;
}
EXPORT void delete_MKTGswExpSampleFFT_array(int32_t nbelts, MKTGswExpSampleFFT* obj) {
    destroy_MKTGswExpSampleFFT_array(nbelts,obj);
    free_MKTGswExpSampleFFT_array(nbelts,obj);
}






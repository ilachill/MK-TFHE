#include <cstdlib>
#include <new>
#include <stdlib.h>

#include "lwekey.h"
#include "lweparams.h"

using namespace std;

LweKey::LweKey(const LweParams* params) {
    this->params = params;
    this->key = new int32_t[params->n];
}

LweKey::~LweKey() {
    delete[] key;
}








// allocate memory space 
EXPORT LweKey* alloc_LweKey() {
    return (LweKey*) malloc(sizeof(LweKey));
}
EXPORT LweKey* alloc_LweKey_array(int32_t nbelts) {
    return (LweKey*) malloc(nbelts*sizeof(LweKey));
}

// free memory space 
EXPORT void free_LweKey(LweKey* ptr) {
    free(ptr);
}
EXPORT void free_LweKey_array(int32_t nbelts, LweKey* ptr) {
    free(ptr);
}

// initialize the structure
EXPORT void init_LweKey(LweKey* obj, const LweParams* params) {
    new(obj) LweKey(params);
}
EXPORT void init_LweKey_array(int32_t nbelts, LweKey* obj, const LweParams* params)
    {
    for (int i = 0; i < nbelts; i++) {
        new(obj+i) LweKey(params);
    }
}

// destroys the structure
EXPORT void destroy_LweKey(LweKey* obj) {
    obj->~LweKey();
}
EXPORT void destroy_LweKey_array(int32_t nbelts, LweKey* obj) {
    for (int i = 0; i < nbelts; i++) {
        (obj+i)->~LweKey();
    }
}
 
// new = alloc + init
EXPORT LweKey* new_LweKey(const LweParams* params) {
    return new LweKey(params);
}
EXPORT LweKey* new_LweKey_array(int32_t nbelts, const LweParams* params) {
    LweKey* obj = alloc_LweKey_array(nbelts);
    init_LweKey_array(nbelts,obj,params);
    return obj;
}

// delete = destroy + free
EXPORT void delete_LweKey(LweKey* obj) {
    delete obj;
}
EXPORT void delete_LweKey_array(int32_t nbelts, LweKey* obj) {
    destroy_LweKey_array(nbelts,obj);
    free_LweKey_array(nbelts,obj);
}




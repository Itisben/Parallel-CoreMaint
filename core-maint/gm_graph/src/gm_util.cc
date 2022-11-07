#include "gm_util.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>

/*
 * This is an implementation of the Tokenizer class.
 * This implementation uses the methods in the C++ string class;
 * specifically, find_first_not_of, find_first_of, and substr.
 */
void GM_Tokenizer::reset() {
    current_pos_ = 0;
    next_token_start_ = 0;
    found_next_token_start_ = false;
    number_of_tokens_ = 0;
}

void GM_Tokenizer::findNextTokenStart() {
    if (!found_next_token_start_) {
        next_token_start_ = str_.find_first_not_of(delim_, current_pos_);
        found_next_token_start_ = true;
    }
}

bool GM_Tokenizer::isEscaped(size_t pos) {
    //we need this function to determine if a character is escaped
    //if there is an odd number of \ in front of the character -> escaped
    //even number -> not escaped (just a bunch of \ in front of it)
    int count = 0;
    while (pos > 0) {
        pos--;
        if (str_[pos] == '\\')
            count++;
        else
            break;
    }
    return count % 2 == 1;
}

std::string GM_Tokenizer::getNextToken() {
    if (!found_next_token_start_) {
        findNextTokenStart();
    }
    size_t next_token_end = str_.find_first_of(delim_, next_token_start_);
    current_pos_ = next_token_end; // setting beyond the next token
    found_next_token_start_ = false;
    if (str_[next_token_start_] == '"') {
//        printf("A: %s\n", str_.substr(next_token_start_, next_token_end - next_token_start_).c_str());
        /*token is string - move next_token_end to end of string -> first occurrence of " without \ */
        next_token_end = str_.find_first_of("\"", next_token_start_ + 1);
//        printf("C: %s\n", str_.substr(next_token_start_, next_token_end - next_token_start_).c_str());
        assert(next_token_end < str_.length());
        while (isEscaped(next_token_end)) {
            next_token_end = str_.find_first_of("\"", next_token_end + 1);
//            printf("D: %s\n", str_.substr(next_token_start_, next_token_end - next_token_start_).c_str());
            assert(next_token_end < str_.length());
        }
        next_token_end++;
        assert(next_token_end <= str_.length());
//        printf("E: %s\n", str_.substr(next_token_start_, next_token_end - next_token_start_).c_str());
    }
//    printf("token: %s\n", str_.substr(next_token_start_, next_token_end - next_token_start_).c_str());
    return str_.substr(next_token_start_, next_token_end - next_token_start_);
}

bool GM_Tokenizer::hasNextToken() {
    if (!found_next_token_start_) {
        findNextTokenStart();
    }
    return next_token_start_ != std::string::npos;
}

long GM_Tokenizer::countNumberOfTokens() {
    if (number_of_tokens_ == 0) {
        for (size_t pos = 0; pos != std::string::npos; number_of_tokens_++) {
            pos = str_.find_first_not_of(delim_, pos);
            if (pos == std::string::npos) break;
            pos = str_.find_first_of(delim_, pos);
        }
    }
    return number_of_tokens_;
}

int gmutil_parseIntFromString(std::string s) {
    return atoi(s.c_str());
}

long gmutil_parseLongFromString(std::string s) {
    return atol(s.c_str());
}

bool gmutil_parseBoolFromString(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s == "1" || s == "true";
}

float gmutil_parseFloatFromString(std::string s) {
    return strtof(s.c_str(), NULL);
}

double gmutil_parseDoubleFromString(std::string s) {
    return strtod(s.c_str(), NULL);
}

#ifdef GM_NODE64
node_t gmutil_parseNodeFromString(std::string s) {
    return (node_t) gmutil_parseLongFromString(s);
}
#else
node_t gmutil_parseNodeFromString(std::string s) {
    return (node_t) gmutil_parseIntFromString(s);
}
#endif

#ifdef GM_EDGE64
edge_t gmutil_parseEdgeFromString(std::string s) {
    return (edge_t) gmutil_parseLongFromString(s);
}
#else
edge_t gmutil_parseEdgeFromString(std::string s) {
    return (edge_t) gmutil_parseIntFromString(s);
}
#endif

std::string gmutil_parseStringFromString(const std::string s) {
    //make sure that you checked that the input is ok before with 'gmutil_check_input_string'
    if (s.length() == 0) return std::string("");

    if (s[0] == '"') {
        assert(s.length() >= 2);
        //remove leading and tailing '"'
        return gmutil_parseStringFromString(s.substr(0, s.length() - 1).substr(1, s.length() - 2));
    }

    std::string result = s;
    for (size_t i = 0; i < result.length(); i++) {
        if (result[i] == '\\') {
            result.erase(i, 1);
        }
    }
    return result;
}


/*
 * Method to create arrays based on a given value type and size.
 */
void *gmutil_getArrayType(VALUE_TYPE vt, int size) {
    switch (vt) {
        case GMTYPE_BOOL:
            return (void *) new bool[size];
        case GMTYPE_INT:
            return (void *) new int[size];
        case GMTYPE_LONG:
            return (void *) new long[size];
        case GMTYPE_FLOAT:
            return (void *) new float[size];
        case GMTYPE_DOUBLE:
            return (void *) new double[size];
        case GMTYPE_NODE:
            return (void *) new node_t[size];
        case GMTYPE_EDGE:
            return (void *) new edge_t[size];
        case GMTYPE_END:
        default:
            assert(false);
            return NULL; // Control should never reach this case.
    }
    return NULL;
}

int gmutil_getSizeOfType(VALUE_TYPE vt) {
    switch (vt) {
        case GMTYPE_BOOL:
            return 1;
        case GMTYPE_INT:
            return 4;
        case GMTYPE_LONG:
            return 8;
        case GMTYPE_FLOAT:
            return 4;
        case GMTYPE_DOUBLE:
            return 8;
        case GMTYPE_NODE:
            return sizeof(node_t);
        case GMTYPE_EDGE:
            return sizeof(edge_t);
        case GMTYPE_END:
        default:
            assert(false);
            return -1; // Control should never reach this case.
    }
    return -1;
}

void* gmutil_createVectorType(VALUE_TYPE vt) {
    switch (vt) {
        case GMTYPE_BOOL:
            return (void *) new GM_BVECT();
        case GMTYPE_INT:
            return (void *) new GM_IVECT();
        case GMTYPE_LONG:
            return (void *) new GM_LVECT();
        case GMTYPE_FLOAT:
            return (void *) new GM_FVECT();
        case GMTYPE_DOUBLE:
            return (void *) new GM_DVECT();
        case GMTYPE_NODE:
            return (void *) new GM_NVECT();
        case GMTYPE_EDGE:
            return (void *) new GM_EVECT();
        case GMTYPE_END:
        default:
            assert(false);
            return NULL; // Control should never reach this case.
    }
}
void gmutil_deleteVectorType(void* vector, VALUE_TYPE vt) {
    switch (vt) {
        case GMTYPE_BOOL:
            delete ((GM_BVECT*) vector);
            break;
        case GMTYPE_INT:
            delete ((GM_IVECT*) vector);
            break;
        case GMTYPE_LONG:
            delete ((GM_LVECT*) vector);
            break;
        case GMTYPE_FLOAT:
            delete ((GM_FVECT*) vector);
            break;
        case GMTYPE_DOUBLE:
            delete ((GM_DVECT*) vector);
            break;
        case GMTYPE_NODE:
            delete ((GM_NVECT*) vector);
            break;
        case GMTYPE_EDGE:
            delete ((GM_EVECT*) vector);
            break;
            break;
        case GMTYPE_END:
        default:
            assert(false);
            break; // Control should never reach this case.
    }
    return;
}

void gmutil_loadValueIntoVector(void *vector, std::string val, VALUE_TYPE vt) {
    switch (vt) {
        case GMTYPE_BOOL: {
            bool b = gmutil_parseBoolFromString(val);
            GM_BVECT* V = (GM_BVECT*) vector;
            V->push_back(b);
            break;
        }

        case GMTYPE_INT: {
            int i = gmutil_parseIntFromString(val);
            GM_IVECT* V = (GM_IVECT*) vector;
            V->push_back(i);
            break;
        }

        case GMTYPE_LONG: {
            long l = gmutil_parseLongFromString(val);
            GM_LVECT* V = (GM_LVECT*) vector;
            V->push_back(l);
            break;
        }

        case GMTYPE_FLOAT: {
            float f = gmutil_parseFloatFromString(val);
            GM_FVECT* V = (GM_FVECT*) vector;
            V->push_back(f);
            break;
        }

        case GMTYPE_DOUBLE: {
            double d = gmutil_parseDoubleFromString(val);
            GM_DVECT* V = (GM_DVECT*) vector;
            V->push_back(d);
            break;
        }

        case GMTYPE_NODE: {
            node_t n = gmutil_parseNodeFromString(val);
            GM_NVECT* V = (GM_NVECT*) vector;
            V->push_back(n);
            break;
        }

        case GMTYPE_EDGE: {
            edge_t e = gmutil_parseEdgeFromString(val);
            GM_EVECT* V = (GM_EVECT*) vector;
            V->push_back(e);
            break;
        }


        case GMTYPE_END:
        default:
            assert(false);
            return; // Control should never reach this case.
    }
}

void gmutil_loadValueIntoVectorAtPosition(void *vector, std::string val, VALUE_TYPE vt, size_t position) {
    switch (vt) {
        case GMTYPE_BOOL: {
            bool b = gmutil_parseBoolFromString(val);
            GM_BVECT* V = (GM_BVECT*) vector;
            if (V->size() <= position) V->resize(position + 1, false);
            V->at(position) = b;
            break;
        }

        case GMTYPE_INT: {
            int i = gmutil_parseIntFromString(val);
            GM_IVECT* V = (GM_IVECT*) vector;
            if (V->size() <= position) V->resize(position + 1, 0);
            V->at(position) = i;
            break;
        }

        case GMTYPE_LONG: {
            long l = gmutil_parseLongFromString(val);
            GM_LVECT* V = (GM_LVECT*) vector;
            if (V->size() <= position) V->resize(position + 1, 0);
            V->at(position) = l;
            break;
        }

        case GMTYPE_FLOAT: {
            float f = gmutil_parseDoubleFromString(val);
            GM_FVECT* V = (GM_FVECT*) vector;
            if (V->size() <= position) V->resize(position + 1, 0.0);
            V->at(position) = f;
            break;
        }

        case GMTYPE_DOUBLE: {
            double d = gmutil_parseDoubleFromString(val);
            GM_DVECT* V = (GM_DVECT*) vector;
            if (V->size() <= position) V->resize(position + 1, 0.0);
            V->at(position) = d;
            break;
        }

        case GMTYPE_NODE: {
            node_t n = gmutil_parseNodeFromString(val);
            GM_NVECT* V = (GM_NVECT*) vector;
            if (V->size() <= position) V->resize(position + 1, -1);
            V->at(position) = n;
            break;
        }

        case GMTYPE_EDGE: {
            edge_t e = gmutil_parseEdgeFromString(val);
            GM_EVECT* V = (GM_EVECT*) vector;
            if (V->size() <= position) V->resize(position + 1, -1);
            V->at(position) = e;
            break;
        }

        case GMTYPE_END:
        default:
            assert(false);
            return; // Control should never reach this case.
    }
}

void gmutil_loadDummyValueIntoVector(void *vector, VALUE_TYPE vt) {
    switch (vt) {
        case GMTYPE_BOOL: {
            GM_BVECT* V = (GM_BVECT*) vector;
            V->push_back(false);
            break;
        }
        case GMTYPE_INT: {
            GM_IVECT* V = (GM_IVECT*) vector;
            V->push_back(0);
            break;
        }
        case GMTYPE_LONG: {
            GM_LVECT* V = (GM_LVECT*) vector;
            V->push_back(0);
            break;
        }
        case GMTYPE_FLOAT: {
            GM_FVECT* V = (GM_FVECT*) vector;
            V->push_back(0);
            break;
        }
        case GMTYPE_DOUBLE: {
            GM_DVECT* V = (GM_DVECT*) vector;
            V->push_back(0);
            break;
        }
        case GMTYPE_NODE: {
            GM_NVECT* V = (GM_NVECT*) vector;
            V->push_back(0);
            break;
        }
        case GMTYPE_EDGE: {
            GM_EVECT* V = (GM_EVECT*) vector;
            V->push_back(0);
            break;
        }
        case GMTYPE_END:
        default:
            assert(false);
            return; // Control should never reach this case.
    }
}

static void gmutil_copyVectorIntoArray(GM_BVECT& SRC, bool* dest, edge_t* ind) {
#pragma omp parallel for
    //for(size_t i = 0; i < SRC.size(); i++) dest[(ind==NULL)?i:ind[i]] = SRC[i];
    for (size_t i = 0; i < SRC.size(); i++)
        dest[i] = SRC[(ind == NULL) ? i : ind[i]];
}

static void gmutil_copyVectorIntoArray(GM_IVECT& SRC, int* dest, edge_t* ind) {
#pragma omp parallel for
    //for(size_t i = 0; i < SRC.size(); i++) dest[(ind==NULL)?i:ind[i]] = SRC [i];
    for (size_t i = 0; i < SRC.size(); i++)
        dest[i] = SRC[(ind == NULL) ? i : ind[i]];
}

static void gmutil_copyVectorIntoArray(GM_LVECT& SRC, long* dest, edge_t* ind) {
#pragma omp parallel for
    //for(size_t i = 0; i < SRC.size(); i++) dest[(ind==NULL)?i:ind[i]] = SRC [i];
    for (size_t i = 0; i < SRC.size(); i++)
        dest[i] = SRC[(ind == NULL) ? i : ind[i]];
}

static void gmutil_copyVectorIntoArray(GM_FVECT& SRC, float* dest, edge_t* ind) {
#pragma omp parallel for
    //for(size_t i = 0; i < SRC.size(); i++) dest[(ind==NULL)?i:ind[i]] = SRC [i];
    for (size_t i = 0; i < SRC.size(); i++)
        dest[i] = SRC[(ind == NULL) ? i : ind[i]];
}

static void gmutil_copyVectorIntoArray(GM_DVECT& SRC, double* dest, edge_t* ind) {
#pragma omp parallel for
    //for(size_t i = 0; i < SRC.size(); i++) dest[(ind==NULL)?i:ind[i]] = SRC [i];
    for (size_t i = 0; i < SRC.size(); i++)
        dest[i] = SRC[(ind == NULL) ? i : ind[i]];
}

void gmutil_copyVectorIntoArray(void* vector, void* array, VALUE_TYPE vt, edge_t* indirection) {
    // should test copy vs array pointer
    switch (vt) {
        case GMTYPE_BOOL:
            gmutil_copyVectorIntoArray(*((GM_BVECT*) vector), (bool*) array, indirection);
            break;
        case GMTYPE_INT:
            gmutil_copyVectorIntoArray(*((GM_IVECT*) vector), (int*) array, indirection);
            break;
        case GMTYPE_LONG:
            gmutil_copyVectorIntoArray(*((GM_LVECT*) vector), (long*) array, indirection);
            break;
        case GMTYPE_FLOAT:
            gmutil_copyVectorIntoArray(*((GM_FVECT*) vector), (float*) array, indirection);
            break;
        case GMTYPE_DOUBLE:
            gmutil_copyVectorIntoArray(*((GM_DVECT*) vector), (double*) array, indirection);
            break;
        case GMTYPE_NODE:
            gmutil_copyVectorIntoArray(*((GM_NVECT*) vector), (node_t*) array, indirection);
            break;
        case GMTYPE_EDGE:
            gmutil_copyVectorIntoArray(*((GM_EVECT*) vector), (edge_t*) array, indirection);
            break;

        case GMTYPE_END:
        default:
            assert(false);
            return; // Control should never reach this case.
    }
}

/*
 * Method to load a value from the given string 'val'
 * and write to a given location in an array based on the given value type.
 */
void loadValueBasedOnType(void *arr, long pos, std::string val, VALUE_TYPE vt) {
    switch (vt) {
        case GMTYPE_BOOL:
            ((bool *) arr)[pos] = gmutil_parseBoolFromString(val);
            break;
        case GMTYPE_INT:
            ((int *) arr)[pos] = gmutil_parseIntFromString(val);
            break;
        case GMTYPE_LONG:
            ((long *) arr)[pos] = gmutil_parseLongFromString(val);
            break;
        case GMTYPE_FLOAT:
            ((float *) arr)[pos] = gmutil_parseFloatFromString(val);
            break;
        case GMTYPE_DOUBLE:
            ((double *) arr)[pos] = gmutil_parseDoubleFromString(val);
            break;
        case GMTYPE_NODE:
            ((node_t *) arr)[pos] = gmutil_parseNodeFromString(val);
            break;
        case GMTYPE_EDGE:
            ((edge_t *) arr)[pos] = gmutil_parseEdgeFromString(val);
            break;
        case GMTYPE_END:
        default:
            assert(false);
            return; // Control should never reach this case.
    }
}

const char* gmutil_getTypeString(VALUE_TYPE v) {
    return (v == GMTYPE_BOOL) ? "bool" : (v == GMTYPE_INT) ? "int" : (v == GMTYPE_LONG) ? "long" : (v == GMTYPE_FLOAT) ? "float" :
           (v == GMTYPE_DOUBLE) ? "double" : (v == GMTYPE_EDGE) ? "edge" : (v == GMTYPE_NODE) ? "node" : "??";
}

/*
 * Method to read a value from the given location in an array based on the given value type
 * and store it in the file
 */
void storeValueBasedOnType(void *arr, long pos, std::ofstream& file, VALUE_TYPE vt) {
    switch (vt) {
        case GMTYPE_BOOL:
            file << std::boolalpha << ((bool *) arr)[pos];
            break;
        case GMTYPE_INT:
            file << ((int *) arr)[pos];
            break;
        case GMTYPE_LONG:
            file << ((long *) arr)[pos];
            break;
        case GMTYPE_FLOAT:
            file << ((float *) arr)[pos];
            break;
        case GMTYPE_DOUBLE:
            file << ((double *) arr)[pos];
            break;
        case GMTYPE_NODE:
            file << ((node_t *) arr)[pos];
            break;
        case GMTYPE_EDGE:
            file << ((edge_t *) arr)[pos];
            break;
        
        case GMTYPE_END:
        default:
            assert(false);
            return; // Control should never reach this case.
    }
}

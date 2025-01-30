/*  Test functions for tabulating carmichaels project.
 *  Andrew Shallue and Jonathan Webster, January 2025
 */

#include "test_functions.h"

// call CN_factorization on several Preproduct objects
bool test_factor(){
    Preproduct PPP;
    uint64_t p = 3;
    uint64_t pl = 2;
    uint64_t bound = 1000;
    
    // initialize with P = 3, L = 2,  
    PPP.initializing(p, pl, bound);
    // now factorize 561 = 3 * 187
    mpz_t n;
    mpz_init(n);
    mpz_set_ui(n, 561);
    mpz_t r;
    mpz_init(r);
    mpz_set_ui(r, 187);

    PPP.CN_factorization(n, r);

    return true;
}

// main for testing
int main(){
    bool t1 = test_factor();

}

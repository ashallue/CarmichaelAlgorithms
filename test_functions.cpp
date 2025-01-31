/*  Test functions for tabulating carmichaels project.
 *  Andrew Shallue and Jonathan Webster, January 2025
 */

#include "test_functions.h"

// call CN_factorization on several Preproduct objects
bool test_factor(){
    // primes for building numbers to factor
    // 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97 
    uint64_t num_small_primes = 20;
    uint64_t* small_primes = new uint64_t[ num_small_primes ];
    small_primes[0] = 5;  small_primes[1] = 7; small_primes[2] = 11; small_primes[3] = 13; 
    small_primes[4] = 17;  small_primes[5] = 19; small_primes[6] = 23; small_primes[7] = 29; 
    small_primes[8] = 31;  small_primes[9] = 37; small_primes[10] = 41; small_primes[11] = 43; 
    small_primes[12] = 47;  small_primes[13] = 53; small_primes[14] = 59; small_primes[15] = 61;
    small_primes[16] = 67;  small_primes[17] = 71; small_primes[18] = 73; small_primes[19] = 79; 
    
    Preproduct PPP;
    uint64_t p = 3;
    uint64_t pl = 2;
    uint64_t bound = 1000;

    mpz_t n;
    mpz_init(n);
    mpz_t r;
    mpz_init(r);
    std::vector<uint64_t> r_primes;

    uint64_t factor_bound = 7;
        
    // we will factor multiples of 3 where we know the factors
    for(uint64_t i1 = 0; i1 < factor_bound; ++i1){
        // reset r, n and r_primes
        mpz_set_ui(r, 1);
        mpz_set_ui(n, p);
        r_primes.clear();

        for(uint64_t i2 = i1 + 1; i2 < factor_bound; ++i2){
            for(uint64_t i3 = i2 + 1; i3 < factor_bound; ++i3){
                PPP.initializing(p, pl, bound);

                
                // form r and n
                mpz_mul_ui( r, r, small_primes[i1] );
                mpz_mul_ui( r, r, small_primes[i2] );
                mpz_mul_ui( r, r, small_primes[i3] );
                mpz_mul(n, n, r);

                // factor
                bool is_carmichael = PPP.CN_factorization(n, r, r_primes);

                if( is_carmichael ){
                    gmp_printf ("%Zd = ", n );
                    std::cout << " is carmichael \n";
                }
            }
        }
    }

    delete[] small_primes;
    
    return true;
}

// main for testing
int main(){
    bool t1 = test_factor();

    

}

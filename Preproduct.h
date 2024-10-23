#ifndef PREPRODUCT_H
#define PREPRODUCT_H

#include <algorithm>
#include <cstdint>
#include <stdio.h>
#include <gmp.h>

// http://www.s369624816.websitehome.co.uk/rgep/cartable.html
// the least carmichael number with a fixed number of prime factors is known
// array lengths set at fixed length to accomodate this
// for B = 10^24 there is 1 CN with 14 prime factors
#define MAX_PRIME_FACTORS 14

// the intended use of this is with a precomputation that limits
// the total number of primes that can be appended in a prime-by-prime way
// for this computation, we choose that bound to be 5
#define APPEND_LIMIT 5

// contains a prime p and the factorization information for p-1 = Lambda(p)
struct primes_stuff
{
    uint32_t prime;
    uint32_t pm1_distinct_primes[10];
    uint16_t pm1_exponents[10];
    uint16_t pm1_len;
};

class Preproduct{
    
	
public:

	mpz_t P;
    uint64_t P_primes[ MAX_PRIME_FACTORS ];
    uint16_t P_len;
    uint64_t append_bound;  // primes appended to P need to exceed this bound

    // Information about L = CarmichaelLambda(P)
    uint64_t L;
    uint64_t L_distinct_primes[ MAX_PRIME_FACTORS ];
    uint16_t L_exponents[ MAX_PRIME_FACTORS ];   
    uint16_t L_len;

    uint16_t len_appended_primes;

    // these arrays are used to avoid gcd computations for admissibility checks
    // they are updated assuming primes are tested for admissibility in increasing order
    // These three arrays are the only data members that can change after initialization
    uint64_t next_inadmissible[ APPEND_LIMIT ]; 
    uint16_t mod_three_status[ APPEND_LIMIT ];  
    uint64_t appended_primes[ APPEND_LIMIT ];   

	// constructor and destructor
	Preproduct();
	~Preproduct();
    
	// initializing call
	// has to factor init_preproduct and init_LofP
	void initialization( uint64_t init_preproduct, uint64_t init_LofP, uint64_t init_append_bound );
	// appending call
	// assume we have an admissible prime to append. 
	// contains a merge computation of LCM( lambda(PP), p-1 )	
	void appending( Preproduct PP, primes_stuff p );

	// member functions
	// done with no gcd check
    bool is_admissible( uint64_t prime_to_append );


};

#endif

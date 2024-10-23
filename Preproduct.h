#ifndef PREPRODUCT_H
#define PREPRODUCT_H

#include <algorithm>
#include <cstdint>
#include <stdio.h>
#include <gmp.h>

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
    uint64_t P_primes[20];
    uint16_t P_len;
    uint64_t append_bound;  // primes appended to P need to exceed this bound

    // Information about L = CarmichaelLambda(P)
    uint64_t L;
    uint64_t L_distinct_primes[20];
    uint16_t L_exponents[20];   
    uint16_t L_len;

    // 5 is chosen because that's the limit of the recursion depth.
    uint16_t len_appended_primes;

    // these arrays are used to avoid gcd computations for admissibility checks
    // they are updated assuming primes are tested for admissibility in increasing order
    // These three arrays are the only data members that can change after initialization
    uint64_t next_inadmissible[5]; 
    uint16_t mod_three_status[5];  
    uint64_t appended_primes[5];   


	// initializing constructor
    Preproduct( uint64_t init_preproduct, uint64_t init_LofP, uint64_t init_append_bound );
	// appending constructor
	Preproduct( Preproduct PP, primes_stuff p );

	~Preproduct();
    
	// assume we have an admissible prime to append. 
	// contains a merge computation of LCM( lambda(PP), p-1 )	

	// member functions
    bool is_admissible( uint64_t prime_to_append );


};

#endif
